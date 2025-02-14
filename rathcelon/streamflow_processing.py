#Code to download GEOGLOWS retrospective datasets
#GEOGLOWS data can be downloaded from http://geoglows-v2.s3-website-us-west-2.amazonaws.com/

# third-party imports
import geopandas as gpd
import pandas as pd
from shapely.geometry import box
import s3fs
import xarray as xr


def Process_and_Write_Retrospective_Data_for_Dam(StrmShp_gdf, rivid_field, dam_csv, dam_id_field, dam_id, CSV_File_Name, OutShp_File_Name):

    # Load the dam data in as a geodataframe
    print('Process_and_Write_Retrospective_Data_for_Dam: Load the dam data in as a geodataframe')
    dam_gdf = pd.read_csv(dam_csv)
    dam_gdf = gpd.GeoDataFrame(dam_gdf, geometry=gpd.points_from_xy(dam_gdf['longitude'], dam_gdf['latitude']),
                               crs="EPSG:4269")

    # Filter the dam data to the dam of interest
    print('Process_and_Write_Retrospective_Data_for_Dam: Filter the dam data to the dam of interest')
    dam_gdf = dam_gdf[dam_gdf[dam_id_field] == dam_id]

    # Ensure there is at least one row remaining
    print('Process_and_Write_Retrospective_Data_for_Dam: Ensure there is at least one row remaining')
    if dam_gdf.empty:
        raise ValueError("No matching dam found for the given dam_id.")
    
    # Reset index to avoid index errors
    print('Process_and_Write_Retrospective_Data_for_Dam: Reset index to avoid index errors')
    dam_gdf = dam_gdf.reset_index(drop=True)

    # Ensure StrmShp_gdf has the same CRS as dam_gdf
    print('Process_and_Write_Retrospective_Data_for_Dam: Ensure StrmShp_gdf has the same CRS as dam_gdf')
    if StrmShp_gdf.crs != dam_gdf.crs:
        dam_gdf = dam_gdf.to_crs(StrmShp_gdf.crs)

    # Find the closest stream using distance calculation
    print('Process_and_Write_Retrospective_Data_for_Dam: Find the closest stream using distance calculation')
    StrmShp_gdf['distance'] = StrmShp_gdf.distance(dam_gdf.geometry.iloc[0])
    StrmShp_gdf = StrmShp_gdf.sort_values('distance')
    StrmShp_filtered_gdf = StrmShp_gdf.head(1)

    # # Use the 'LINKNO' and 'DSLINKNO' fields to find the stream upstream and downstream of the dam
    # current_rivid = StrmShp_filtered_gdf['LINKNO'].values[0]
    # downstream_rivid = StrmShp_filtered_gdf['DSLINKNO'].values[0]
    # upstream_StrmShp_gdf = StrmShp_gdf[StrmShp_gdf['DSLINKNO'] == current_rivid]
    # downstream_StrmShp_gdf = StrmShp_gdf[StrmShp_gdf['LINKNO'] == downstream_rivid]

    # Use the 'LINKNO' and 'DSLINKNO' fields to find the stream upstream and downstream of the dam
    print('Process_and_Write_Retrospective_Data_for_Dam: Use the LINKNO and DSLINKNO fields to find the stream upstream and downstream of the dam')
    current_rivid = StrmShp_filtered_gdf['LINKNO'].values[0]
    downstream_rivid = StrmShp_filtered_gdf['DSLINKNO'].values[0]

    # Find the upstream segment (if needed)
    print('Process_and_Write_Retrospective_Data_for_Dam: Find the upstream segment (if needed)')
    upstream_StrmShp_gdf = StrmShp_gdf[StrmShp_gdf['DSLINKNO'] == current_rivid]

    # Initialize a list to store the downstream segments.
    downstream_segments = []

    # Start with the dam's downstream segment.
    current_downstream_rivid = downstream_rivid

    # Loop to find up to 10 downstream segments.
    print('Process_and_Write_Retrospective_Data_for_Dam: Loop to find up to 10 downstream segments.')
    for i in range(11):
        # Find the stream segment whose LINKNO matches the current downstream rivid.
        segment = StrmShp_gdf[StrmShp_gdf['LINKNO'] == current_downstream_rivid]
        
        # If no segment is found, break the loop.
        if segment.empty:
            print(f"No downstream segment found after {i} segments.")
            break
        
        # Append the found segment to our list.
        downstream_segments.append(segment)
        
        # Update the current_downstream_rivid to the DSLINKNO of the found segment.
        # This will be used to find the next downstream segment.
        current_downstream_rivid = segment['DSLINKNO'].values[0]

    # Combine the downstream segments into one GeoDataFrame.
    print('Process_and_Write_Retrospective_Data_for_Dam: Combine the downstream segments into one GeoDataFrame.')
    if downstream_segments:
        downstream_StrmShp_gdf = pd.concat(downstream_segments, ignore_index=True)
    else:
        # If no downstream segments were found, create an empty GeoDataFrame.
        downstream_StrmShp_gdf = gpd.GeoDataFrame()

    # merge the StrmShp_filtered_gdf, upstream_StrmShp_gdf, and downstream_StrmShp_gdf into a single geodataframe
    print('Process_and_Write_Retrospective_Data_for_Dam: merge the StrmShp_filtered_gdf, upstream_StrmShp_gdf, and downstream_StrmShp_gdf into a single geodataframe')
    StrmShp_filtered_gdf = pd.concat([StrmShp_filtered_gdf, upstream_StrmShp_gdf, downstream_StrmShp_gdf])
    
    StrmShp_filtered_gdf.to_file(OutShp_File_Name)
    StrmShp_filtered_gdf[rivid_field] = StrmShp_filtered_gdf[rivid_field].astype(int)

    # create a list of river IDs to throw to AWS
    rivids_str = StrmShp_filtered_gdf[rivid_field].astype(str).to_list()
    rivids_int = StrmShp_filtered_gdf[rivid_field].astype(int).to_list()

    # Set up the S3 connection
    ODP_S3_BUCKET_REGION = 'us-west-2'
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name=ODP_S3_BUCKET_REGION))

    # Load FDC data from S3 using Dask
    # Convert to a list of integers
    print('Process_and_Write_Retrospective_Data_for_Dam: Load FDC data from S3 using Dask')
    fdc_s3_uri = 's3://geoglows-v2-retrospective/fdc.zarr'
    fdc_s3store = s3fs.S3Map(root=fdc_s3_uri, s3=s3, check=False)
    p_exceedance = [float(50.0), float(0.0)]
    fdc_ds = xr.open_zarr(fdc_s3store).sel(p_exceed=p_exceedance, river_id=rivids_str)
    # Convert Xarray to Dask DataFrame
    fdc_df = fdc_ds.to_dataframe().reset_index()

    # Check if fdc_df is empty
    print('Process_and_Write_Retrospective_Data_for_Dam: Check if fdc_df is empty')
    if fdc_df.empty:
        print(f"Skipping processing for {DEM_Tile} because fdc_df is empty.")
        CSV_File_Name = None
        OutShp_File_Name = None
        rivids_int = None
        StrmShp_filtered_gdf = None
        return (CSV_File_Name, OutShp_File_Name, rivids_int, StrmShp_filtered_gdf)

    print('Process_and_Write_Retrospective_Data_for_Dam: qout_median qout_max river_id')
    # Create 'qout_median' column where 'p_exceed' is 50.0
    fdc_df.loc[fdc_df['p_exceed'] == 50.0, 'qout_median'] = fdc_df['fdc']
    # Create 'qout_max' column where 'p_exceed' is 100.0
    fdc_df.loc[fdc_df['p_exceed'] == 0.0, 'qout_max'] = fdc_df['fdc']
    # Group by 'river_id' and aggregate 'qout_median' and 'qout_max' by taking the non-null value
    fdc_df = fdc_df.groupby('river_id').agg({
        'qout_median': 'max',  # or use 'max' as both approaches would work
        'qout_max': 'max'
    }).reset_index()

    # making our index for this dataframe match the recurrence interval index 
    fdc_df['rivid'] = fdc_df['river_id'].astype(int)
    # Drop two columns from the DataFrame
    fdc_df = fdc_df.drop(['river_id'], axis=1)
    fdc_df = fdc_df.set_index('rivid')

    # round the values
    fdc_df['qout_median'] = fdc_df['qout_median'].round(3)
    fdc_df['qout_max'] = fdc_df['qout_max'].round(3)
    
    # Load return periods data from S3 using Dask
    rp_s3_uri = 's3://geoglows-v2-retrospective/return-periods.zarr'
    rp_s3store = s3fs.S3Map(root=rp_s3_uri, s3=s3, check=False)
    rp_ds = xr.open_zarr(rp_s3store).sel(rivid=rivids_int)
    
    # Convert Xarray to Dask DataFrame and pivot
    rp_df = rp_ds.to_dataframe().reset_index()

    # Check if rp_df is empty
    print('Process_and_Write_Retrospective_Data_for_Dam: Check if rp_df is empty')
    if rp_df.empty:
        print(f"Skipping processing for {DEM_Tile} because rp_df is empty.")
        CSV_File_Name = None
        OutShp_File_Name = None
        rivids_int = None
        StrmShp_filtered_gdf = None
        return (CSV_File_Name, OutShp_File_Name, rivids_int, StrmShp_filtered_gdf)

    # Convert 'return_period' to category dtype
    rp_df['return_period'] = rp_df['return_period'].astype('category')
    
    # Pivot the table
    rp_pivot_df = rp_df.pivot_table(index='rivid', columns='return_period', values='return_period_flow', aggfunc='mean')

    # Rename columns to indicate return periods
    rp_pivot_df = rp_pivot_df.rename(columns={col: f'rp{int(col)}' for col in rp_pivot_df.columns})

    # Combine the results from retrospective and return periods data
    # final_df = pd.concat([combined_df, rp_pivot_df], axis=1)
    final_df = pd.concat([fdc_df, rp_pivot_df], axis=1)
    final_df['COMID'] = final_df.index

    # Column to move to the front
    target_column = 'COMID'

    # Reorder the DataFrame
    columns = [target_column] + [col for col in final_df.columns if col != target_column]
    final_df = final_df[columns]

    # Add a safety factor to one of the columns we could use to run the ARC model
    print('Process_and_Write_Retrospective_Data_for_Dam: Add a safety factor to one of the columns we could use to run the ARC model')
    for col in final_df.columns:
        if col in ['qout_max','rp100']:
            final_df[f'{col}_premium'] = round(final_df[col]*1.5, 3)
    
    print(final_df)

    # Write the final Dask DataFrame to CSV
    print('Process_and_Write_Retrospective_Data_for_Dam: Write the final Dask DataFrame to CSV')
    final_df.to_csv(CSV_File_Name, index=False)
    
    # Return the combined DataFrame as a Dask DataFrame
    return (CSV_File_Name, OutShp_File_Name, rivids_int, StrmShp_filtered_gdf)

if __name__ == "__main__":
    pass