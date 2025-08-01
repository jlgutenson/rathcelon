#Code to download GEOGLOWS retrospective datasets
#GEOGLOWS data can be downloaded from http://geoglows-v2.s3-website-us-west-2.amazonaws.com/

# third-party imports
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import s3fs
import xarray as xr
import requests
import io
import time

from .classes import Dam


def get_nwm_rp(comids: list[int]):
    rp_url = 'https://nwm-api.ciroh.org/return-period'

    header = {'x-api-key': 'AIzaSyC4BXXMQ9KIodnLnThFi5Iv4y1fDR4U1II'}
    params = {'comids': ','.join(map(str, comids)),
              'output_format': 'csv',
              'order_by_comid': False,}

    response = requests.get(rp_url, params=params, headers=header, timeout=60)

    if response.status_code == 200:
        return_period_df = pd.read_csv(io.StringIO(response.text))
    else:
        raise requests.exceptions.HTTPError(response.text)
    return_period_df = return_period_df.set_index("feature_id")
    return_period_df.index.name = "river_id"
    return_period_df.columns = ['rp2', 'rp5', 'rp10', 'rp25', 'rp50', 'rp100']

    return return_period_df


def get_stream_coords(stream_gdf: gpd.GeoDataFrame, rivid_field: str, rivids: list[int|str], method='centroid'):
    """

    Parameters
    ----------
    stream_gdf -
    rivid_field - either 'LINKNO' or 'hydroseq'
    rivids - list of all rivids...
    method - where we pull the coords from

    Returns
    -------
    results - each rivid and its coords {rivid: [lat, lon]}
    """
    # Ensure original CRS is EPSG:4326 for lat/lon output
    if stream_gdf.crs != "EPSG:4326":
        stream_gdf = stream_gdf.to_crs("EPSG:4326")

    result = {}

    for rivid in rivids:
        match = stream_gdf[stream_gdf[rivid_field] == rivid]

        if match.empty:
            print(f"No match found for {rivid_field} = {rivid}" )
            continue

        geom = match.geometry.iloc[0]

        # extract the point from the geom
        if method == "centroid":
            point = geom.centroid
        elif method == "start":
            point = Point(geom.coords[0])
        elif method == "end":
            point = Point(geom.coords[-1])
        else:
            raise ValueError(f"Unknown method '{method}'")

        result[rivid] = [point.y, point.x]

    return result


def Process_and_Write_Retrospective_Data_for_Dam(dam: Dam):
    """

    Parameters
    ----------
    dam: A Dam object... i'll write more later
    -------

    """
    # Load the dam data in as a geodataframe
    print('Process_and_Write_Retrospective_Data_for_Dam: Load the dam data in as a geodataframe')
    dam_gdf = pd.read_csv(dam.csv_path)
    dam_gdf = gpd.GeoDataFrame(dam_gdf, geometry=gpd.points_from_xy(dam_gdf['longitude'], dam_gdf['latitude']),
                               crs="EPSG:4269")

    # Filter the dam data to the dam of interest
    print('Process_and_Write_Retrospective_Data_for_Dam: Filter the dam data to the dam of interest')
    print(dam_gdf.tail())
    print(f'dam_id = {dam.dam_id}')
    dam_gdf = dam_gdf[dam_gdf[dam.id_field] == dam.dam_id]

    # Ensure there is at least one row remaining
    print('Process_and_Write_Retrospective_Data_for_Dam: Ensure there is at least one row remaining')
    if dam_gdf.empty:
        raise ValueError("No matching dam found for the given dam_id.")
    
    # Reset index to avoid index errors
    print('Process_and_Write_Retrospective_Data_for_Dam: Reset index to avoid index errors')
    dam_gdf = dam_gdf.reset_index(drop=True)

    # Print stage info
    print('Process_and_Write_Retrospective_Data_for_Dam: Convert both GeoDataFrames to a common projected CRS')

    # save the StrmShp CRS to convert StrmShp_filtered_gdf back to after the distance calculation
    flowline_crs = dam.flowline_gdf.crs

    # Determine an appropriate UTM zone using GeoPandas
    dam.flowline_gdf = dam.flowline_gdf[dam.flowline_gdf.geometry.notnull()]
    dam.flowline_gdf = dam.flowline_gdf[~dam.flowline_gdf.geometry.is_empty]


    if not dam.flowline_gdf.empty:
        utm_crs = dam.flowline_gdf.estimate_utm_crs()
    else:
        raise ValueError(f"{dam.flowline_gdf} has no valid geometries for UTM estimation.")

    # Reproject both GeoDataFrames to the UTM CRS
    dam_gdf = dam_gdf.to_crs(utm_crs)
    dam.flowline_gdf = dam.flowline_gdf.to_crs(utm_crs)

    # Distance calculation
    print('Process_and_Write_Retrospective_Data_for_Dam: Find the closest stream using distance calculation')
    dam.flowline_gdf['distance'] = dam.flowline_gdf.distance(dam_gdf.geometry.iloc[0])
    dam.flowline_gdf = dam.flowline_gdf.sort_values('distance')
    filtered_flowline_gdf = dam.flowline_gdf.head(1)

    # convert flowline_gdf and filtered_flowline_gdf back to its original CRS for ARC to use
    dam.flowline_gdf = dam.flowline_gdf.to_crs(flowline_crs)
    filtered_flowline_gdf = filtered_flowline_gdf.to_crs(flowline_crs)

    # # Use the 'LINKNO' and 'DSLINKNO' fields to find the stream upstream and downstream of the dam

    # current_rivid = StrmShp_filtered_gdf['LINKNO'].values[0]
    # downstream_rivid = StrmShp_filtered_gdf['DSLINKNO'].values[0]
    # upstream_gdf = dam.flowline_gdf[dam.flowline_gdf['DSLINKNO'] == current_rivid]
    # downstream_gdf = dam.flowline_gdf[dam.flowline_gdf['LINKNO'] == downstream_rivid]

    # set the field names equal to the values as they appear in GEOGLOWS or NHDPlus
    if dam.rivid_field == 'LINKNO':
        ds_rivid_field = 'DSLINKNO'
        stream_order_field = 'strmOrder'

    else: #  rivid_field == 'hydroseq'
        ds_rivid_field = 'dnhydroseq'
        stream_order_field = 'streamorde'

    # Use the 'LINKNO' and 'DSLINKNO' fields to find the stream upstream and downstream of the dam
    print('Process_and_Write_Retrospective_Data_for_Dam: Use the LINKNO and DSLINKNO (hydroseq and dnhydroseq) '
          'fields to find the stream upstream and downstream of the dam')
    current_rivid = filtered_flowline_gdf[dam.rivid_field].values[0]
    downstream_rivid = filtered_flowline_gdf[ds_rivid_field].values[0]

    # Find the upstream segment (if needed)
    print('Process_and_Write_Retrospective_Data_for_Dam: Find the upstream segment (if needed)')
    upstream_gdf = dam.flowline_gdf[dam.flowline_gdf[ds_rivid_field] == current_rivid]
    # Select the upstream segment with the highest Stream Order, if needed

    if not upstream_gdf.empty:
        if stream_order_field in upstream_gdf.columns:
            # Use the highest stream order if available
            upstream_gdf = upstream_gdf.loc[[upstream_gdf[stream_order_field].idxmax()]]
        else:
            # Use all matching upstream segments
            print(f"Optional field '{stream_order_field}' not found — using all upstream matches.")
    else:
        print("No upstream segments found.")

    # Initialize a list to store the downstream segments.
    downstream_segments = []

    # Start with the dam's downstream segment.
    current_downstream_rivid = downstream_rivid

    # Loop to find up to 2 downstream segments. lol
    print('Process_and_Write_Retrospective_Data_for_Dam: Loop to find up to 3 downstream segments.')
    for i in range(4):
        # Find the stream segment whose LINKNO matches the current downstream rivid.
        segment = dam.flowline_gdf[dam.flowline_gdf[dam.rivid_field] == current_downstream_rivid]
        
        # If no segment is found, break the loop.
        if segment.empty:
            print(f"No downstream segment found after {i} segments.")
            break
        
        # Append the found segment to our list.
        downstream_segments.append(segment)
        
        # Update the current_downstream_rivid to the DSLINKNO of the found segment.
        # This will be used to find the next downstream segment.
        current_downstream_rivid = segment[ds_rivid_field].values[0]

    # Combine the downstream segments into one GeoDataFrame.
    print('Process_and_Write_Retrospective_Data_for_Dam: Combine the downstream segments into one GeoDataFrame.')
    if downstream_segments:
        downstream_gdf = pd.concat(downstream_segments, ignore_index=True)
    else:
        # If no downstream segments were found, create an empty GeoDataFrame.
        downstream_gdf = gpd.GeoDataFrame()

    # merge the StrmShp_filtered_gdf, upstream_gdf, and downstream_gdf into a single geodataframe
    print('Process_and_Write_Retrospective_Data_for_Dam: merge the filtered_flowline_gdf, upstream_gdf, and downstream_gdf into a single geodataframe')
    StrmShp_filtered_gdf = pd.concat([filtered_flowline_gdf, upstream_gdf, downstream_gdf])
    
    StrmShp_filtered_gdf.to_file(dam.dam_shp)
    StrmShp_filtered_gdf[dam.rivid_field] = StrmShp_filtered_gdf[dam.rivid_field].astype(int)

    # create a list of river IDs to throw to AWS
    # rivids_str = StrmShp_filtered_gdf[rivid_field].astype(str).to_list()
    rivids_int = StrmShp_filtered_gdf[dam.rivid_field].astype(int).to_list()

    if dam.rivid_field == 'LINKNO':
        # Set up the S3 connection
        ODP_S3_BUCKET_REGION = 'us-west-2'
        s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name=ODP_S3_BUCKET_REGION))

        # # Load FDC data from S3 using Dask
        # # Convert to a list of integers
        fdc_s3_uri = 's3://geoglows-v2/retrospective/fdc.zarr'
        fdc_s3store = s3fs.S3Map(root=fdc_s3_uri, s3=s3, check=False)
        p_exceedance = [float(50.0), float(0.0)]
        fdc_ds = xr.open_zarr(fdc_s3store).sel(p_exceed=p_exceedance, river_id=rivids_int)


        # Convert Xarray to Dask DataFrame
        fdc_df = fdc_ds.to_dataframe().reset_index()


        # Check if fdc_df is empty
        if fdc_df.empty:
            # print(f"Skipping processing for {DEM_Tile} because fdc_df is empty.")
            dam.reanalysis_csv = None
            dam.dam_shp = None
            rivids_int = None
            StrmShp_filtered_gdf = None
            return dam.reanalysis_csv, dam.dam_shp, rivids_int, StrmShp_filtered_gdf

        # Create 'qout_median' column where 'p_exceed' is 50.0
        fdc_df.loc[fdc_df['p_exceed'] == 50.0, 'qout_median'] = fdc_df['hourly_annual']
        # Create 'qout_max' column where 'p_exceed' is 100.0
        fdc_df.loc[fdc_df['p_exceed'] == 0.0, 'qout_max'] = fdc_df['hourly_annual']
        # Group by 'river_id' and aggregate 'qout_median' and 'qout_max' by taking the non-null value
        fdc_df = fdc_df.groupby('river_id').agg({
            'qout_median': 'max',  # or use 'max' as both approaches would work
            'qout_max': 'max'
        }).reset_index()

        # set the dataframe index
        fdc_df = fdc_df.set_index('river_id')

        # round the values
        fdc_df['qout_median'] = fdc_df['qout_median'].round(3)
        fdc_df['qout_max'] = fdc_df['qout_max'].round(3)

        # drop all the columns except for the river_id, qout_median, and qout_max
        fdc_df = fdc_df[['qout_median', 'qout_max']]

        # Load return periods data from S3 using Dask
        rp_s3_uri = 's3://geoglows-v2/retrospective/return-periods.zarr'
        rp_s3store = s3fs.S3Map(root=rp_s3_uri, s3=s3, check=False)
        rp_ds = xr.open_zarr(rp_s3store).sel(river_id=rivids_int)

        # Convert Xarray to Dask DataFrame and pivot
        rp_df = rp_ds.to_dataframe().reset_index()

        # find the maximum between the gumbel and logpearson3 return periods and label this new column 'return_period_flow'
        rp_df['return_period_flow'] = rp_df[['gumbel', 'logpearson3']].mean(axis=1).round(3)

        # keep just the column 'return_period_flow'
        rp_df = rp_df[['river_id', 'return_period', 'return_period_flow']]

        # Check if rp_df is empty
        if rp_df.empty:
            # print(f"Skipping processing for {DEM_Tile} because rp_df is empty.")
            dam.reanalysis_csv = None
            dam.dam_shp = None
            rivids_int = None
            StrmShp_filtered_gdf = None
            return dam.reanalysis_csv, dam.dam_shp, rivids_int, StrmShp_filtered_gdf

        # Convert 'return_period' to category dtype
        rp_df['return_period'] = rp_df['return_period'].astype('category')

        # Pivot the table
        rp_pivot_df = rp_df.pivot_table(index='river_id', columns='return_period', values='return_period_flow', aggfunc='mean')

        # Rename columns to indicate return periods
        rp_pivot_df = rp_pivot_df.rename(columns={col: f'rp{int(col)}' for col in rp_pivot_df.columns})

        # Combine the results from retrospective and return periods data
        # final_df = pd.concat([combined_df, rp_pivot_df], axis=1)
        final_df = pd.concat([fdc_df, rp_pivot_df], axis=1)

    else:
        strm_coords = get_stream_coords(dam.flowline_gdf, dam.rivid_field, rivids_int)
        print(f'strm_coords: {strm_coords}')

        hydroseqs = []
        reach_ids = []

        for key, value in strm_coords.items():
            lat, lon = value
            url = f"https://nwm-api.ciroh.org/geometry?lat={lat}&lon={lon}&output_format=csv&key=AIzaSyC4BXXMQ9KIodnLnThFi5Iv4y1fDR4U1II"

            retries = 3
            for attempt in range(retries):
                try:
                    r = requests.get(url, timeout=30)
                    r.raise_for_status()

                    df = pd.read_csv(io.StringIO(r.text))
                    reach_id = df['station_id'].values[0]

                    hydroseqs.append(key)
                    reach_ids.append(reach_id)
                    break

                except (requests.exceptions.Timeout, requests.exceptions.HTTPError) as e:
                    print(f"Request failed for key={key} lat={lat} lon={lon}, attempt {attempt + 1}/{retries}")
                    print(f"Error: {e}")
                    if attempt < retries - 1:
                        wait = 2 ** attempt
                        print(f"Retrying in {wait} seconds...")
                        time.sleep(wait)
                    else:
                        print("Giving up on this point.\n")
                        continue

        # Fetch return periods (rp2, rp100, etc.)
        rp_df = get_nwm_rp(reach_ids)

        # Add derived flows directly to rp_df without dropping anything
        rp_df["qout_median"] = (rp_df["rp2"] / 2).round(3)
        rp_df["qout_max"] = (rp_df["rp100"] * 1.5).round(3)

        # Reorder columns so qout_median and qout_max come first
        cols = ["qout_median", "qout_max"] + [col for col in rp_df.columns if col.startswith("rp")]
        rp_df = rp_df[cols]
        print(rp_df)

        # Map hydroseqs to river IDs
        map_df = pd.DataFrame({'hydroseq': hydroseqs, 'river_id': reach_ids})
        print(map_df)

        # Merge and reshape
        final_df = map_df.merge(rp_df, on='river_id') \
            .set_index('hydroseq') \
            .drop(columns='river_id')

        # Set index name for clarity
        final_df.index.name = 'river_id'

        print(final_df)


    # else:
    #     strm_coords = get_stream_coords(dam.flowline_gdf, dam.rivid_field, rivids_int)
    #     print(f'strm_coords: {strm_coords}')
    #     # use the nwm api to get reach_id based on lat lon...
    #     hydroseqs = []
    #     reach_ids = []
    #
    #     for key, value in strm_coords.items():
    #         lat, lon = value
    #         url = f"https://nwm-api.ciroh.org/geometry?lat={lat}&lon={lon}&output_format=csv&key=AIzaSyC4BXXMQ9KIodnLnThFi5Iv4y1fDR4U1II"
    #
    #         retries = 3
    #         for attempt in range(retries):
    #             try:
    #                 r = requests.get(url, timeout=30)  # 30s max wait per request
    #                 r.raise_for_status()
    #
    #                 df = pd.read_csv(io.StringIO(r.text))
    #                 reach_id = df['station_id'].values[0]
    #
    #                 hydroseqs.append(key)
    #                 reach_ids.append(reach_id)
    #                 break  # success, exit retry loop
    #
    #             except (requests.exceptions.Timeout, requests.exceptions.HTTPError) as e:
    #                 print(f"Request failed for key={key} lat={lat} lon={lon}, attempt {attempt + 1}/{retries}")
    #                 print(f"Error: {e}")
    #
    #                 if attempt < retries - 1:
    #                     wait = 2 ** attempt
    #                     print(f"Retrying in {wait} seconds...")
    #                     time.sleep(wait)
    #                 else:
    #                     print("Giving up on this point.\n")
    #                     continue  # skip this point after final failure
    #
    #     # now hydroseq_reach_id has each hydroseq matched with its reach_id
    #
    #     # we don't really need the median or max flow... so this may be useful later, but i'll comment it out for now
    #     # if known_baseflow is None:
    #     s3_path = 's3://noaa-nwm-retrospective-3-0-pds/CONUS/zarr/chrtout.zarr'
    #     fs = s3fs.S3FileSystem(anon=True)
    #     nwm_ds = xr.open_zarr(fs.get_mapper(s3_path), consolidated=True)
    #     reach_ds = nwm_ds['streamflow'].sel(feature_id=reach_ids)
    #     df = reach_ds.to_dataframe().reset_index()
    #     df = df.groupby("feature_id").agg(
    #         qout_median=("streamflow", "median"),
    #         qout_max=("streamflow", "max"),
    #     ).round(3)
    #     df.index.name = 'river_id'
    #
    #     # else:
    #     #     df = pd.DataFrame({'river_id': reach_ids,
    #     #                        'qout_median': 1,
    #     #                        'qout_max': 100,}).set_index('river_id')
    #
    #     # Get return period DataFrame
    #     rp_df = get_nwm_rp(reach_ids)
    #
    #     # Merge it with df on index (river_id)
    #     final_df = pd.concat([df, rp_df], axis=1)
    #
    #     # Create a DataFrame to map hydroseq to river_id and join it in
    #     map_df = pd.DataFrame({'hydroseq': hydroseqs, 'river_id': reach_ids})
    #
    #     # Merge map_df with final_df on river_id, set hydroseq as index, and drop river_id
    #     final_df = map_df.merge(final_df, on='river_id') \
    #         .set_index('hydroseq') \
    #         .drop(columns='river_id')
    #
    #     # Rename the index to river_id
    #     final_df.index.name = 'river_id'
    #     print(final_df)

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

    # add the known_baseflow column to the final_df
    if dam.known_baseflow is not None:
        final_df['known_baseflow'] = dam.known_baseflow
    
    # add the known_channel_forming_discharge column to the final_df
    if dam.known_channel_forming_discharge is not None:
        final_df['known_channel_forming_discharge'] = dam.known_channel_forming_discharge

    # print(final_df)

    # Write the final Dask DataFrame to CSV
    print('Process_and_Write_Retrospective_Data_for_Dam: Write the final Dask DataFrame to CSV')
    final_df.to_csv(dam.reanalysis_csv, index=False)
    
    # Return the combined DataFrame as a Dask DataFrame
    return rivids_int, StrmShp_filtered_gdf
