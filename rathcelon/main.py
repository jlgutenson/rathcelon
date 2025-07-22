"""
Here is a list of where Low Head Dams are located.  What we're wanting to do:

1.) Identify the stream reach(es) that are associated with the Low-Head Dam.
2.) Determine the approximate Top-Width (TW) of the Dam (or the stream itself)
3.) Determine the average (base) flow for the stream as well as the seasonal high flow (not a flood, just a normal high-flow).
4.) Go downstream approximately TW distance and pull a perpendicular cross-section.
5.) Go downstream another TW distance (2*TW from the dam) and pull another perpendicular cross-section.
6.) Go downstream another TW distance (3*TW from the dam) and pull another perpendicular cross-section.
7.) For each cross-section estimate the bathymetry.
8.) For each cross-section calculate the rating curve of the cross-section.  Slope can likely be calculated from steps 4-6.

"""

# build-in imports
import os       # working with the operating system
import sys      # working with the interpreter
import argparse # parses command-line args
# import multiprocessing as mp

# third-party imports
from arc import Arc                 # automated rating curve generator (arc)
from pandas import DataFrame        # object for handling tabular data
from geopandas import GeoDataFrame  # object for geospatial data + tabular data
# import re
# import subprocess
# from datetime import datetime, timedelta

try:
    import gdal                         # geospatial data abstraction library (gdal)
except ImportError:
    from osgeo import gdal, ogr, osr    # import from osgeo if direct import doesn't work

import json         # json encoding/decoding
import stat         # interpreting file status
import fiona        # reading/writing vector spatial data
import platform     # access information about system
import rasterio     # raster data access/manipulation
import numpy as np      # numerical + py = numpy
import pandas as pd     # tabular data manipulation
import networkx as nx       # network/graph analysis
import geopandas as gpd     # geospatial + pandas = geopandas
from pathlib import Path    # OOP - programming + filesystem paths
from shapely.geometry import box, shape     # geometric shapes + manipulation
from pyproj import CRS, Geod, Transformer   # CRS and transformations
from shapely.geometry import Point, LineString #, MultiLineString, shape, mapping
from shapely.ops import nearest_points, linemerge, transform #, split
# import shutil
from rasterio.features import rasterize
# from rasterio.warp import calculate_default_transform, reproject, Resampling

# local imports
from . import streamflow_processing
from . import esa_download_processing as esa


def Process_Geospatial_Data(ARC_FileName_Bathy: str,  DEM_File: str, LandCoverFile: str,
                            STRM_File: str, STRM_File_Clean: str, LAND_File: str, ManningN: str, VDT_File: str, Curve_File: str,
                            ARC_BathyFile: str, Dam_StrmShp: str, Dam_Reanalsyis_FlowFile: str, XS_File_Out: str,
                            bathy_use_banks: bool, find_banks_based_on_landcover: bool,
                            create_reach_average_curve_file: bool, dam_csv: str, dam_id_field: str, dam_id: int|str,
                            known_baseflow: float, known_channel_forming_discharge: float,
                            StrmShp_gdf: gpd.GeoDataFrame):
    """

    Parameters
    ----------
    ARC_FileName_Bathy
    DEM_File
    LandCoverFile
    STRM_File
    STRM_File_Clean
    LAND_File
    ManningN: .txt file path w/ Manning's n based on LU
    VDT_File
    Curve_File
    ARC_BathyFile
    Dam_StrmShp
    Dam_Reanalsyis_FlowFile
    XS_File_Out
    bathy_use_banks
    find_banks_based_on_landcover
    create_reach_average_curve_file
    dam_csv: .csv file w/ all dam info
    dam_id_field: name of ID field, e.g., 'ID'
    dam_id: ...
    known_baseflow: baseflow for provided DEM
    known_channel_forming_discharge: discharge at bank-full depth
    StrmShp_gdf

    Returns
    -------

    """

    #Get the Spatial Information from the DEM Raster
    (minx, miny, maxx, maxy, dx, dy, ncols, nrows, dem_geoTransform, dem_projection) = Get_Raster_Details(DEM_File)
    projWin_extents = [minx, maxy, maxx, miny]
    # outputBounds = [minx, miny, maxx, maxy]  #https://gdal.org/api/python/osgeo.gdal.html
    rivid_field = None

    #Create Land Dataset
    if os.path.isfile(LAND_File):
        print(f'{LAND_File} Already Exists')
    else: 
        print(f'Creating {LAND_File}')
        # Let's make sure all the GIS data is using the same coordinate system as the DEM
        LandCoverFile = Check_and_Change_Coordinate_Systems(DEM_File, LandCoverFile)
        Create_AR_LandRaster(LandCoverFile, LAND_File, projWin_extents, ncols, nrows)

    # now we need to figure out if our Dam_StrmShp and DEM_Reanalysis_Flowfile exists and if not, create it
    if os.path.isfile(Dam_StrmShp) and os.path.isfile(Dam_Reanalsyis_FlowFile):
        print(f'{Dam_StrmShp} Already Exists\n'
              f'{Dam_Reanalsyis_FlowFile} Already Exists')
        Dam_StrmShp_gdf = gpd.read_file(Dam_StrmShp)

        if 'LINKNO' in Dam_StrmShp_gdf.columns:
            rivid_field = 'LINKNO'
        else:
            rivid_field = 'hydroseq'

        rivids = Dam_StrmShp_gdf[rivid_field].values
        print(rivids)

    elif StrmShp_gdf is not None and os.path.isfile(Dam_StrmShp) is False and os.path.isfile(Dam_Reanalsyis_FlowFile) is False:
        if 'LINKNO' in StrmShp_gdf.columns:
            rivid_field = 'LINKNO'
        else:
            rivid_field = 'hydroseq'
        print('Running Function: Process_and_Write_Retrospective_Data_for_Dam')
        (DEM_Reanalsyis_FlowFile, Dam_StrmShp, rivids, Dam_StrmShp_gdf) = streamflow_processing.Process_and_Write_Retrospective_Data_for_Dam(StrmShp_gdf, rivid_field, dam_csv,
                                                                                                                                 dam_id_field, dam_id, known_baseflow,
                                                                                                                                 known_channel_forming_discharge,
                                                                                                                                 Dam_Reanalsyis_FlowFile,
                                                                                                                                 Dam_StrmShp)

    #Create Stream Raster
    if os.path.isfile(STRM_File):
        print(f'{STRM_File} Already Exists')
    else:
        print(f'Creating {STRM_File}\n'
              f'\tby rasterizing {Dam_StrmShp}')
        print(f'rivid_field: {rivid_field}')
        Create_AR_StrmRaster(Dam_StrmShp, STRM_File, DEM_File, rivid_field)
    
    #Clean Stream Raster
    if os.path.isfile(STRM_File_Clean):
        print(f'{STRM_File_Clean} Already Exists')
    else:
        print(f'Creating {STRM_File_Clean}')
        Clean_STRM_Raster(STRM_File, STRM_File_Clean)

    #Get the unique values for all the stream ids
    (S, ncols, nrows, cellsize, yll, yur, xll, xur, lat, dem_geotransform, dem_projection) = Read_Raster_GDAL(STRM_File_Clean)
    (RR,CC) = S.nonzero()
    num_strm_cells = len(RR)
    print(f'num_strm_cells: {num_strm_cells}')
    COMID_Unique = np.unique(S)
    # COMID_Unique = np.delete(COMID_Unique, 0)  #We don't need the first entry of zero
    COMID_Unique = COMID_Unique[np.where(COMID_Unique > 0)]
    COMID_Unique = np.sort(COMID_Unique).astype(int)
    num_comids = len(COMID_Unique)
    print(f'num_comids: {num_comids}')


    #Create the Bathy Input File

    # Let's extract the weir length real quick
    dam_gdf = pd.read_csv(dam_csv)
    # Filter to the specific dam
    dam_gdf = dam_gdf[dam_gdf[dam_id_field] == dam_id]
    weir_length = dam_gdf['weir_length'].values[0]

    print(f'Creating ARC Input File: {ARC_FileName_Bathy}')

    if bathy_use_banks is False and known_baseflow is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'qout_median'
        Q_Param = 'rp2'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True and known_baseflow is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'rp2'
        Q_Param = 'rp2'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is False and known_baseflow is not None and known_channel_forming_discharge is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_baseflow'
        Q_Param = 'rp2'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True and known_channel_forming_discharge is not None and known_baseflow is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_channel_forming_discharge'
        Q_Param = 'rp2'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is False and known_baseflow is not None and known_channel_forming_discharge is not None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_baseflow'
        Q_Param = 'rp2'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True and known_channel_forming_discharge is not None and known_baseflow is not None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_channel_forming_discharge'
        Q_Param = 'rp2'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    else:
        print('Error: bathy_use_banks and known_baseflow are both set to True.  Please check your inputs.\n')
        print('You want to pair known_baseflow with bathy_use_banks set to False.')
        sys.exit("Terminating: Invalid input combination.")

    
    return ARC_FileName_Bathy, Dam_Reanalsyis_FlowFile, Dam_StrmShp


def Create_FlowFile(MainFlowFile: str, FlowFileName: str, OutputID: int|str, Qparam: str) -> None:
    with open(MainFlowFile, 'r') as infile:
        lines = infile.readlines()

    headers = lines[0].strip().split(',')

    try:
        q_index = headers.index(Qparam)
        id_index = headers.index(OutputID)
    except ValueError as e:
        raise ValueError(f"Missing excepted column in input file: {e}")

    with open(FlowFileName, 'w') as outfile:
        outfile.write(f"{OutputID},{Qparam}")
        for line in lines[1:]:
            values = line.strip().split(',')
            outfile.write(f"{values[id_index]},{values[q_index]}\n")


def Create_Folder(F: str) -> None:
    if not os.path.exists(F): 
        os.makedirs(F)


def Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy: str, DEM_File: str, COMID_Param: str, Q_BF_Param: str,
                                      Q_Param: str, STRM_File_Clean: str, LAND_File: str, DEM_Reanalsyis_FlowFile: str,
                                      VDT_File: str, Curve_File: str, XS_Out_File: str, ManningN: str,
                                      ARC_BathyFile: str, DEM_StrmShp: str, bathy_use_banks: bool, weir_length: float,
                                      find_banks_based_on_landcover: bool, create_reach_average_curve_file: bool)\
        -> None:

    x_section_dist = int(10 * weir_length)

    with open(ARC_FileName_Bathy,'w') as out_file:
        out_file.write('#ARC_Inputs\n')
        out_file.write(f'DEM_File\t{DEM_File}\n')
        out_file.write(f'Stream_File\t{STRM_File_Clean}\n')
        out_file.write(f'LU_Raster_SameRes\t{LAND_File}\n')
        out_file.write(f'LU_Manning_n\t{ManningN}\n')
        out_file.write(f'Flow_File\t{DEM_Reanalsyis_FlowFile}\n')
        out_file.write(f'Flow_File_ID\t{COMID_Param}\n')
        out_file.write(f'Flow_File_BF\t{Q_BF_Param}\n')
        out_file.write(f'Flow_File_QMax\t{Q_Param}\n')
        out_file.write(f'Spatial_Units\tdeg\n')
        out_file.write(f'X_Section_Dist\t{x_section_dist}\n')
        out_file.write(f'Degree_Manip\t6.1\n')
        out_file.write(f'Degree_Interval\t1.5\n')
        out_file.write(f'Low_Spot_Range\t2\n')
        out_file.write(f'Str_Limit_Val\t1\n')
        out_file.write(f'Gen_Dir_Dist\t10\n')
        out_file.write(f'Gen_Slope_Dist\t10\n\n')

        out_file.write('#VDT_Output_File_and_CurveFile\n')
        out_file.write('VDT_Database_NumIterations\t30\n')
        out_file.write(f'VDT_Database_File\t{VDT_File}\n')
        out_file.write(f'Print_VDT_Database\t{VDT_File}\n')
        out_file.write(f'Print_Curve_File\t{Curve_File}\n')
        out_file.write(f'Reach_Average_Curve_File\t{create_reach_average_curve_file}\n\n')

        out_file.write('#Bathymetry_Information\n')
        out_file.write('Bathy_Trap_H\t0.20\n')
        out_file.write(f'Bathy_Use_Banks\t{bathy_use_banks}\n')
        if find_banks_based_on_landcover is True:
            out_file.write(f'FindBanksBasedOnLandCover\t{find_banks_based_on_landcover}\n')
        out_file.write(f'AROutBATHY\t{ARC_BathyFile}\n')
        out_file.write(f'BATHY_Out_File\t{ARC_BathyFile}\n')
        out_file.write(f'XS_Out_File\t{XS_Out_File}\n')
    

def Create_BaseLine_Manning_n_File(ManningN: str):
    with open(ManningN, 'w') as out_file:
        out_file.write('LC_ID\tDescription\tManning_n\n')
        out_file.write('11\tWater\t0.030\n')
        out_file.write('21\tDev_Open_Space\t0.013\n')
        out_file.write('22\tDev_Low_Intensity\t0.050\n')
        out_file.write('23\tDev_Med_Intensity\t0.075\n')
        out_file.write('24\tDev_High_Intensity\t0.100\n')
        out_file.write('31\tBarren_Land\t0.030\n')
        out_file.write('41\tDecid_Forest\t0.120\n')
        out_file.write('42\tEvergreen_Forest\t0.120\n')
        out_file.write('43\tMixed_Forest\t0.120\n')
        out_file.write('52\tShrub\t0.050\n')
        out_file.write('71\tGrass_Herb\t0.030\n')
        out_file.write('81\tPasture_Hay\t0.040\n')
        out_file.write('82\tCultivated_Crops\t0.035\n')
        out_file.write('90\tWoody_Wetlands\t0.100\n')
        out_file.write('95\tEmergent_Herb_Wet\t0.100')


def Create_BaseLine_Manning_n_File_ESA(ManningN: str):
    with open(ManningN, 'w') as out_file:
        out_file.write('LC_ID\tDescription\tManning_n\n')
        out_file.write('10\tTree Cover\t0.120\n')
        out_file.write('20\tShrubland\t0.050\n')
        out_file.write('30\tGrassland\t0.030\n')
        out_file.write('40\tCropland\t0.035\n')
        out_file.write('50\tBuiltup\t0.075\n')
        out_file.write('60\tBare\t0.030\n')
        out_file.write('70\tSnowIce\t0.030\n')
        out_file.write('80\tWater\t0.030\n')
        out_file.write('90\tEmergent_Herb_Wet\t0.100\n')
        out_file.write('95\tMangroves\t0.100\n')
        out_file.write('100\tMossLichen\t0.100')


def Create_AR_LandRaster(LandCoverFile: str, LAND_File: str, projWin_extents, ncols, nrows):
    ds = gdal.Open(LandCoverFile)
    ds = gdal.Translate(LAND_File, ds, projWin = projWin_extents, width=ncols, height = nrows)
    del ds


# def Create_AR_StrmRaster(StrmSHP: str, STRM_File: str, outputBounds, ncols, nrows, Param):
#     print(StrmSHP)
#     source_ds = gdal.OpenEx(StrmSHP)
#
#     gdal.Rasterize(STRM_File, source_ds, format='GTiff', outputType=gdal.GDT_Int32, outputBounds = outputBounds, width = ncols, height = nrows, noData = -9999, attribute = Param)
#     del source_ds



def Create_AR_StrmRaster(StrmSHP: str, output_raster_path: str, DEM_File: str, value_field: str):
    # Load the shapefile
    gdf = gpd.read_file(StrmSHP)

    # Load the reference raster to get transform, shape, and CRS
    with rasterio.open(DEM_File) as ref:
        meta = ref.meta.copy()
        ref_transform = ref.transform
        out_shape = (ref.height, ref.width)
        crs = ref.crs

    # Reproject the shapefile to match reference CRS
    if gdf.crs != crs:
        gdf = gdf.to_crs(crs)

    # Prepare shapes (geometry, value)
    shapes = [(geom, value if value is not None else 0)
              for geom, value in zip(gdf.geometry, gdf[value_field])]

    # Rasterize
    raster = rasterize(shapes=shapes, out_shape=out_shape, fill=0, transform=ref_transform, dtype='int64')

    # Update metadata
    meta.update({"driver": "GTiff", "dtype": "int64", "count": 1, "compress": "lzw"})

    # Write to file
    with rasterio.open(output_raster_path, 'w', **meta) as dst:
        dst.write(raster, 1)

    print(f"Raster written to {output_raster_path}")


def Write_Output_Raster(s_output_filename: str, raster_data, ncols: int, nrows: int, dem_geotransform, dem_projection, s_file_format, s_output_type):
    o_driver = gdal.GetDriverByName(s_file_format)  #Typically will be a GeoTIFF "GTiff"
    
    # Construct the file with the appropriate data shape
    o_output_file = o_driver.Create(s_output_filename, xsize=ncols, ysize=nrows, bands=1, eType=s_output_type)

    # Get the first band (assuming a single-band raster)
    band = o_output_file.GetRasterBand(1)

    # Initialize the band with zeros
    band.Fill(0)

    # Write to disk
    band.FlushCache()
    
    # Set the geotransform
    o_output_file.SetGeoTransform(dem_geotransform)
    
    # Set the spatial reference
    o_output_file.SetProjection(dem_projection)
    
    # Write the data to the file
    o_output_file.GetRasterBand(1).WriteArray(raster_data)
    
    # Once we're done, close properly the dataset
    del o_output_file


def Get_Raster_Details(DEM_File: str):
    print(DEM_File)
    with rasterio.open(DEM_File) as dataset:
        geoTransform = dataset.transform # affine transform
        ncols = dataset.width
        nrows = dataset.height
        minx = geoTransform.c
        dx = geoTransform.a
        maxy = geoTransform.f
        dy = geoTransform.e
        maxx = minx + dx * ncols
        miny = maxy + dy * nrows
        Rast_Projection = dataset.crs.to_wkt() # wkt format

    return minx, miny, maxx, maxy, dx, dy, ncols, nrows, geoTransform, Rast_Projection


def Read_Raster_GDAL(InRAST_Name: str):
    try:
        dataset = gdal.Open(InRAST_Name, gdal.GA_ReadOnly)
        if dataset is None:
            raise RuntimeError("Field Raster File cannot be read?")
    except RuntimeError as e:
        sys.exit(f"ERROR: {e}")

    # Retrieve dimensions of cell size and cell count then close DEM dataset
    geotransform = dataset.GetGeoTransform()
    band = dataset.GetRasterBand(1)
    RastArray = band.ReadAsArray()

    #global ncols, nrows, cellsize, yll, yur, xll, xur
    ncols=dataset.RasterXSize
    nrows=dataset.RasterYSize

    cellsize = geotransform[1]
    yll = geotransform[3] - nrows * abs(geotransform[5])
    yur = geotransform[3]
    xll = geotransform[0]
    xur = xll + ncols * geotransform[1]
    lat = abs((yll + yur) / 2.0)

    Rast_Projection = dataset.GetProjectionRef()

    # close datasets
    del dataset
    del band

    print(
        f'Spatial Data for Raster File:\n'
        f'\tncols = {ncols}\n'
        f'\tnrows = {nrows}\n'
        f'\tcellsize = {cellsize}\n'
        f'\tyll = {yll}\n'
        f'\tyur = {yur}\n'
        f'\txll = {xll}\n'
        f'\txur = {xur}'
    )

    return RastArray, ncols, nrows, cellsize, yll, yur, xll, xur, lat, geotransform, Rast_Projection


def Clean_STRM_Raster(STRM_File: str, STRM_File_Clean: str) -> None:
    print('\nCleaning up the Stream File.')
    (SN, ncols, nrows, cellsize, yll, yur, xll, xur, lat, dem_geotransform, dem_projection) = Read_Raster_GDAL(STRM_File)
    
    #Create an array that is slightly larger than the STRM Raster Array
    B = np.zeros((nrows + 2, ncols + 2))
    
    #Imbed the STRM Raster within the Larger Zero Array
    B[1:(nrows + 1), 1:(ncols + 1)] = SN
    
    #Added this because sometimes the non-stream values end up as -9999
    B = np.where(B > 0, B, 0)
    #(RR,CC) = B.nonzero()
    (RR, CC) = np.where(B > 0)
    num_nonzero = len(RR)
    
    for filterpass in range(2):
        #First pass is just to get rid of single cells hanging out not doing anything
        p_count = 0
        p_percent = (num_nonzero + 1) / 100.0
        n = 0
        for x in range(num_nonzero):
            if x >= p_count * p_percent:
                p_count = p_count + 1
                print(f' {p_count}', end =" ")
            r = RR[x]
            c = CC[x]
            V = B[r, c]
            if V > 0:
                #Left and Right cells are zeros
                if B[r, c + 1] == 0 and B[r, c - 1] == 0:
                    #The bottom cells are all zeros as well, but there is a cell directly above that is legit
                    if (B[r + 1, c - 1] + B[r + 1, c] + B[r + 1, c + 1]) == 0 and B[r - 1, c] > 0:
                        B[r, c] = 0
                        n += 1
                    #The top cells are all zeros as well, but there is a cell directly below that is legit
                    elif (B[r - 1, c - 1] + B[r - 1, c] + B[r - 1, c + 1]) == 0 and B[r + 1, c] > 0:
                        B[r, c] = 0
                        n += 1
                #top and bottom cells are zeros
                if B[r, c] > 0 and B[r + 1, c] == 0 and B[r - 1, c] == 0:
                    #All cells on the right are zero, but there is a cell to the left that is legit
                    if (B[r + 1, c + 1] + B[r, c + 1] + B[r - 1, c + 1]) == 0 and B[r, c- 1] > 0:
                        B[r, c] = 0
                        n += 1
                    elif (B[r + 1, c - 1] + B[r, c - 1] + B[r - 1, c - 1]) == 0 and B[r, c + 1] > 0:
                        B[r, c] = 0
                        n += 1
        print(f'\nFirst pass removed {n} cells')
        
        
        #This pass is to remove all the redundant cells
        n = 0
        p_count = 0
        p_percent = (num_nonzero + 1) / 100.0
        for x in range(num_nonzero):
            if x >= p_count * p_percent:
                p_count = p_count + 1
                print(f' {p_count}', end =" ")
            r = RR[x]
            c = CC[x]
            V = B[r, c]
            if V > 0:
                if B[r + 1, c] == V and (B[r + 1, c + 1] == V or B[r + 1, c - 1] == V):
                    if sum(B[r + 1, c - 1:c + 2]) == 0:
                        B[r + 1, c] = 0
                        n += 1
                elif B[r - 1, c] == V and (B[r - 1, c + 1] == V or B[r - 1, c - 1] == V):
                    if sum(B[r - 1,c - 1:c + 2]) == 0:
                        B[r - 1, c] = 0
                        n += 1
                elif B[r, c + 1] == V and (B[r + 1, c + 1] == V or B[r - 1, c + 1] == V):
                    if sum(B[r - 1:r + 1, c + 2]) == 0:
                        B[r, c + 1] = 0
                        n += 1
                elif B[r, c - 1] == V and (B[r + 1, c - 1] == V or B[r - 1, c - 1] == V):
                    if sum(B[r - 1:r + 1, c - 2]) == 0:
                            B[r, c - 1] = 0
                            n += 1

        print(f'\nSecond pass removed {n} redundant cells')
    
    print(f'Writing Output File {STRM_File_Clean}')
    Write_Output_Raster(STRM_File_Clean, B[1:nrows+1,1:ncols+1], ncols, nrows, dem_geotransform, dem_projection, "GTiff", gdal.GDT_Int64)
    #return B[1:nrows+1,1:ncols+1], ncols, nrows, cellsize, yll, yur, xll, xur


def Check_and_Change_Coordinate_Systems(DEM_File: str, LandCoverFile: str) -> str:
    # Load the projection of the DEM file
    with rasterio.open(DEM_File) as src:
        dem_projection = src.crs
    src.close()

    # re-project the LAND file raster, if necessary
    with rasterio.open(LandCoverFile) as src:
        current_crs = src.crs
    src.close()
        
    if current_crs != dem_projection:
        input_raster = gdal.Open(LandCoverFile)
        LandCoverFile_Update = f"{LandCoverFile[:-4]}_new.tif"
        output_raster = LandCoverFile_Update
        gdal.Warp(output_raster, input_raster, dstSRS=dem_projection)
        # Closes the files
        del input_raster

        # delete the old LAND raster, if it was replaced and change the name
        os.remove(LandCoverFile)
        LandCoverFile = LandCoverFile_Update

    return LandCoverFile


def find_nearest_idx(point, tree, gdf: GeoDataFrame):
    """Find the nearest index and corresponding data for a given point using spatial indexing."""
    nearest_idx = tree.nearest(point)  # Directly query with the geometry
    return nearest_idx, gdf.iloc[nearest_idx]


# Function to get a point strictly on the network
def get_point_on_stream(line, target_distance) -> Point:
    """
    Walks along a LineString and finds a point at a given distance,
    ensuring it follows the stream path without deviating.
    """
    current_distance = 0  # Start distance tracking

    for i in range(len(line.coords) - 1):
        start_pt = Point(line.coords[i])
        end_pt = Point(line.coords[i + 1])

        segment = LineString([start_pt, end_pt])
        segment_length = segment.length

        # If target_distance falls within this segment
        if current_distance + segment_length >= target_distance:
            # Compute exact location on this segment
            remaining_distance = target_distance - current_distance
            return segment.interpolate(remaining_distance)

        current_distance += segment_length  # Update walked distance

    return Point(line.coords[-1])  # If we exceed the length, return last point


def walk_stream_for_point(line, target_distance: float) -> Point:
    """
        Walks along a LineString segment-by-segment, stopping at the exact
        cumulative distance to ensure the point stays on the stream network.
    """
    current_distance = 0

    for i in range(len(line.coords) - 1):
        start_pt = Point(line.coords[i])
        end_pt = Point(line.coords[i + 1])
        segment = LineString([start_pt, end_pt])
        segment_length = segment.length

        # If our target point is within this segment
        if current_distance + segment_length >= target_distance:
            remaining_distance = target_distance - current_distance
            return segment.interpolate(remaining_distance)  # Stay exactly on the stream

        current_distance += segment_length  # Move forward

    return Point(line.coords[-1])  # Return last point if over length


def move_upstream(point, current_link, distance, G):
    # same logic as your original while-block, but for a fixed 'distance'
    remaining = distance
    curr_pt = point
    link = current_link

    while remaining > 0:
        # retrieve segment geometry for this link
        if 'geometry' in G.nodes[link]:
            seg_geom = G.nodes[link]['geometry']
        else:
            # fallback to edge data
            for u, v, data in list(G.out_edges(link, data=True)) + list(G.in_edges(link, data=True)):
                if 'geometry' in data:
                    seg_geom = data['geometry']
                    break
            else:
                raise RuntimeError(f"No geometry found for link {link}")

        # flatten MultiLine if needed
        if hasattr(seg_geom, 'geom_type') and seg_geom.geom_type.startswith("Multi"):
            seg_line = linemerge(seg_geom)
        else:
            seg_line = seg_geom

        # ensure seg_line endpoints with current point at end
        coords = list(seg_line.coords)
        if not Point(coords[-1]).equals_exact(curr_pt, tolerance=1e-6):
            coords.reverse()
            seg_line = LineString(coords)

        # how far along this segment we are
        proj = seg_line.project(curr_pt)
        available = proj

        if available >= remaining:
            # can step within this segment
            new_pt = seg_line.interpolate(proj - remaining)
            remaining = 0
        else:
            # go to upstream end of this segment
            remaining -= available
            new_pt = Point(seg_line.coords[0])
            # step to upstream link
            ups = list(G.in_edges(link))
            if not ups:
                raise RuntimeError("Ran out of upstream network")
            # pick first, or implement smarter selection
            link = ups[0][0]
        curr_pt = new_pt
    return curr_pt, link


def find_stream_cells_at_increments_above_and_below_dam(CurveParam_File: str, VDT_File: str, XS_Out_File: str,
                                                        dam_csv: str, dam_id_field: str, dam_id: str|int,
                                                        Dam_StrmShp: str, dam_reanalysis_flowfile: str, DEM_File: str,
                                                        upstream_elevation_change_threshold: float,
                                                        distance_upstream_increment: float, known_baseflow: float,
                                                        Rast_Projection: str, number_of_cross_sections=4)\
        -> tuple[GeoDataFrame, DataFrame, DataFrame, GeoDataFrame]:
    """
    Finds a location on the stream network that is downstream of the dam at
        specified increments and saves it as a shapefile.

    Parameters:
    - CurveParam_File: str, path to CurveParam file output by ARC
    - VDT_File: str, path to VDT file output by ARC
    - XS_Out_File: str, path to output cross-section file output by ARC
    - dam_csv: str, path to the CSV containing dam locations
    - dam_id_field: str, name of the column identifying dams
    - dam_id: int/str, ID of the dam of interest
    - Dam_StrmShp: str, path to the stream shapefile
    - dam_reanalysis_flowfile: str, path to the flow reanalysis file
    - output_shapefile: str, path to save the output shapefile
    - STRM_Raster_File: str, path to the stream raster file
    - Rast_Projection: str, original projection of the DEM raster
    - upstream_elevation_change_threshold: float, elevation change threshold to use when
        identifying upstream cross-section
    - distance_upstream_increment: float, increment distance upstream from the upstream dam cross-section
    - known_baseflow: float, known baseflow value to use for top width calculation 
    - number_of_cross_sections: int, number of distances at which to extract points (optional)

    Returns:
    - downstream_points: List of shapely.geometry.Point, the calculated downstream locations.
    """

    # Read the VDT and Curve data into DataFrames
    vdt_df = pd.read_csv(VDT_File)
    curve_data_df = pd.read_csv(CurveParam_File)
    
    # Read the dam reanalysis flow file
    dam_reanalysis_df = pd.read_csv(dam_reanalysis_flowfile)

    # Merge flow parameters
    merged_df = pd.merge(curve_data_df, dam_reanalysis_df, on='COMID', how='left')
    if known_baseflow is None:
        merged_df['tw_rp2'] = merged_df['tw_a'] * (merged_df['rp2'] ** merged_df['tw_b'])
    else:
        merged_df['tw_known_baseflow'] = merged_df['tw_a'] * (known_baseflow ** merged_df['tw_b'])


    # Read stream shapefile and dam locations
    Dam_StrmShp_gdf = gpd.read_file(Dam_StrmShp, engine='pyogrio')
    dam_gdf = pd.read_csv(dam_csv)
    dam_gdf = gpd.GeoDataFrame(dam_gdf, geometry=gpd.points_from_xy(dam_gdf['longitude'], dam_gdf['latitude']),
                               crs="EPSG:4269")

    # Convert to a projected CRS for accurate distance calculations
    projected_crs = Dam_StrmShp_gdf.estimate_utm_crs()
    Dam_StrmShp_gdf = Dam_StrmShp_gdf.to_crs(projected_crs)
    dam_gdf = dam_gdf.to_crs(projected_crs)

    # Filter to the specific dam
    dam_gdf = dam_gdf[dam_gdf[dam_id_field] == dam_id]
    if dam_gdf.empty:
        raise ValueError("No matching dam found for the given dam_id.")

    dam_gdf = dam_gdf.reset_index(drop=True)
    dam_point = dam_gdf.geometry.iloc[0]

    # **1. Build a Directed Graph Using LINKNO and DSLINKNO**
    G = nx.DiGraph()

    if 'LINKNO' in Dam_StrmShp_gdf.columns:
        rivid_field = 'LINKNO'
        ds_rivid_field = 'DSLINKNO'
    else:
        rivid_field = 'hydroseq'
        ds_rivid_field = 'dnhydroseq'


    
    for _, row in Dam_StrmShp_gdf.iterrows():
        link_id = row[rivid_field]
        ds_link_id = row[ds_rivid_field]
        geometry = row.geometry

        if ds_link_id > 0:  # Ignore terminal reaches
            G.add_edge(link_id, ds_link_id, geometry=geometry, weight=geometry.length)

    # Find the Closest Stream to the Dam**
    Dam_StrmShp_gdf['distance'] = Dam_StrmShp_gdf.distance(dam_point)
    closest_stream = Dam_StrmShp_gdf.loc[Dam_StrmShp_gdf['distance'].idxmin()]
    start_link = closest_stream[rivid_field]

    """
        none of the following variables were used later, so i just commented them out... i'll see if it still runs...
    """
    # Filter the top-width data to the closest stream and find the median
    # merged_df = merged_df[merged_df['COMID']==start_link]
    # if known_baseflow is None:
    #     tw_median = merged_df.groupby('COMID')['tw_rp2'].median()
    # else:
    #     tw_median = merged_df.groupby('COMID')['tw_known_baseflow'].median()
    

    # Calculate the Top Width (tw) for the Closest Stream and find the locations of interest on the stream network**
    # Select tw_median based on closest stream 'COMID'
    # tw = tw_median.get(start_link, 100)  # Default to 50m if not found

    # use a minimum of 100 meters for the top width
    # if tw < 100:
    #     tw = 100


    # Get the dam’s intersection point on the stream
    current_point = nearest_points(closest_stream.geometry, dam_point)[0]
    # Assume current_link is the stream segment containing current_point (the dam intersection)
    current_link = start_link

    dam_ids = []
    link_nos = []
    points_of_interest = []

    weir_length = dam_gdf['weir_length'].values[0]

    # Loop for each downstream cross-section
    for i in range(1, number_of_cross_sections + 1):
        # We want to move downstream the length of the weir from the current point
        remaining_distance_to_travel = weir_length

        print(f"Calculating point {i*remaining_distance_to_travel} meters downstream of the dam.")
        
        # Traverse the network until we've moved the full tw meters
        while remaining_distance_to_travel > 0:
            # Get downstream edges from the current link
            downstream_edges = list(G.out_edges(current_link, data=True))
            if not downstream_edges:
                raise ValueError(
                    f"Not enough downstream stream length for cross-section {i} (link {current_link})."
                )
            # Select the first valid downstream edge (adjust this logic as needed)
            next_link = None
            for edge in downstream_edges:
                _, candidate_next_link, _ = edge
                if candidate_next_link in G:
                    next_link = candidate_next_link
                    break
            if next_link is None:
                raise ValueError(f"No valid downstream link from link {current_link}.")

            # Retrieve geometry for the current segment (from current_link to next_link)
            seg_geom = shape(G[current_link][next_link]['geometry'])
            # seg_coords = list(seg_geom.coords)

            if seg_geom.geom_type.startswith('Multi'):
                merged_geom = linemerge(seg_geom)
                seg_coords = list(merged_geom.coords)
            else:
                seg_coords = list(seg_geom.coords)


            # Ensure the segment’s coordinates are ordered so that its start is near our current_point.
            # (Using almost_equals avoids floating-point comparison issues.)
            if not Point(seg_coords[0]).equals_exact(current_point, tolerance=0.2):
                seg_coords.reverse()
            seg_line = LineString(seg_coords)

            # Find how far along seg_line our current_point lies
            proj_distance = seg_line.project(current_point)
            distance_remaining_in_seg = seg_line.length - proj_distance

            if distance_remaining_in_seg >= remaining_distance_to_travel:
                # The target point lies within the current segment.
                new_point = seg_line.interpolate(proj_distance + remaining_distance_to_travel)
                # Update our current point and finish the "walk" for this cross-section.
                current_point = new_point
                remaining_distance_to_travel = 0
            else:
                # Use up the remainder of this segment and continue into the next one.
                remaining_distance_to_travel -= distance_remaining_in_seg
                # Set current_point to the end of the segment.
                current_point = Point(seg_coords[-1])
                # Update current_link to the downstream link we just traversed.
                current_link = next_link

        # At this point, current_point has moved exactly tw meters from its previous location.
        # Record the result (convert to EPSG:4326 if needed).
        downstream_point = (
            gpd.GeoSeries([current_point], crs=projected_crs)
            .geometry.iloc[0]
        )
        dam_ids.append(dam_id)
        link_nos.append(current_link)
        points_of_interest.append(downstream_point)

    # downstream_points now contains cross-sections spaced exactly tw meters apart along the stream.
    # Compute lat/lon for curve data
    (minx, miny, maxx, maxy, dx, dy, _, _, _, _) = Get_Raster_Details(DEM_File)
    cellsize_x, cellsize_y = abs(float(dx)), abs(float(dy))
    lat_base, lon_base = float(maxy) - 0.5 * cellsize_y, float(minx) + 0.5 * cellsize_x

    curve_data_df['Lat'] = lat_base - curve_data_df['Row'] * cellsize_y
    curve_data_df['Lon'] = lon_base + curve_data_df['Col'] * cellsize_x
    vdt_df['Lat'] = lat_base - vdt_df['Row'] * cellsize_y
    vdt_df['Lon'] = lon_base + vdt_df['Col'] * cellsize_x

    curve_data_gdf = gpd.GeoDataFrame(curve_data_df, geometry=gpd.points_from_xy(curve_data_df['Lon'], curve_data_df['Lat']), crs="EPSG:4269")
    vdt_gdf = gpd.GeoDataFrame(vdt_df, geometry=gpd.points_from_xy(vdt_df['Lon'], vdt_df['Lat']), crs="EPSG:4269")

    # Convert all to projected CRS
    curve_data_gdf, vdt_gdf = (gdf.to_crs(projected_crs) for gdf in [curve_data_gdf, vdt_gdf])

    # use upstream_elevation_change_threshold and distance_upstream_increment to find upstream cross-section
    # initial intersection point at dam
    current_point_geom = nearest_points(closest_stream.geometry, dam_point)[0]
    current_link = start_link


    # store dam origin so we can always measure from it
    origin_point_geom = current_point_geom
    origin_link       = current_link

    # sample initial BaseElev
    init_idx = curve_data_gdf.geometry.distance(current_point_geom).idxmin()
    last_base_elev = curve_data_gdf.at[init_idx, 'BaseElev']

    increment = 1
    while True:
        # compute cumulative distance from the dam
        distance = distance_upstream_increment * increment  
        # always start from the original dam point/link
        current_point_geom, current_link = move_upstream(
            origin_point_geom,
            origin_link,
            distance,
            G
        )

        # reproject to lat/lon if needed, then record
        upstream_point = (
            gpd.GeoSeries([current_point_geom], crs=projected_crs)
            .geometry.iloc[0]
        )

        # find nearest BaseElev in curve_data_gdf
        idx = curve_data_gdf.geometry.distance(current_point_geom).idxmin()
        current_base_elev = curve_data_gdf.at[idx, 'BaseElev']


        # check elevation change
        elev_diff = abs(current_base_elev - last_base_elev)
        if elev_diff >= upstream_elevation_change_threshold:
            print(f"Reached threshold to find upstream cross-section:"
                  f"distance = {distance:.3f} and Δ elevation = {elev_diff:.3f}")
            break

        # update for next iteration
        last_base_elev = current_base_elev
        increment += 1

    # record the upstream that matches the elevation change threshold we were looking for
    dam_ids.append(dam_id)
    link_nos.append(current_link)
    points_of_interest.append(upstream_point)

    # **Save the Downstream Points as a Shapefile**
    downstream_gdf = gpd.GeoDataFrame(data={'dam_id': dam_ids, rivid_field: link_nos, 'geometry': points_of_interest},
                                      crs=projected_crs)


    # **4. Find Nearest VDT and Curve Data Points to the points of interest**
    vdt_gdfs = []
    curve_data_gdfs = []
    # Extract target point from the GeoDataFrame
    for pt in downstream_gdf.geometry:
        # compute distances
        dists_curve = curve_data_gdf.geometry.distance(pt)
        dists_vdt   = vdt_gdf.geometry.distance(pt)

        # grab the *index* of the nearest cell
        idx_curve = dists_curve.idxmin()
        idx_vdt   = dists_vdt.idxmin()

        # select exactly that one row
        nearest_curves = curve_data_gdf.loc[[idx_curve]]
        nearest_vdt    = vdt_gdf.loc[[idx_vdt]]

        curve_data_gdfs.append(nearest_curves)
        vdt_gdfs.append(nearest_vdt)
    
    # combine the VDT gdfs and curve data gdfs into one a piece
    vdt_gdf = pd.concat(vdt_gdfs)
    curve_data_gdf = pd.concat(curve_data_gdfs)

    # Open the cross-section file and read the end point lat and lons
    XS_Out_File_df = pd.read_csv(XS_Out_File, sep='\t')

    (minx, miny, maxx, maxy, dx, dy, _, _, _, _) = Get_Raster_Details(DEM_File)
    cellsize_x, cellsize_y = abs(float(dx)), abs(float(dy))
    lat_base, lon_base = float(maxy) - 0.5 * cellsize_y, float(minx) + 0.5 * cellsize_x

    XS_Out_File_df['Lat1'] = lat_base - XS_Out_File_df['r1'] * cellsize_y
    XS_Out_File_df['Lon1'] = lon_base + XS_Out_File_df['c1'] * cellsize_x
    XS_Out_File_df['Lat2'] = lat_base - XS_Out_File_df['r2'] * cellsize_y
    XS_Out_File_df['Lon2'] = lon_base + XS_Out_File_df['c2'] * cellsize_x

    curve_data_df['Lat'] = lat_base - curve_data_df['Row'] * cellsize_y
    curve_data_df['Lon'] = lon_base + curve_data_df['Col'] * cellsize_x
    vdt_df['Lat'] = lat_base - vdt_df['Row'] * cellsize_y
    vdt_df['Lon'] = lon_base + vdt_df['Col'] * cellsize_x

    # filter the XS_Out_File_df to only include cross-sections with 'row' and 'col' combinations that match the curve_data_gdf
    # keep only rows whose (Row,Col) is in curve_data_gdf
    XS_Out_File_df = (XS_Out_File_df.merge(curve_data_gdf[['COMID','Row','Col']].drop_duplicates(),
                                           on=['COMID','Row','Col'], how='inner'))

    # use lat lon pairs to create a line vector shapefile
    xs_lines = []
    for index, row in XS_Out_File_df.iterrows():
        start_point = Point(row['Lon1'], row['Lat1'])
        end_point = Point(row['Lon2'], row['Lat2'])
        xs_lines.append(LineString([start_point, end_point]))

    # Create a GeoDataFrame from the lines, set the CRS to the same as the dem
    XS_Out_File_df['geometry'] = xs_lines
    # Build the CRS object as you already do
    crs = CRS.from_wkt(Rast_Projection)
    # Create the GeoDataFrame using all columns from XS_Out_File_df, including geometry
    xs_gdf = gpd.GeoDataFrame(XS_Out_File_df, geometry='geometry', crs=crs)
    # # convert the crs of the xs_gdf to match the vdt_gdf and curve_data_gdf
    xs_gdf = xs_gdf.to_crs(projected_crs)

    return downstream_gdf, vdt_gdf, curve_data_gdf, xs_gdf


def Dam_Assessment(DEM_Folder: str, DEM: str, ESA_LC_Folder: str, STRM_Folder: str, LAND_Folder: str,
                   FLOW_Folder: str, VDT_Folder: str, ARC_Folder: str, BathyFileFolder: str, XS_Folder: str,
                   ManningN: str, bathy_use_banks: bool, find_banks_based_on_landcover: bool,
                   create_reach_average_curve_file: bool, dam_csv: str, dam_id_field: str, dam_id: int|str,
                   known_baseflow: float, known_channel_forming_discharge: float,
                   upstream_elevation_change_threshold=1.0, StrmShp_gdf: GeoDataFrame=None):

    if DEM.endswith((".tif", ".img")):
        DEM_Name = DEM
        # FileName = os.path.splitext(DEM_Name)[0] # grabs just name without extension
        DEM_File = os.path.join(DEM_Folder, DEM_Name)
        
        #Input Dataset
        ARC_FileName_Bathy = os.path.join(ARC_Folder, f'ARC_Input_{dam_id}.txt')
        
        #Datasets to be Created
        Dam_StrmShp = os.path.join(STRM_Folder, f"{dam_id}_StrmShp.shp")
        Dam_Reanalsyis_FlowFile = os.path.join(FLOW_Folder,f"{dam_id}_Reanalysis.csv")

        STRM_File = os.path.join(STRM_Folder, f'{dam_id}_STRM_Raster.tif')
        STRM_File_Clean = STRM_File.replace('.tif','_Clean.tif')
        LAND_File = os.path.join(LAND_Folder, f'{dam_id}_LAND_Raster.tif')
        
        VDT_File = os.path.join(VDT_Folder, f'{dam_id}_VDT_Database.txt')
        Curve_File = os.path.join(VDT_Folder, f'{dam_id}_CurveFile.csv')

        # these are the files that will be created by the code
        Local_VDT_File = os.path.join(VDT_Folder, f'{dam_id}_Local_VDT_Database.shp')
        Local_Curve_File = os.path.join(VDT_Folder, f'{dam_id}_Local_CurveFile.shp')
        
        ARC_BathyFile = os.path.join(BathyFileFolder, f'{dam_id}_ARC_Bathy.tif')

        XS_Out_File = os.path.join(XS_Folder, f'{dam_id}_XS_Out.txt')
        Local_XS_File = os.path.join(XS_Folder, f'{dam_id}_Local_XS_Lines.shp')

        # Get the details of the DEM file
        lon_1, lat_1, lon_2, lat_2, dx, dy, ncols, nrows, geoTransform, Rast_Projection = Get_Raster_Details(DEM_File)

        # set up a WGS84 ellipsoid
        geod = Geod(ellps='WGS84')

        # horizontal distance: from (lon1,lat1) to (lon1+dx, lat1)
        _, _, res_x_m = geod.inv(lon_1, lat_1, lon_1 + dx, lat_1)

        # vertical distance: from (lon1,lat1) to (lon1, lat1+dy)
        _, _, res_y_m = geod.inv(lon_1, lat_1, lon_1, lat_1 + dy)

        print(f"Dam Assessment: Exact pixel size of DEM: {res_x_m:.3f} m × {res_y_m:.3f} m")

        # Calculate the average distance increment for upstream processing
        distance_upstream_increment = (res_x_m + res_y_m)/2 
        
        #Download and Process Land Cover Data
        LandCoverFile = ''
        if not os.path.exists(LAND_File):

            # Get geometry in original projection
            geom = esa.Get_Polygon_Geometry(lon_1, lat_1, lon_2, lat_2)

            # Check if raster projection is WGS 84
            raster_crs = CRS.from_wkt(Rast_Projection)
            wgs84_crs = CRS.from_epsg(4326)

            if raster_crs != wgs84_crs:
                # Convert geometry to WGS 84
                transformer = Transformer.from_crs(raster_crs, wgs84_crs, always_xy=True)
                geom = transform(transformer.transform, geom)

            LandCoverFile = esa.Download_ESA_WorldLandCover(ESA_LC_Folder, geom, 2021)

        # This function sets-up the Input files for ARC and FloodSpreader
        # It also does some of the geospatial processing
        ARC_FileName_Bathy, DEM_Reanalsyis_FlowFile, Dam_StrmShp \
            = Process_Geospatial_Data(ARC_FileName_Bathy, DEM_File, LandCoverFile, STRM_File, STRM_File_Clean,
                                      LAND_File, ManningN, VDT_File, Curve_File, ARC_BathyFile, Dam_StrmShp,
                                      Dam_Reanalsyis_FlowFile, XS_Out_File, bathy_use_banks,
                                      find_banks_based_on_landcover, create_reach_average_curve_file, dam_csv,
                                      dam_id_field, dam_id, known_baseflow, known_channel_forming_discharge, StrmShp_gdf)

        # read in the reanalysis streamflow and break the code if the dataframe is empty or if the streamflow is all 0
        DEM_Reanalsyis_FlowFile_df = pd.read_csv(DEM_Reanalsyis_FlowFile)
        if (DEM_Reanalsyis_FlowFile_df.empty
                or DEM_Reanalsyis_FlowFile_df['qout_max'].mean() <= 0
                or len(DEM_Reanalsyis_FlowFile_df.index) == 0):

            print(f"Dam_Assessment: Results for {DEM} are not possible because we don't have streamflow estimates...")
            return

        
        # Create our Curve and VDT Database data
        if not os.path.exists(ARC_BathyFile) or not os.path.exists(VDT_File) or not os.path.exists(Curve_File):
            print(f'Dam_Assessment: Cannot find bathy file, so creating {ARC_BathyFile}')
            arc = Arc(ARC_FileName_Bathy)
            arc.run() # Runs ARC
        
        # Now we need to use the Dam_StrmShp and VDT data to find the stream cells at distance increments below the dam
        downstream_gdf, vdt_gdf, curve_data_gdf, xs_gdf = find_stream_cells_at_increments_above_and_below_dam(Curve_File, VDT_File, XS_Out_File, dam_csv, dam_id_field, dam_id,
                                                                                                      Dam_StrmShp, Dam_Reanalsyis_FlowFile,
                                                                                                      DEM_File, upstream_elevation_change_threshold, 
                                                                                                      distance_upstream_increment, known_baseflow, Rast_Projection)
        
        # output the results to shapefiles
        vdt_gdf.to_file(Local_VDT_File)
        curve_data_gdf.to_file(Local_Curve_File)
        xs_gdf.to_file(Local_XS_File)


    return


def process_dam(dam_dict):

    # set if the system will use banks of water surface elevation to estimate bathymetry
    bathy_use_banks = dam_dict['bathy_use_banks']

    # if you don't have the stream network preprocessed, do so
    process_stream_network = dam_dict['process_stream_network']

    # set if you want to use the landcover data to find the banks of the river, instead of the flat water surface elevation in the DEM
    find_banks_based_on_landcover = dam_dict['find_banks_based_on_landcover']

    # let's tell ARC whether we want the curvefile parameters to be the same for each reach or vary by stream cell
    create_reach_average_curve_file = dam_dict['create_reach_average_curve_file']

    # if we have a known baseflow, we will add it to the streamflow file and use it to estimate bathymetry
    known_baseflow = dam_dict['known_baseflow']

    # if we have a known channel forming discharge, we will add it to the streamflow file and use it to estimate bathymetry
    known_channel_forming_discharge = dam_dict['known_channel_forming_discharge']

    # the elevation change threshold to use when identifying upstream cross-section
    upstream_elevation_change_threshold = dam_dict['upstream_elevation_change_threshold']


    #Folder Management
    Output_Dir = dam_dict['output_dir']
    dam = dam_dict['name']
    dam_csv = dam_dict['dam_csv']
    dam_id_field = dam_dict['dam_id_field']
    dam_id = dam_dict['dam_id']
    DEM_Folder = dam_dict['dem_dir']
    ARC_Folder = os.path.join(Output_Dir, f'{dam}', 'ARC_InputFiles')
    BathyFileFolder = os.path.join(Output_Dir, f'{dam}', 'Bathymetry')
    STRM_Folder = os.path.join(Output_Dir, f'{dam}', 'STRM')
    LAND_Folder = os.path.join(Output_Dir, f'{dam}', 'LAND')
    FLOW_Folder = os.path.join(Output_Dir, f'{dam}', 'FLOW')
    VDT_Folder = os.path.join(Output_Dir, f'{dam}', 'VDT')
    ESA_LC_Folder = os.path.join(Output_Dir, f'{dam}', 'ESA_LC')
    XS_Folder = os.path.join(Output_Dir, f'{dam}', 'XS')

    # Create Folders
    Create_Folder(ESA_LC_Folder)
    Create_Folder(STRM_Folder)
    Create_Folder(LAND_Folder)
    Create_Folder(FLOW_Folder)
    Create_Folder(VDT_Folder)
    Create_Folder(ARC_Folder)
    Create_Folder(BathyFileFolder)
    Create_Folder(XS_Folder)

    # Datasets that can be good for a large domain
    StrmSHP = dam_dict['flowline']
    ManningN = os.path.join(LAND_Folder, 'AR_Manning_n_MED.txt')

    # Create a Baseline Manning N File
    print(f'Creating Manning n file: {ManningN}')
    Create_BaseLine_Manning_n_File_ESA(ManningN)

    # This is the list of all the DEM files we will go through
    DEM_List = os.listdir(DEM_Folder)

    # Before we get too far ahead, let's make sure that our DEMs and Flowlines have the same coordinate system
    # we will assume that all DEMs in the DEM list have the same coordinate system
    StrmShp_gdf = gpd.GeoDataFrame()

    if process_stream_network:
        print(f'Reading in stream file: {StrmSHP}')

        # Get the first DEM file to determine extent and CRS
        test_dem = next((file for file in DEM_List if file.endswith('.tif')), None)
        if test_dem is None:
            raise ValueError("No .tif files found in DEM directory")

        test_dem_path = os.path.join(DEM_Folder, test_dem)

        # Get DEM bounds and CRS using GDAL (more efficient than opening full raster)
        raster_dataset = gdal.Open(test_dem_path)
        gt = raster_dataset.GetGeoTransform()

        # Get the bounds of the raster (xmin, ymin, xmax, ymax)
        xmin = gt[0]
        xmax = xmin + gt[1] * raster_dataset.RasterXSize
        ymin = gt[3] + gt[5] * raster_dataset.RasterYSize
        ymax = gt[3]

        # Get DEM CRS
        dem_proj = raster_dataset.GetProjection()
        dem_spatial_ref = osr.SpatialReference()
        dem_spatial_ref.ImportFromWkt(dem_proj)
        dem_spatial_ref.AutoIdentifyEPSG()
        dem_epsg_code = dem_spatial_ref.GetAuthorityCode(None)

        # Close the raster dataset
        del raster_dataset

        # Create bounding box in DEM CRS
        dem_bbox = (xmin, ymin, xmax, ymax)

        dem_bounds_geom = box(*dem_bbox)
        dem_bbox_gdf = gpd.GeoDataFrame(geometry=[dem_bounds_geom], crs=f"EPSG:{dem_epsg_code}")


        print(f'DEM bounds: {dem_bbox}')
        print(f'DEM CRS: {dem_epsg_code}')

        # Read stream file with bbox filter (much more efficient!)
        if StrmSHP.endswith(".gdb"):
            # Read the layer from the geodatabase with bbox
            with fiona.Env():
                with fiona.open(StrmSHP, layer="geoglowsv2") as src:
                    stream_crs = src.crs
        elif StrmSHP.endswith(".shp") or StrmSHP.endswith(".gpkg"):
            stream_meta = gpd.read_file(StrmSHP, rows=1)
            stream_crs = stream_meta.crs
        else:
            stream_crs = f"EPSG:{dem_epsg_code}"

        dem_bbox_gdf = dem_bbox_gdf.to_crs(stream_crs)
        bbox_geom = dem_bbox_gdf.geometry.iloc[0].bounds

        if StrmSHP.endswith(".gdb"):
            layer_name = "geoglowsv2"
        elif "NHD" in StrmSHP:
            layer_name = "NHDFlowline"
        else:
            layer_name = None

        StrmShp_gdf = gpd.read_file(StrmSHP, layer=layer_name, engine='pyogrio', bbox=bbox_geom)

        print(f'Filtered stream features: {len(StrmShp_gdf)} (vs reading entire file)')

        # removing any lingering NoneType geometries
        StrmShp_gdf = StrmShp_gdf[~StrmShp_gdf.geometry.isna()]

        print('Converting the coordinate system of the stream file to match the DEM files, if necessary')

        # Check if stream CRS matches DEM CRS
        if StrmShp_gdf.crs != f"EPSG:{dem_epsg_code}":
            print("DEM and Stream Network have different coordinate systems...")
            print(f"Stream CRS: {StrmShp_gdf.crs}")
            print(f"DEM CRS: EPSG:{dem_epsg_code}")
            # Reproject the stream data to match DEM CRS
            StrmShp_gdf = StrmShp_gdf.to_crs(f"EPSG:{dem_epsg_code}")

    # Now go through each DEM dataset
    for DEM in DEM_List:
        Dam_Assessment(DEM_Folder, DEM, ESA_LC_Folder, STRM_Folder, LAND_Folder, FLOW_Folder,
                       VDT_Folder, ARC_Folder, BathyFileFolder, XS_Folder, ManningN, bathy_use_banks,
                       find_banks_based_on_landcover, create_reach_average_curve_file,
                       dam_csv, dam_id_field, dam_id, known_baseflow, known_channel_forming_discharge,
                       upstream_elevation_change_threshold, StrmShp_gdf)

    # delete the ESA_LC_Folder and the data in it
    # Loop through all files in the directory and remove them
    for file in Path(ESA_LC_Folder).glob("*"):
        try:
            if file.is_file():
                # Adjust file permissions before deletion
                if platform.system() == "Windows":
                    os.chmod(file, stat.S_IWRITE)  # Remove read-only attribute on Windows
                else:
                    os.chmod(file, stat.S_IWUSR)  # Give user write permission on Unix systems
                os.remove(file)
                print(f"process_dam: Deleted file: {file}")
        except Exception as e:
            print(f"Error deleting file {file}: {e}")
    if os.path.exists(ESA_LC_Folder):
        # Adjust file permissions before deletion
        if platform.system() == "Windows":
            os.chmod(ESA_LC_Folder, stat.S_IWRITE)  # Remove read-only attribute on Windows
        else:
            os.chmod(ESA_LC_Folder, stat.S_IWUSR)  # Give user write permission on Unix systems
        os.rmdir(ESA_LC_Folder)
        print(f"process_dam: Deleted empty folder: {ESA_LC_Folder}")
    else:
        print(f"process_dam: Folder {ESA_LC_Folder} does not exist.")

    
def process_json_input(json_file):
    """
        Process input from a JSON file.
    """
    with open(json_file, 'r') as file:
        print(f'Opening {file}')
        data = json.load(file)
        print(data)
    
    dams = data.get("dams", [])
    for dam in dams:
        dam_name = dam.get("name")
        dam_csv = os.path.normpath(dam.get("dam_csv"))
        flowline = os.path.normpath(dam.get("flowline"))
        dem_dir = os.path.normpath(dam.get("dem_dir"))
        output_dir = os.path.normpath(dam.get("output_dir"))

        dam_dict = {
            "name": dam_name,
            "dam_csv": dam_csv,
            "dam_id_field": dam.get("dam_id_field"),
            "dam_id": int(dam.get("dam_id")), 
            "flowline": flowline,
            "dem_dir": dem_dir,
            "output_dir": output_dir,
            "bathy_use_banks": dam.get("bathy_use_banks", False),
            "flood_waterlc_and_strm_cells":dam.get("flood_waterlc_and_strm_cells", False),
            "process_stream_network": dam.get("process_stream_network", False),
            "find_banks_based_on_landcover": dam.get("find_banks_based_on_landcover", True),
            "create_reach_average_curve_file": dam.get("create_reach_average_curve_file", False),
            "known_baseflow": dam.get("known_baseflow", None),
            "known_channel_forming_discharge": dam.get("known_channel_forming_discharge", None),
            "upstream_elevation_change_threshold": dam.get("upstream_elevation_change_threshold", 1.0)
        }

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        print(f"Processing dam: {dam_name} with parameters: {dam_dict}")

        # Call your existing processing logic here
        process_dam(dam_dict)


def normalize_path(path):
    return os.path.normpath(path)


def process_cli_arguments(args):
    """Process input from CLI arguments."""
    output_dir = args.output_dir
    dam_name = args.watershed
    dam_dict = {
        "name": dam_name,
        "dam_csv": normalize_path(args.dam_csv),
        "dam_id_field": args.dam_id_field,
        "dam_id": args.dam_id,
        "flowline": normalize_path(args.flowline),
        "dem_dir": normalize_path(args.dem_dir),
        "bathy_use_banks": args.bathy_use_banks,
        "output_dir": normalize_path(output_dir),
        "process_stream_network": args.process_stream_network,
        "find_banks_based_on_landcover": args.find_banks_based_on_landcover,
        "create_reach_average_curve_file": args.create_reach_average_curve_file,
        "known_baseflow": args.known_baseflow,
        "known_channel_forming_discharge": args.known_channel_forming_discharge,
    }

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing dam: {args.watershed} with parameters: {dam_dict}")
    print(f"Results will be saved in: {output_dir}")

    # Call the existing processing logic here
    process_dam(dam_dict)


def main():
    parser = argparse.ArgumentParser(description="Process rating curves on streams below a dam.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for JSON input
    json_parser = subparsers.add_parser("json", help="Process watersheds from a JSON file")
    json_parser.add_argument("json_file", type=str, help="Path to the JSON file")

    # Subcommand for CLI input
    cli_parser = subparsers.add_parser("cli", help="Process watershed parameters via CLI")
    cli_parser.add_argument("dam", type=str, help="Dam name")
    cli_parser.add_argument("dam_csv", type=str, help="Path to the dam csv file")
    cli_parser.add_argument("dam_id_field", type=str, help="Name of the csv field with the dam id")   
    cli_parser.add_argument("dam_id", type=int, help="ID of the dam in the damn_id_field")  
    cli_parser.add_argument("flowline", type=str, help="Path to the flowline shapefile")
    cli_parser.add_argument("dem_dir", type=str, help="Directory containing DEM files")
    cli_parser.add_argument("output_dir", type=str, help="Directory where results will be saved")
    cli_parser.add_argument("--bathy_use_banks", action="store_true", help="Use bathy banks for processing")
    cli_parser.add_argument("--process_stream_network", action="store_true", help="Clean DEM data before processing")
    cli_parser.add_argument("--find_banks_based_on_landcover", action="store_true", help="Use landcover data for finding banks when estimating bathymetry")
    cli_parser.add_argument("--create_reach_average_curve_file", action="store_true", help="Create a reach average curve file instead of one that varies for each stream cell")
    cli_parser.add_argument("--known_baseflow", type=float, default=None, help="Known baseflow value")
    cli_parser.add_argument("--known_channel_forming_discharge", type=float, default=None, help="Known channel forming discharge value")
    cli_parser.add_argument("--upstream_elevation_change_threshold", type=float, help="The upstream elevation change used to identify the appropriate upstream cross-section, default is 0.5 meters", default=1.0)

    args = parser.parse_args()

    if args.command == "json":
        print(f'Processing {args.json_file}')
        process_json_input(args.json_file)
    elif args.command == "cli":
        process_cli_arguments(args)

if __name__ == "__main__":
    main()








