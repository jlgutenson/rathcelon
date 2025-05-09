'''
Here is a list of where Low Head Dams are located.  What we're wanting to do:

1.) Identify the stream reach(es) that are associated with the Low-Head Dam.
2.) Determine the approximate Top-Width (TW) of the Dam (or the stream itself)
3.) Determine the average (base) flow for the stream as well as the seasonal high flow (not a flood, just a normal high-flow).
4.) Go downstream approximately TW distance and pull a perpendicular cross-section.
5.) Go downstream another TW distance (2*TW from the dam) and pull another perpendicular cross-section.
6.) Go downstream another TW distance (3*TW from the dam) and pull another perpendicular cross-section.
7.) For each cross-section estimate the bathymetry.
8.) For each cross-section calculate the rating curve of the cross-section.  Slope can likely be calculated from steps 4-6.

'''

# build-in imports
import argparse
import multiprocessing as mp
import sys
import os
import re
import subprocess
from datetime import datetime, timedelta

# third-party imports
from arc import Arc
try:
    import gdal 
except: 
    from osgeo import gdal, ogr, osr
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import rasterize
import geopandas as gpd
import fiona
from shapely.geometry import shape
import networkx as nx
from pathlib import Path
import platform
from pyproj import CRS, Transformer
from shapely.geometry import Point, LineString, MultiLineString, shape, mapping
from shapely.ops import nearest_points, linemerge, split, transform
import shutil
import stat
import numpy as np

# local imports
from . import streamflow_processing as HistFlows
from . import esa_download_processing as ESA


def Process_Geospatial_Data(ARC_Folder,
                            ARC_FileName_Bathy, 
                            DEM_File, LandCoverFile, STRM_File, 
                            STRM_File_Clean, LAND_File, 
                            BathyFileFolder, FLOW_Folder, XS_Folder, ManningN,
                            VDT_File, Curve_File, ARC_BathyFile, Dam_StrmShp, 
                            Dam_Reanalsyis_FlowFile, XS_File_Out, bathy_use_banks,
                            find_banks_based_on_landcover, create_reach_average_curve_file,
                            dam_csv, dam_id_field, dam_id, known_baseflow, known_channel_forming_discharge,
                            StrmShp_gdf=None):

    
  
    
    #Get the Spatial Information from the DEM Raster
    (minx, miny, maxx, maxy, dx, dy, ncols, nrows, dem_geoTransform, dem_projection) = Get_Raster_Details(DEM_File)
    projWin_extents = [minx, maxy, maxx, miny]
    outputBounds = [minx, miny, maxx, maxy]  #https://gdal.org/api/python/osgeo.gdal.html
   
    #Create Land Dataset
    if os.path.isfile(LAND_File):
        print(LAND_File + ' Already Exists')
    else: 
        print('Creating ' + LAND_File) 
        # Let's make sure all the GIS data is using the same coordinate system as the DEM
        LandCoverFile = Check_and_Change_Coordinate_Systems(DEM_File, LandCoverFile)
        Create_AR_LandRaster(LandCoverFile, LAND_File, projWin_extents, dem_projection, ncols, nrows)

    # now we need to figure out if our Dam_StrmShp and DEM_Reanalysis_Flowfile exists and if not, create it
    if os.path.isfile(Dam_StrmShp) and os.path.isfile(Dam_Reanalsyis_FlowFile):
        print(Dam_StrmShp + ' Already Exists')
        print(Dam_Reanalsyis_FlowFile + ' Already Exists')
        Dam_StrmShp_gdf = gpd.read_file(Dam_StrmShp)
        rivids = Dam_StrmShp_gdf['LINKNO'].values
    elif StrmShp_gdf is not None and os.path.isfile(Dam_StrmShp) is False and os.path.isfile(Dam_Reanalsyis_FlowFile) is False:
        print('Running Function: Process_and_Write_Retrospective_Data_for_Dam')
        (DEM_Reanalsyis_FlowFile, Dam_StrmShp, rivids, Dam_StrmShp_gdf) = HistFlows.Process_and_Write_Retrospective_Data_for_Dam(StrmShp_gdf, 'LINKNO', dam_csv, 
                                                                                                                                 dam_id_field, dam_id, known_baseflow,
                                                                                                                                 known_channel_forming_discharge,
                                                                                                                                 Dam_Reanalsyis_FlowFile,
                                                                                                                                 Dam_StrmShp)


    #Create Stream Raster
    if os.path.isfile(STRM_File):
        print(STRM_File + ' Already Exists')
    else:
        print('Creating ' + STRM_File)
        print('   by rasterizing ' + Dam_StrmShp)
        Create_AR_StrmRaster(Dam_StrmShp, STRM_File, outputBounds, minx, miny, maxx, maxy, dx, dy, ncols, nrows, 'LINKNO')
    
    #Clean Stream Raster
    if os.path.isfile(STRM_File_Clean):
        print(STRM_File_Clean + ' Already Exists')
    else:
        print('Creating ' + STRM_File_Clean)
        Clean_STRM_Raster(STRM_File, STRM_File_Clean)
    
    
    #Get the unique values for all the stream ids
    (S, ncols, nrows, cellsize, yll, yur, xll, xur, lat, dem_geotransform, dem_projection) = Read_Raster_GDAL(STRM_File_Clean)
    (RR,CC) = S.nonzero()
    num_strm_cells = len(RR)
    COMID_Unique = np.unique(S)
    # COMID_Unique = np.delete(COMID_Unique, 0)  #We don't need the first entry of zero
    COMID_Unique = COMID_Unique[np.where(COMID_Unique > 0)]
    COMID_Unique = np.sort(COMID_Unique).astype(int)
    num_comids = len(COMID_Unique)


    #Create the Bathy Input File

    # Let's extract the weir length real quick
    dam_gdf = pd.read_csv(dam_csv)
    # Filter to the specific dam
    dam_gdf = dam_gdf[dam_gdf[dam_id_field] == dam_id]
    weir_length = dam_gdf['weir_length'].values[0]

    print('Creating ARC Input File: ' + ARC_FileName_Bathy)

    if bathy_use_banks is False and known_baseflow is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'qout_median'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True and known_baseflow is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'rp2'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is False and known_baseflow is not None and known_channel_forming_discharge is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_baseflow'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True and known_channel_forming_discharge is not None and known_baseflow is None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_channel_forming_discharge'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is False and known_baseflow is not None and known_channel_forming_discharge is not None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_baseflow'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, 
                                          COMID_Param, Q_BF_Param, Q_Param, 
                                          STRM_File_Clean, LAND_File, Dam_Reanalsyis_FlowFile, 
                                          VDT_File, Curve_File, XS_File_Out, ManningN, ARC_BathyFile,
                                          Dam_StrmShp, bathy_use_banks, weir_length,
                                          find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True and known_channel_forming_discharge is not None and known_baseflow is not None:
        COMID_Param = 'COMID'
        Q_BF_Param = 'known_channel_forming_discharge'
        Q_Param = 'rp100_premium'
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

def Create_FlowFile(MainFlowFile, FlowFileName, OutputID, Qparam):
    infile = open(MainFlowFile,'r')
    lines = infile.readlines()
    ls = lines[0].strip().split(',')
    q_val = 0
    c_val = 0
    for i in range(len(ls)):
        if ls[i]==Qparam:
            q_val=i
        if ls[i]==OutputID:
            c_val=i
    
    outfile = open(FlowFileName, 'w')
    outfile.write(OutputID + ',' + Qparam)
    
    for r in range(1,len(lines)):
        ls = lines[r].strip().split(',')
        out_str = '\n' + ls[c_val] + ',' + ls[q_val]
        outfile.write(out_str)
    outfile.close()
    return

def Create_Folder(F):
    if not os.path.exists(F): 
        os.makedirs(F)
    return

def Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, COMID_Param, Q_BF_Param, Q_Param, 
                                      STRM_File_Clean, LAND_File, DEM_Reanalsyis_FlowFile, VDT_File, 
                                      Curve_File, XS_Out_File, ManningN, ARC_BathyFile,
                                      DEM_StrmShp, bathy_use_banks, weir_length, find_banks_based_on_landcover,
                                      create_reach_average_curve_file):
    
    out_file = open(ARC_FileName_Bathy,'w')
    out_file.write('#ARC_Inputs')
    out_file.write('\n' + 'DEM_File	' + DEM_File)
    out_file.write('\n' + 'Stream_File	' + STRM_File_Clean)
    out_file.write('\n' + 'LU_Raster_SameRes	' + LAND_File)
    out_file.write('\n' + 'LU_Manning_n	' + ManningN)
    out_file.write('\n' + 'Flow_File	' + DEM_Reanalsyis_FlowFile)
    out_file.write('\n' + 'Flow_File_ID	' + COMID_Param)
    out_file.write('\n' + 'Flow_File_BF	' + Q_BF_Param)
    out_file.write('\n' + 'Flow_File_QMax	' + Q_Param)
    out_file.write('\n' + 'Spatial_Units	deg')
    x_section_dist = round(2 * weir_length,0) # set the cross-section distance to two x weir length
    out_file.write('\n' + 'X_Section_Dist	' + str(x_section_dist))
    out_file.write('\n' + 'Degree_Manip	6.1')
    out_file.write('\n' + 'Degree_Interval	1.5')
    out_file.write('\n' + 'Low_Spot_Range	2')
    out_file.write('\n' + 'Str_Limit_Val	1')
    out_file.write('\n' + 'Gen_Dir_Dist	10')
    out_file.write('\n' + 'Gen_Slope_Dist	10')
    
    out_file.write('\n\n#VDT_Output_File_and_CurveFile')
    out_file.write('\n' + 'VDT_Database_NumIterations	30')
    out_file.write('\n' + 'VDT_Database_File	' + VDT_File)
    out_file.write('\n' + 'Print_VDT_Database' + '\t' + VDT_File)
    out_file.write('\n' + 'Print_Curve_File' + '\t' + Curve_File)
    out_file.write('\n' + 'Reach_Average_Curve_File' + '\t' + f'{create_reach_average_curve_file}')

    out_file.write('\n\n#Bathymetry_Information')
    out_file.write('\n' + 'Bathy_Trap_H	0.20')
    out_file.write('\n' + 'Bathy_Use_Banks' + '\t' + str(bathy_use_banks))
    if find_banks_based_on_landcover is True:
        out_file.write('\n' + 'FindBanksBasedOnLandCover' + '\t' + str(find_banks_based_on_landcover))
    out_file.write('\n' + 'AROutBATHY	' + ARC_BathyFile)
    out_file.write('\n' + 'BATHY_Out_File	' + ARC_BathyFile)
    out_file.write('\n' + 'XS_Out_File' + '\t' + XS_Out_File)
    out_file.close()
    

def Create_BaseLine_Manning_n_File(ManningN):
    out_file = open(ManningN,'w')
    out_file.write('LC_ID	Description	Manning_n')
    out_file.write('\n' + '11	Water	0.030')
    out_file.write('\n' + '21	Dev_Open_Space	0.013')
    out_file.write('\n' + '22	Dev_Low_Intesity	0.050')
    out_file.write('\n' + '23	Dev_Med_Intensity	0.075')
    out_file.write('\n' + '24	Dev_High_Intensity	0.100')
    out_file.write('\n' + '31	Barren_Land	0.030')
    out_file.write('\n' + '41	Decid_Forest	0.120')
    out_file.write('\n' + '42	Evergreen_Forest	0.120')
    out_file.write('\n' + '43	Mixed_Forest	0.120')
    out_file.write('\n' + '52	Shrub	0.050')
    out_file.write('\n' + '71	Grass_Herb	0.030')
    out_file.write('\n' + '81	Pasture_Hay	0.040')
    out_file.write('\n' + '82	Cultivated_Crops	0.035')
    out_file.write('\n' + '90	Woody_Wetlands	0.100')
    out_file.write('\n' + '95	Emergent_Herb_Wet	0.100')
    out_file.close()

def Create_BaseLine_Manning_n_File_ESA(ManningN):
    out_file = open(ManningN,'w')
    out_file.write('LC_ID	Description	Manning_n')
    out_file.write('\n' + '10	Tree Cover	0.120')
    out_file.write('\n' + '20	Shrubland	0.050')
    out_file.write('\n' + '30	Grassland	0.030')
    out_file.write('\n' + '40	Cropland	0.035')
    out_file.write('\n' + '50	Builtup	0.075')
    out_file.write('\n' + '60	Bare	0.030')
    out_file.write('\n' + '70	SnowIce	0.030')
    out_file.write('\n' + '80	Water	0.030')
    out_file.write('\n' + '90	Emergent_Herb_Wet	0.100')
    out_file.write('\n' + '95	Mangroves	0.100')
    out_file.write('\n' + '100	MossLichen	0.100')
    out_file.close()


def Create_AR_LandRaster(LandCoverFile, LAND_File, projWin_extents, out_projection, ncols, nrows):
    ds = gdal.Open(LandCoverFile)
    ds = gdal.Translate(LAND_File, ds, projWin = projWin_extents, width=ncols, height = nrows)
    ds = None
    return

def Create_AR_StrmRaster(StrmSHP, STRM_File, outputBounds, minx, miny, maxx, maxy, dx, dy, ncols, nrows, Param):
    print(StrmSHP)
    source_ds = gdal.OpenEx(StrmSHP)

    gdal.Rasterize(STRM_File, source_ds, format='GTiff', outputType=gdal.GDT_Int32, outputBounds = outputBounds, width = ncols, height = nrows, noData = -9999, attribute = Param)
    source_ds = None

    return

def Write_Output_Raster(s_output_filename, raster_data, ncols, nrows, dem_geotransform, dem_projection, s_file_format, s_output_type):   
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
    o_output_file = None


def Get_Raster_Details(DEM_File):
    print(DEM_File)
    gdal.Open(DEM_File, gdal.GA_ReadOnly)
    data = gdal.Open(DEM_File)
    geoTransform = data.GetGeoTransform()
    ncols = int(data.RasterXSize)
    nrows = int(data.RasterYSize)
    minx = geoTransform[0]
    dx = geoTransform[1]
    maxy = geoTransform[3]
    dy = geoTransform[5]
    maxx = minx + dx * ncols
    miny = maxy + dy * nrows
    Rast_Projection = data.GetProjectionRef()
    data = None
    return minx, miny, maxx, maxy, dx, dy, ncols, nrows, geoTransform, Rast_Projection


def Read_Raster_GDAL(InRAST_Name):
    try:
        dataset = gdal.Open(InRAST_Name, gdal.GA_ReadOnly)     
    except RuntimeError:
        sys.exit(" ERROR: Field Raster File cannot be read!")
    # Retrieve dimensions of cell size and cell count then close DEM dataset
    geotransform = dataset.GetGeoTransform()
    # Continue grabbing geospatial information for this use...
    band = dataset.GetRasterBand(1)
    RastArray = band.ReadAsArray()
    #global ncols, nrows, cellsize, yll, yur, xll, xur
    ncols=band.XSize
    nrows=band.YSize
    band = None
    cellsize = geotransform[1]
    yll = geotransform[3] - nrows * np.fabs(geotransform[5])
    yur = geotransform[3]
    xll = geotransform[0];
    xur = xll + (ncols)*geotransform[1]
    lat = np.fabs((yll+yur)/2.0)
    Rast_Projection = dataset.GetProjectionRef()
    dataset = None
    print('Spatial Data for Raster File:')
    print('   ncols = ' + str(ncols))
    print('   nrows = ' + str(nrows))
    print('   cellsize = ' + str(cellsize))
    print('   yll = ' + str(yll))
    print('   yur = ' + str(yur))
    print('   xll = ' + str(xll))
    print('   xur = ' + str(xur))
    return RastArray, ncols, nrows, cellsize, yll, yur, xll, xur, lat, geotransform, Rast_Projection

def Clean_STRM_Raster(STRM_File, STRM_File_Clean):
    print('\nCleaning up the Stream File.')
    (SN, ncols, nrows, cellsize, yll, yur, xll, xur, lat, dem_geotransform, dem_projection) = Read_Raster_GDAL(STRM_File)
    
    #Create an array that is slightly larger than the STRM Raster Array
    B = np.zeros((nrows+2,ncols+2))
    
    #Imbed the STRM Raster within the Larger Zero Array
    B[1:(nrows+1), 1:(ncols+1)] = SN
    
    #Added this because sometimes the non-stream values end up as -9999
    B = np.where(B>0,B,0)
    #(RR,CC) = B.nonzero()
    (RR,CC) = np.where(B>0)
    num_nonzero = len(RR)
    
    for filterpass in range(2):
        #First pass is just to get rid of single cells hanging out not doing anything
        p_count = 0
        p_percent = (num_nonzero+1)/100.0
        n=0
        for x in range(num_nonzero):
            if x>=p_count*p_percent:
                p_count = p_count + 1
                print(' ' + str(p_count), end =" ")
            r=RR[x]
            c=CC[x]
            V = B[r,c]
            if V>0:
                #Left and Right cells are zeros
                if B[r,c+1]==0 and B[r,c-1]==0:
                    #The bottom cells are all zeros as well, but there is a cell directly above that is legit
                    if (B[r+1,c-1]+B[r+1,c]+B[r+1,c+1])==0 and B[r-1,c]>0:
                        B[r,c] = 0
                        n=n+1
                    #The top cells are all zeros as well, but there is a cell directly below that is legit
                    elif (B[r-1,c-1]+B[r-1,c]+B[r-1,c+1])==0 and B[r+1,c]>0:
                        B[r,c] = 0
                        n=n+1
                #top and bottom cells are zeros
                if B[r,c]>0 and B[r+1,c]==0 and B[r-1,c]==0:
                    #All cells on the right are zero, but there is a cell to the left that is legit
                    if (B[r+1,c+1]+B[r,c+1]+B[r-1,c+1])==0 and B[r,c-1]>0:
                        B[r,c] = 0
                        n=n+1
                    elif (B[r+1,c-1]+B[r,c-1]+B[r-1,c-1])==0 and B[r,c+1]>0:
                        B[r,c] = 0
                        n=n+1
        print('\nFirst pass removed ' + str(n) + ' cells')
        
        
        #This pass is to remove all the redundant cells
        n=0
        p_count = 0
        p_percent = (num_nonzero+1)/100.0
        for x in range(num_nonzero):
            if x>=p_count*p_percent:
                p_count = p_count + 1
                print(' ' + str(p_count), end =" ")
            r=RR[x]
            c=CC[x]
            V = B[r,c]
            if V>0:
                if B[r+1,c]==V and (B[r+1,c+1]==V or B[r+1,c-1]==V):
                    if sum(B[r+1,c-1:c+2])==0:
                        B[r+1,c] = 0
                        n=n+1
                elif B[r-1,c]==V and (B[r-1,c+1]==V or B[r-1,c-1]==V):
                    if sum(B[r-1,c-1:c+2])==0:
                        B[r-1,c] = 0
                        n=n+1
                elif B[r,c+1]==V and (B[r+1,c+1]==V or B[r-1,c+1]==V):
                    if sum(B[r-1:r+1,c+2])==0:
                        B[r,c+1] = 0
                        n=n+1
                elif B[r,c-1]==V and (B[r+1,c-1]==V or B[r-1,c-1]==V):
                    if sum(B[r-1:r+1,c-2])==0:
                            B[r,c-1] = 0
                            n=n+1
        print('\nSecond pass removed ' + str(n) + ' redundant cells')
    
    print('Writing Output File ' + STRM_File_Clean)
    Write_Output_Raster(STRM_File_Clean, B[1:nrows+1,1:ncols+1], ncols, nrows, dem_geotransform, dem_projection, "GTiff", gdal.GDT_Int32)
    #return B[1:nrows+1,1:ncols+1], ncols, nrows, cellsize, yll, yur, xll, xur
    return


def Check_and_Change_Coordinate_Systems(DEM_File, LandCoverFile):

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
        warp = gdal.Warp(output_raster,input_raster,dstSRS=dem_projection)
        warp = None # Closes the files
        input_raster = None

    # delete the old LAND raster, if it was replaced and change the name
    if current_crs != dem_projection:
        os.remove(LandCoverFile)
        LandCoverFile = LandCoverFile_Update

    return (LandCoverFile)

def find_nearest_idx(point, tree, gdf):
    """Find the nearest index and corresponding data for a given point using spatial indexing."""
    nearest_idx = tree.nearest(point)  # Directly query with the geometry
    return nearest_idx, gdf.iloc[nearest_idx]

# Function to get a point strictly on the network
def get_point_on_stream(line, target_distance):
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

def walk_stream_for_point(line, target_distance):
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

def find_stream_cells_at_increments_above_and_below_dam(CurveParam_File, VDT_File, dam_csv, dam_id_field, dam_id, Dam_StrmShp, dam_reanalysis_flowfile, STRM_Raster_File, known_baseflow, number_of_cross_sections=4):
    """
    Finds a location on the stream network that is downstream of the dam at specified increments and saves it as a shapefile.

    Parameters:
    - CurveParam_File: str, path to CurveParam file
    - VDT_File: str, path to VDT file
    - dam_csv: str, path to the CSV containing dam locations
    - dam_id_field: str, name of the column identifying dams
    - dam_id: int/str, ID of the dam of interest
    - Dam_StrmShp: str, path to the stream shapefile
    - dam_reanalysis_flowfile: str, path to the flow reanalysis file
    - output_shapefile: str, path to save the output shapefile
    - number_of_cross_sections: int, number of distances at which to extract points

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
    Dam_StrmShp_gdf = gpd.read_file(Dam_StrmShp)
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
    
    for _, row in Dam_StrmShp_gdf.iterrows():
        link_id = row['LINKNO']
        ds_link_id = row['DSLINKNO']
        geometry = row.geometry

        if ds_link_id > 0:  # Ignore terminal reaches
            G.add_edge(link_id, ds_link_id, geometry=geometry, weight=geometry.length)

    # **2. Find the Closest Stream to the Dam**
    Dam_StrmShp_gdf['distance'] = Dam_StrmShp_gdf.distance(dam_point)
    closest_stream = Dam_StrmShp_gdf.loc[Dam_StrmShp_gdf['distance'].idxmin()]
    start_link = closest_stream['LINKNO']

    # Filter the top-width data to the closest stream and find the median
    merged_df = merged_df[merged_df['COMID']==start_link]
    if known_baseflow is None:
        tw_median = merged_df.groupby('COMID')['tw_rp2'].median()
    else:
        tw_median = merged_df.groupby('COMID')['tw_known_baseflow'].median()
    

    # **3. Calculate the Top Width (tw) for the Closest Stream and find the locations of interest on the stream network**
    # Select tw_median based on closest stream 'COMID'
    tw = tw_median.get(start_link, 100)  # Default to 50m if not found

    # use a minimum of 100 meters for the top width
    if tw < 100:
        tw = 100


    # Get the dam’s intersection point on the stream
    current_point = nearest_points(closest_stream.geometry, dam_point)[0]
    # Assume current_link is the stream segment containing current_point (the dam intersection)
    current_link = start_link

    dam_ids = []
    link_nos = []
    points_of_interest = []

    # Loop for each downstream cross-section
    for i in range(1, number_of_cross_sections + 1):
        # We want to move exactly tw meters from the current point
        remaining_distance_to_travel = tw

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
            seg_geom = G[current_link][next_link]['geometry']
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
            .to_crs("EPSG:4326")
            .geometry.iloc[0]
        )
        dam_ids.append(dam_id)
        link_nos.append(current_link)
        points_of_interest.append(downstream_point)

    # downstream_points now contains cross-sections spaced exactly tw meters apart along the stream.

 
    # now let's find the upstream point that is 1/8 of the top width upstream of the dam location
    # Get the dam’s intersection point on the stream
    current_point = nearest_points(closest_stream.geometry, dam_point)[0]
    # Assume current_link is the stream segment containing current_point (the dam intersection)
    current_link = start_link
    # calculate 1/8 tw to look upstream of the dam location to use upstream water surface elevation and known discharge to estimate dam height
    tw_upstream = tw / 8
    remaining_distance_to_travel = tw_upstream

    print(f"Calculating point {remaining_distance_to_travel} meters upstream of the dam.")

    while remaining_distance_to_travel > 0:
        # Get geometry of the current link segment (self-contained)
        # Assumes you are storing geometry from current_link to one of its neighbors — adjust if needed
        # You may need to maintain a mapping or use a convention to fetch a self-link geometry if it's not available in edges
        if 'geometry' in G.nodes[current_link]:
            seg_geom = G.nodes[current_link]['geometry']
        else:
            # Fallback: get one of the out_edges or in_edges that contains this geometry
            # For now, we use a simple assumption: pull from an out_edge or in_edge
            candidate_edges = list(G.out_edges(current_link, data=True)) + list(G.in_edges(current_link, data=True))
            for u, v, data in candidate_edges:
                if 'geometry' in data:
                    seg_geom = data['geometry']
                    break
            else:
                raise ValueError(f"No geometry found for link {current_link}")

        # Prepare geometry
        if seg_geom.geom_type.startswith("Multi"):
            seg_coords = list(linemerge(seg_geom).coords)
        else:
            seg_coords = list(seg_geom.coords)

        # Ensure the geometry direction ends at current_point (reverse if needed)
        if not Point(seg_coords[-1]).equals_exact(current_point, tolerance=0.2):
            seg_coords.reverse()

        seg_line = LineString(seg_coords)

        # STEP 1: Measure how far upstream we can go in this segment
        proj_distance = seg_line.project(current_point)
        if proj_distance >= remaining_distance_to_travel:
            # We have enough upstream length within current_link itself
            current_point = seg_line.interpolate(proj_distance - remaining_distance_to_travel)
            remaining_distance_to_travel = 0
        else:
            # Move to the start of the segment
            remaining_distance_to_travel -= proj_distance
            current_point = Point(seg_coords[0])  # Start of the upstream segment

            # Now move to upstream link
            upstream_edges = list(G.in_edges(current_link, data=True))
            if not upstream_edges:
                raise ValueError(f"Ran out of stream while looking upstream from {current_link}")
            
            # Pick first upstream link — modify if smarter selection needed
            next_link = upstream_edges[0][0]
            current_link = next_link

        # At this point, current_point has moved exactly tw/4 meters upstream from its previous location.
        # Record the result (convert to EPSG:4326 if needed).
        upstream_point = (
            gpd.GeoSeries([current_point], crs=projected_crs)
            .to_crs("EPSG:4326")
            .geometry.iloc[0]
        )
        dam_ids.append(dam_id)
        link_nos.append(current_link)
        points_of_interest.append(upstream_point)


    # **Save the Downstream Points as a Shapefile**
    downstream_gdf = gpd.GeoDataFrame({'dam_id': dam_ids, 'LINKNO': link_nos, 'geometry': points_of_interest}, crs="EPSG:4326")

    # **4. Find Nearest VDT and Curve Data Points to the points of interest**
    # Compute lat/lon for curve data
    (minx, miny, maxx, maxy, dx, dy, _, _, _, _) = Get_Raster_Details(STRM_Raster_File)
    cellsize_x, cellsize_y = abs(float(dx)), abs(float(dy))
    lat_base, lon_base = float(maxy) - 0.5 * cellsize_y, float(minx) + 0.5 * cellsize_x

    curve_data_df['Lat'] = lat_base - curve_data_df['Row'] * cellsize_y
    curve_data_df['Lon'] = lon_base + curve_data_df['Col'] * cellsize_x
    vdt_df['Lat'] = lat_base - vdt_df['Row'] * cellsize_y
    vdt_df['Lon'] = lon_base + vdt_df['Col'] * cellsize_x

    curve_data_gdf = gpd.GeoDataFrame(curve_data_df, geometry=gpd.points_from_xy(curve_data_df['Lon'], curve_data_df['Lat']), crs="EPSG:4269")
    vdt_gdf = gpd.GeoDataFrame(vdt_df, geometry=gpd.points_from_xy(vdt_df['Lon'], vdt_df['Lat']), crs="EPSG:4269")

    # Convert all to projected CRS
    curve_data_gdf, vdt_gdf, downstream_gdf = (gdf.to_crs(projected_crs) for gdf in [curve_data_gdf, vdt_gdf, downstream_gdf])
    
    vdt_gdfs = []
    curve_data_gdfs = []
    # Extract target point from the GeoDataFrame
    for i in downstream_gdf.index:
        # For example, let's use the first point in the GeoDataFrame as the target point
        target_point = downstream_gdf.geometry.iloc[i]
        # Calculate distances to the target point
        curve_data_gdf['distance'] = curve_data_gdf.geometry.apply(lambda x: target_point.distance(x))
        vdt_gdf['distance'] = vdt_gdf.geometry.apply(lambda x: target_point.distance(x))

        # Find the nearest VDT and curve point
        min_distance_curve_data_gdf = curve_data_gdf['distance'].min()
        min_distance_vdt_gdf = vdt_gdf['distance'].min()

        # Filter the GeoDataFrames to only include the nearest points
        nearest_curves_data_gdf = curve_data_gdf[(curve_data_gdf['distance']==min_distance_curve_data_gdf)]
        nearest_vdt_gdf = vdt_gdf[(vdt_gdf['distance']==min_distance_vdt_gdf)]

        vdt_gdfs.append(nearest_vdt_gdf)
        curve_data_gdfs.append(nearest_curves_data_gdf)
    
    # combine the VDT gdfs and curve data gdfs into one a piece
    vdt_gdf = pd.concat(vdt_gdfs)
    curve_data_gdf = pd.concat(curve_data_gdfs)

    # Dropping the 'distance' column
    vdt_gdf = vdt_gdf.drop(columns=['distance'])
    curve_data_gdf = curve_data_gdf.drop(columns=['distance'])
    
    return downstream_gdf, vdt_gdf, curve_data_gdf


def Dam_Assessment(DEM_Folder, DEM, watershed, ESA_LC_Folder, STRM_Folder, LAND_Folder, FLOW_Folder, 
                     VDT_Folder, ARC_Folder, BathyFileFolder, XS_Folder, ManningN, bathy_use_banks,
                     find_banks_based_on_landcover, create_reach_average_curve_file,
                     dam_csv, dam_id_field, dam_id, known_baseflow, known_channel_forming_discharge,
                     StrmShp_gdf=None):

    if DEM.endswith(".tif") or DEM.endswith(".img"):
        DEM_Name = DEM
        FileName = DEM_Name.replace('.tif','')
        FileName = FileName.replace('.img','')
        DEM_File = os.path.join(DEM_Folder, DEM_Name)
        
        #Input Dataset
        ARC_FileName_Bathy = os.path.join(ARC_Folder, 'ARC_Input_' + str(dam_id) + '.txt')
        
        #Datasets to be Created
        Dam_StrmShp = os.path.join(STRM_Folder, f"{str(dam_id)}_StrmShp.shp")
        Dam_Reanalsyis_FlowFile = os.path.join(FLOW_Folder,f"{str(dam_id)}_Reanalysis.csv")


        STRM_File = os.path.join(STRM_Folder, str(dam_id) + '_STRM_Raster.tif')
        STRM_File_Clean = STRM_File.replace('.tif','_Clean.tif')
        LAND_File = os.path.join(LAND_Folder, str(dam_id) + '_LAND_Raster.tif')
        
        VDT_File = os.path.join(VDT_Folder, str(dam_id) + '_VDT_Database.txt')
        Curve_File = os.path.join(VDT_Folder, str(dam_id) + '_CurveFile.csv')

        # these are the files that will be created by the code
        Local_VDT_File = os.path.join(VDT_Folder, str(dam_id) + '_Local_VDT_Database.shp')
        Local_Curve_File = os.path.join(VDT_Folder, str(dam_id) + '_Local_CurveFile.shp')
        
        ARC_BathyFile = os.path.join(BathyFileFolder, str(dam_id) + '_ARC_Bathy.tif')

        XS_Out_File = os.path.join(XS_Folder, str(dam_id) + '_XS_Out.txt')
        
        #Download and Process Land Cover Data
        LandCoverFile = ''
        if not os.path.exists(LAND_File):
            (lon_1, lat_1, lon_2, lat_2, dx, dy, ncols, nrows, geoTransform, Rast_Projection) = Get_Raster_Details(DEM_File)
            
            # Get geometry in original projection
            geom = ESA.Get_Polygon_Geometry(lon_1, lat_1, lon_2, lat_2)

            # Check if raster projection is WGS 84
            raster_crs = CRS.from_wkt(Rast_Projection)
            wgs84_crs = CRS.from_epsg(4326)

            if raster_crs != wgs84_crs:
                # Convert geometry to WGS 84
                transformer = Transformer.from_crs(raster_crs, wgs84_crs, always_xy=True)
                geom = transform(transformer.transform, geom)

            LandCoverFile = ESA.Download_ESA_WorldLandCover(ESA_LC_Folder, geom, 2021)

        # This function sets-up the Input files for ARC and FloodSpreader
        # It also does some of the geospatial processing
        (ARC_FileName_Bathy, DEM_Reanalsyis_FlowFile, Dam_StrmShp) = Process_Geospatial_Data(ARC_Folder,  
                                                                                            ARC_FileName_Bathy, 
                                                                                            DEM_File, LandCoverFile, STRM_File, 
                                                                                            STRM_File_Clean, LAND_File, 
                                                                                            BathyFileFolder, FLOW_Folder, XS_Folder, ManningN,
                                                                                            VDT_File, Curve_File, ARC_BathyFile, Dam_StrmShp, 
                                                                                            Dam_Reanalsyis_FlowFile, XS_Out_File, bathy_use_banks,
                                                                                            find_banks_based_on_landcover, create_reach_average_curve_file,
                                                                                            dam_csv, dam_id_field, dam_id, known_baseflow, known_channel_forming_discharge,
                                                                                            StrmShp_gdf)  

        # read in the reanalysis streamflow and break the code if the dataframe is empty or if the streamflow is all 0
        DEM_Reanalsyis_FlowFile_df = pd.read_csv(DEM_Reanalsyis_FlowFile)
        if DEM_Reanalsyis_FlowFile_df.empty is True or DEM_Reanalsyis_FlowFile_df['qout_max'].mean() <= 0 or len(DEM_Reanalsyis_FlowFile_df.index)==0:
            print(f"Dam_Assessment: Results for {DEM} are not possible because we don't have streamflow estimates...")
            return

        
        # Create our Curve and VDT Database data
        if os.path.exists(ARC_BathyFile) == False or os.path.exists(VDT_File) == False or os.path.exists(Curve_File) == False:
            print('Dam_Assessment: Cannot find bathy file, so creating ' + ARC_BathyFile)
            arc = Arc(ARC_FileName_Bathy)
            arc.run() # Runs ARC
        
        # Now we need to use the Dam_StrmShp and VDT data to find the stream cells at distance increments below the dam
        downstream_gdf, vdt_gdf, curve_data_gdf = find_stream_cells_at_increments_above_and_below_dam(Curve_File, VDT_File, dam_csv, dam_id_field, dam_id,
                                                                                                      Dam_StrmShp, Dam_Reanalsyis_FlowFile,
                                                                                                      STRM_File_Clean, known_baseflow)
        
        # output the results to shapefiles
        vdt_gdf.to_file(Local_VDT_File)
        curve_data_gdf.to_file(Local_Curve_File)


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


    #Folder Management
    Output_Dir = dam_dict['output_dir']
    dam = dam_dict['name']
    dam_csv = dam_dict['dam_csv']
    dam_id_field = dam_dict['dam_id_field']
    dam_id = dam_dict['dam_id']
    DEM_Folder = dam_dict['dem_dir']
    ARC_Folder = os.path.join(Output_Dir, str(dam), 'ARC_InputFiles')
    BathyFileFolder = os.path.join(Output_Dir, str(dam), 'Bathymetry')
    STRM_Folder = os.path.join(Output_Dir, str(dam), 'STRM')
    LAND_Folder = os.path.join(Output_Dir, str(dam), 'LAND')
    FLOW_Folder = os.path.join(Output_Dir, str(dam), 'FLOW')
    VDT_Folder = os.path.join(Output_Dir, str(dam), 'VDT')
    ESA_LC_Folder = os.path.join(Output_Dir, str(dam), 'ESA_LC')
    XS_Folder = os.path.join(Output_Dir, str(dam), 'XS')

    #Create Folders
    # Create_Folder(watershed)
    Create_Folder(ESA_LC_Folder)
    Create_Folder(STRM_Folder)
    Create_Folder(LAND_Folder)
    Create_Folder(FLOW_Folder)
    Create_Folder(VDT_Folder)
    Create_Folder(ARC_Folder)
    Create_Folder(BathyFileFolder)
    Create_Folder(XS_Folder)
    
    #Datasets that can be good for a large domain
    StrmSHP = dam_dict['flowline']
    ManningN = os.path.join(LAND_Folder, 'AR_Manning_n_MED.txt')

    #Create a Baseline Manning N File
    print('Creating Manning n file: ' + ManningN)
    Create_BaseLine_Manning_n_File_ESA(ManningN)
    
    #This is the list of all the DEM files we will go through
    DEM_List = os.listdir(DEM_Folder)

    # Before we get too far ahead, let's make sure that our DEMs and Flowlines have the same coordinate system
    # we will assume that all DEMs in the DEM list have the same coordinate system
    if process_stream_network is True:
        print('Reading in stream file: ' + StrmSHP)
        if StrmSHP.endswith(".gdb"):
            # Specify the layer you want to access
            layer_name = "geoglowsv2"
            # Read the layer from the geodatabase
            StrmShp_gdf = gpd.read_file(StrmSHP, layer=layer_name)    
        elif StrmSHP.endswith(".shp") or StrmSHP.endswith(".gpkg"):
            # Read the layer from the shapefile
            StrmShp_gdf = gpd.read_file(StrmSHP)
        elif StrmSHP.endswith(".parquet"):
            # Read the layer from the shapefile
            StrmShp_gdf = gpd.read_parquet(StrmSHP)

        # removing any lingering NoneType geometries
        StrmShp_gdf = StrmShp_gdf[~StrmShp_gdf.geometry.isna()]

        print('Converting the coordinate system of the stream file to match the DEM files, if necessary')
        test_dem = next((file for file in DEM_List if file.endswith('.tif')), None)
        test_dem_path = os.path.join(DEM_Folder,test_dem)
        # Load the DEM file and get its CRS using gdal
        dem_dataset = gdal.Open(test_dem_path)
        dem_proj = dem_dataset.GetProjection()  # Get the projection as a WKT string
        dem_spatial_ref = osr.SpatialReference()
        dem_spatial_ref.ImportFromWkt(dem_proj)
        # Get the EPSG code
        dem_spatial_ref.AutoIdentifyEPSG()
        dem_epsg_code = dem_spatial_ref.GetAuthorityCode(None)  # This extracts the EPSG code as a string
        # Check if the CRS of the shapefile matches the DEM's CRS
        if int(str(StrmShp_gdf.crs)[5:]) != int(dem_epsg_code):
            print("DEM and Stream Network have different coordinate systems...")
            print(f"Stream CRS: {str(StrmShp_gdf.crs)[5:]}")
            print(f"DEM CRS: {dem_epsg_code}")
            # Reproject the shapefile to match the DEM's CRS
            StrmShp_gdf = StrmShp_gdf.to_crs(dem_epsg_code)
        dem_dataset = None
        dem_proj = None 
        dem_spatial_ref = None
    elif process_stream_network is False:
        StrmShp_gdf = None   

    #Now go through each DEM dataset
    for DEM in DEM_List:
            
        Dam_Assessment(DEM_Folder, DEM, dam, ESA_LC_Folder, STRM_Folder, LAND_Folder, FLOW_Folder, 
                     VDT_Folder, ARC_Folder, BathyFileFolder, XS_Folder, ManningN, bathy_use_banks,
                     find_banks_based_on_landcover, create_reach_average_curve_file,
                     dam_csv, dam_id_field, dam_id, known_baseflow, known_channel_forming_discharge,
                     StrmShp_gdf)
        

    # delete the ESA_LC_Folder and the data in it
    # Loop through all files in the directory and remove them
    for file in Path(ESA_LC_Folder).glob("*"):
        try:
            if file.is_file():
                # Adjust file permissions before deletion
                if platform.system() == "Windows":
                    os.chmod(file, stat.S_IWRITE)  # Remove read-only attribute on Windows
                else:
                    os.chmod(file, stat.S_IWUSR)   # Give user write permission on Unix systems
                os.remove(file)
                print(f"process_dam: Deleted file: {file}")
        except Exception as e:
            print(f"Error deleting file {file}: {e}")
    if os.path.exists(ESA_LC_Folder):
        # Adjust file permissions before deletion
        if platform.system() == "Windows":
            os.chmod(ESA_LC_Folder, stat.S_IWRITE)  # Remove read-only attribute on Windows
        else:
            os.chmod(ESA_LC_Folder, stat.S_IWUSR)   # Give user write permission on Unix systems
        os.rmdir(ESA_LC_Folder)
        print(f"process_dam: Deleted empty folder: {ESA_LC_Folder}")
    else:
        print(f"process_dam: Folder {ESA_LC_Folder} does not exist.")
    
    return
    
def process_json_input(json_file):
    """Process input from a JSON file."""
    with open(json_file, 'r') as file:
        print('Opening ' + str(file))
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

    args = parser.parse_args()

    if args.command == "json":
        print('Processing ' + str(args.json_file))
        process_json_input(args.json_file)
    elif args.command == "cli":
        process_cli_arguments(args)

if __name__ == "__main__":
    main()








