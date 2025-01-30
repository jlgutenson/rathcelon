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
import numpy as np

# local imports
from . import streamflow_processing as HistFlows
from . import esa_download_processing as ESA


def Process_FloodForecasting_Geospatial_Data(ARC_Folder, ARC_FileName_Initial, 
                                            ARC_FileName_Bathy, ARC_FileName_FloodForecast, 
                                            DEM_File, LandCoverFile, 
                                            VDT_Test_File, STRM_File, 
                                            STRM_File_Clean, LAND_File, 
                                            BathyFileFolder, FLOW_Folder, ManningN, 
                                            VDT_File, Curve_File, ARC_BathyFile, FS_BathyFile, DEM_StrmShp, 
                                            DEM_Reanalsyis_FlowFile, bathy_use_banks, 
                                            find_banks_based_on_landcover, create_reach_average_curve_file,
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

    # now we need to figure out if our DEM_StrmShp and DEM_Reanalysis_Flowfile exists and if not, create it
    if os.path.isfile(DEM_StrmShp) and os.path.isfile(DEM_Reanalsyis_FlowFile):
        print(DEM_StrmShp + ' Already Exists')
        print(DEM_Reanalsyis_FlowFile + ' Already Exists')
        DEM_StrmShp_gdf = gpd.read_file(DEM_StrmShp)
        rivids = DEM_StrmShp_gdf['LINKNO'].values
    elif StrmShp_gdf is not None and os.path.isfile(DEM_StrmShp) is False and os.path.isfile(DEM_Reanalsyis_FlowFile) is False:
        (DEM_Reanalsyis_FlowFile, DEM_StrmShp, rivids, DEM_StrmShp_gdf) = HistFlows.Process_and_Write_Retrospective_Data_for_DEM_Tile(StrmShp_gdf, 'LINKNO', DEM_File, DEM_Reanalsyis_FlowFile, DEM_StrmShp)

    #Create Stream Raster
    if os.path.isfile(STRM_File):
        print(STRM_File + ' Already Exists')
    else:
        print('Creating ' + STRM_File)
        Create_AR_StrmRaster(DEM_StrmShp, STRM_File, outputBounds, minx, miny, maxx, maxy, dx, dy, ncols, nrows, 'LINKNO')
    
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
    
    #Create a Starting AutoRoute Input File
    print('Creating ARC Input File: ' + ARC_FileName_Initial)
    #Create the Initial Flow


    #Create the Bathy Input File
    print('Creating ARC Input File: ' + ARC_FileName_Bathy)

    if bathy_use_banks is False:
        COMID_Param = 'COMID'
        Q_BF_Param = 'qout_median'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, COMID_Param, Q_BF_Param, Q_Param, STRM_File_Clean, LAND_File, DEM_Reanalsyis_FlowFile, VDT_File, Curve_File, ManningN, ARC_BathyFile, FS_BathyFile, VDT_Test_File, DEM_StrmShp, bathy_use_banks, find_banks_based_on_landcover, create_reach_average_curve_file)
    elif bathy_use_banks is True:
        COMID_Param = 'COMID'
        Q_BF_Param = 'rp2'
        Q_Param = 'rp100_premium'
        Create_ARC_Model_Input_File_Bathy(ARC_FileName_Bathy, DEM_File, COMID_Param, Q_BF_Param, Q_Param, STRM_File_Clean, LAND_File, DEM_Reanalsyis_FlowFile, VDT_File, Curve_File, ManningN, ARC_BathyFile, FS_BathyFile, VDT_Test_File, DEM_StrmShp, bathy_use_banks, find_banks_based_on_landcover, create_reach_average_curve_file)

    
    return ARC_FileName_Initial, ARC_FileName_Bathy, ARC_FileName_FloodForecast, Forecast_Flood_Map, DEM_Reanalsyis_FlowFile, ForecastFlowFile, DEM_StrmShp, forecastdate

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
                                      Curve_File, ManningN, ARC_BathyFile, FS_BathyFile, VDT_Test_File, 
                                      DEM_StrmShp, bathy_use_banks, find_banks_based_on_landcover, 
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
    out_file.write('\n' + 'X_Section_Dist	5000.0')
    out_file.write('\n' + 'Degree_Manip	6.1')
    out_file.write('\n' + 'Degree_Interval	1.5')
    out_file.write('\n' + 'Low_Spot_Range	2')
    out_file.write('\n' + 'Str_Limit_Val	1')
    out_file.write('\n' + 'Gen_Dir_Dist	10')
    out_file.write('\n' + 'Gen_Slope_Dist	10')
    
    out_file.write('\n\n#VDT_Output_File_and_CurveFile')
    out_file.write('\n' + 'VDT_Database_NumIterations	30')
    out_file.write('\n' + 'Print_VDT_Database' + '\t' + VDT_File.replace('.txt', '_Bathy.txt'))
    out_file.write('\n' + 'Print_Curve_File' + '\t' + Curve_File.replace('.csv', '_Bathy.csv'))
    out_file.write('\n' + 'Reach_Average_Curve_File' + '\t' + f'{create_reach_average_curve_file}')

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

def DEM_Forecast(DEM_Folder, DEM, watershed, ESA_LC_Folder, STRM_Folder, LAND_Folder, FLOW_Folder, 
                     VDT_Folder, ARC_Folder, BathyFileFolder, ManningN, bathy_use_banks, 
                     find_banks_based_on_landcover, create_reach_average_curve_file,
                     StrmShp_gdf=None):

    if DEM.endswith(".tif") or DEM.endswith(".img"):
        DEM_Name = DEM
        FileName = DEM_Name.replace('.tif','')
        FileName = FileName.replace('.img','')
        DEM_File = os.path.join(DEM_Folder, DEM_Name)
        
        #Input Dataset
        ARC_FileName_Initial = os.path.join(ARC_Folder, 'ARC_Input_' + FileName + '_InitialFlood.txt')
        ARC_FileName_Bathy = os.path.join(ARC_Folder, 'ARC_Input_' + FileName + '_Bathy.txt')
        ARC_FileName_FloodForecast = os.path.join(ARC_Folder, 'ARC_Input_' + FileName + '_FloodForecast.txt')
        VDT_Test_File = os.path.join(VDT_Folder, FileName + '_VDT_FS.csv')
        
        #Datasets to be Created
        DEM_StrmShp = os.path.join(STRM_Folder, f"{FileName}_StrmShp.shp")
        DEM_Reanalsyis_FlowFile = os.path.join(FLOW_Folder,f"{FileName}_Reanalysis.csv")


        STRM_File = os.path.join(STRM_Folder, FileName + '_STRM_Raster.tif')
        STRM_File_Clean = STRM_File.replace('.tif','_Clean.tif')
        LAND_File = os.path.join(LAND_Folder, FileName + '_LAND_Raster.tif')
        
        VDT_File = os.path.join(VDT_Folder, FileName + '_VDT_Database.txt')
        Curve_File = os.path.join(VDT_Folder, FileName + '_CurveFile.csv')
        
        ARC_BathyFile = os.path.join(BathyFileFolder, FileName + '_ARC_Bathy.tif')
        FS_BathyFile = os.path.join(BathyFileFolder, FileName + '_FS_Bathy.tif')    

        
        #Download and Process Land Cover Data
        LandCoverFile = ''
        if not os.path.exists(LAND_File):
            (lon_1, lat_1, lon_2, lat_2, dx, dy, ncols, nrows, geoTransform, Rast_Projection) = Get_Raster_Details(DEM_File)
            geom = ESA.Get_Polygon_Geometry(lon_1, lat_1, lon_2, lat_2)
            LandCoverFile = ESA.Download_ESA_WorldLandCover(ESA_LC_Folder, geom, 2021)

        # This function sets-up the Input files for ARC and FloodSpreader
        # It also does some of the geospatial processing
        (ARC_FileName_Initial, ARC_FileName_Bathy, ARC_FileName_FloodForecast, Forecast_Flood_Map, DEM_Reanalsyis_FlowFile, ForecastFlowFile, DEM_StrmShp, forecastdate) = Process_FloodForecasting_Geospatial_Data(ARC_Folder, ARC_FileName_Initial, 
                                                                                                                                                                           ARC_FileName_Bathy, ARC_FileName_FloodForecast, 
                                                                                                                                                                           DEM_File, LandCoverFile, 
                                                                                                                                                                           VDT_Test_File, STRM_File, 
                                                                                                                                                                           STRM_File_Clean, LAND_File, 
                                                                                                                                                                           BathyFileFolder, FLOW_Folder, ManningN, 
                                                                                                                                                                           VDT_File, Curve_File, ARC_BathyFile, FS_BathyFile, DEM_StrmShp, 
                                                                                                                                                                           DEM_Reanalsyis_FlowFile, bathy_use_banks, 
                                                                                                                                                                           find_banks_based_on_landcover, create_reach_average_curve_file,
                                                                                                                                                                           StrmShp_gdf)  

        # read in the reanalysis streamflow and break the code if the dataframe is empty or if the streamflow is all 0
        DEM_Reanalsyis_FlowFile_df = pd.read_csv(DEM_Reanalsyis_FlowFile)
        if DEM_Reanalsyis_FlowFile_df.empty is True or DEM_Reanalsyis_FlowFile_df['qout_max'].mean() <= 0 or len(DEM_Reanalsyis_FlowFile_df.index)==0:
            print(f"Results for {DEM} are not possible because we don't have streamflow estimates...")
            return

        
        # Create our Curve and VDT Database data
        if os.path.exists(FS_BathyFile) is False:
            print('Cannot find bathy file, so creating ' + FS_BathyFile)
            if os.path.exists(ARC_BathyFile) is False:
                arc = Arc(ARC_FileName_Bathy)
                arc.run() # Runs ARC

    return

def process_dem(watershed_dict):
    
    # set if the system will use banks of water surface elevation to estimate bathymetry
    bathy_use_banks = watershed_dict['bathy_use_banks']

    # if you don't have the stream network preprocessed, do so
    process_stream_network = watershed_dict['process_stream_network']

    # set if you want to use the landcover data to find the banks of the river, instead of the flat water surface elevation in the DEM
    find_banks_based_on_landcover = watershed_dict['find_banks_based_on_landcover']

    # let's tell ARC whether we want the curvefile parameters to be the same for each reach or vary by stream cell
    create_reach_average_curve_file = watershed_dict['create_reach_average_curve_file']


    #Folder Management
    Output_Dir = watershed_dict['output_dir']
    watershed = watershed_dict['name']
    DEM_Folder = watershed_dict['dem_dir']
    ARC_Folder = os.path.join(Output_Dir, watershed, 'ARC_InputFiles')
    BathyFileFolder = os.path.join(Output_Dir, watershed, 'Bathymetry')
    STRM_Folder = os.path.join(Output_Dir, watershed, 'STRM')
    LAND_Folder = os.path.join(Output_Dir, watershed, 'LAND')
    FLOW_Folder = os.path.join(Output_Dir, watershed, 'FLOW')
    VDT_Folder = os.path.join(Output_Dir, watershed, 'VDT')
    ESA_LC_Folder = os.path.join(Output_Dir, watershed, 'ESA_LC')
    
    #Create Folders
    # Create_Folder(watershed)
    Create_Folder(ESA_LC_Folder)
    Create_Folder(STRM_Folder)
    Create_Folder(LAND_Folder)
    Create_Folder(FLOW_Folder)
    Create_Folder(VDT_Folder)
    Create_Folder(ARC_Folder)
    Create_Folder(BathyFileFolder)
    
    #Datasets that can be good for a large domain
    StrmSHP = watershed_dict['flowline']
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
        dem_crs = dem_spatial_ref.ExportToProj4()  # Export CRS to a Proj4 string (or other formats if needed)
        # Check if the CRS of the shapefile matches the DEM's CRS
        if StrmShp_gdf.crs != dem_crs:
            # Reproject the shapefile to match the DEM's CRS
            StrmShp_gdf = StrmShp_gdf.to_crs(dem_crs)
        dem_dataset = None
        dem_proj = None 
        dem_spatial_ref = None
        dem_crs = None 
    elif process_stream_network is False:
        StrmShp_gdf = None   

    #Now go through each DEM dataset
    for DEM in DEM_List:
            
        DEM_Forecast(DEM_Folder, DEM, watershed, ESA_LC_Folder, STRM_Folder, LAND_Folder, FLOW_Folder, 
                     VDT_Folder, ARC_Folder, BathyFileFolder, ManningN, bathy_use_banks, 
                     find_banks_based_on_landcover, create_reach_average_curve_file,
                     StrmShp_gdf)
    
    return
    
def process_json_input(json_file):
    """Process input from a JSON file."""
    with open(json_file, 'r') as file:
        print('Opening ' + str(file))
        data = json.load(file)
        print(data)
    
    watersheds = data.get("watersheds", [])
    for watershed in watersheds:
        watershed_name = watershed.get("name")
        flowline = os.path.normpath(watershed.get("flowline"))
        dem_dir = os.path.normpath(watershed.get("dem_dir"))
        output_dir = os.path.normpath(watershed.get("output_dir"))

        watershed_dict = {
            "name": watershed_name,
            "flowline": flowline,
            "dem_dir": dem_dir,
            "output_dir": output_dir,
            "bathy_use_banks": watershed.get("bathy_use_banks", False),
            "flood_waterlc_and_strm_cells":watershed.get("flood_waterlc_and_strm_cells", False),
            "process_stream_network": watershed.get("process_stream_network", False),
            "find_banks_based_on_landcover": watershed.get("find_banks_based_on_landcover", True),
            "create_reach_average_curve_file": watershed.get("create_reach_average_curve_file", False)
        }

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        print(f"Processing watershed: {watershed_name} with parameters: {watershed_dict}")

        # Call your existing processing logic here
        process_dem(watershed_dict)

def normalize_path(path):
    return os.path.normpath(path)

def process_cli_arguments(args):
    """Process input from CLI arguments."""
    output_dir = args.output_dir
    watershed_name = args.watershed
    watershed_dict = {
        "name": watershed_name,
        "flowline": normalize_path(args.flowline),
        "dem_dir": normalize_path(args.dem_dir),
        "bathy_use_banks": args.bathy_use_banks,
        "output_dir": normalize_path(output_dir),
        "process_stream_network": args.process_stream_network,
        "find_banks_based_on_landcover": args.find_banks_based_on_landcover,
        "create_reach_average_curve_file": args.create_reach_average_curve_file
    }

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing watershed: {args.watershed} with parameters: {watershed_dict}")
    print(f"Results will be saved in: {output_dir}")

    # Call the existing processing logic here
    process_dem(watershed_dict)

def main():
    parser = argparse.ArgumentParser(description="Flood Mapping Script")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for JSON input
    json_parser = subparsers.add_parser("json", help="Process watersheds from a JSON file")
    json_parser.add_argument("json_file", type=str, help="Path to the JSON file")

    # Subcommand for CLI input
    cli_parser = subparsers.add_parser("cli", help="Process watershed parameters via CLI")
    cli_parser.add_argument("watershed", type=str, help="Watershed name")
    cli_parser.add_argument("flowline", type=str, help="Path to the flowline shapefile")
    cli_parser.add_argument("dem_dir", type=str, help="Directory containing DEM files")
    cli_parser.add_argument("output_dir", type=str, help="Directory where results will be saved")
    cli_parser.add_argument("--bathy_use_banks", action="store_true", help="Use bathy banks for processing")
    cli_parser.add_argument("--process_stream_network", action="store_true", help="Clean DEM data before processing")
    cli_parser.add_argument("--find_banks_based_on_landcover", action="store_true", help="Use landcover data for finding banks when estimating bathymetry")
    cli_parser.add_argument("--create_reach_average_curve_file", action="store_true", help="Create a reach average curve file instead of one that varies for each stream cell")

    args = parser.parse_args()

    if args.command == "json":
        print('Processing ' + str(args.json_file))
        process_json_input(args.json_file)
    elif args.command == "cli":
        process_cli_arguments(args)

if __name__ == "__main__":
    main()








