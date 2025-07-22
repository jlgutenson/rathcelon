#Program downloads esa world land cover datasets and also creates a water mask
#https://esa-worldcover.org/en/data-access

import os
import sys
import geopandas as gpd   #conda install --channel conda-forge geopandas

# use requests library to download them
import requests   #conda install anaconda::requests
from tqdm.auto import tqdm  # provides a progressbar     #conda install conda-forge::tqdm
from pathlib import Path    #conda install anaconda::pathlib
from shapely.geometry import LineString, Polygon    #conda install conda-forge::shapely
import numpy as np

# from pyproj import CRS

try:
    import gdal     #conda install conda-forge::gdal
except ImportError:
    from osgeo import gdal


def Geom_Based_On_Country(country, Shapefile_Use):
    ne = gpd.read_file(Shapefile_Use)
    
    # get AOI geometry (Italy in this case)
    geom = ne[ne.NAME == country].iloc[0].geometry
    return geom


def Download_ESA_WorldLandCover(output_folder, geom, year):
    s3_url_prefix = "https://esa-worldcover.s3.eu-central-1.amazonaws.com"
    # load natural earth low res shapefile
    
    # load worldcover grid
    url = f'{s3_url_prefix}/esa_worldcover_grid.geojson'
    grid = gpd.read_file(url, crs="epsg:4326")
    
    # get grid tiles intersecting AOI
    tiles = grid[grid.intersects(geom)]
    print(tiles)
    
    # select version tag, based on the year
    version = {2020: 'v100',
               2021: 'v200'}[year]
    
    lc_list = []
    for tile in tqdm(tiles.ll_tile):
        url = f"{s3_url_prefix}/{version}/{year}/map/ESA_WorldCover_10m_{year}_{version}_{tile}_Map.tif"
        r = requests.get(url, allow_redirects=True)
        out_fn = Path(output_folder) / Path(url).name
        lc_list.append(str(out_fn))
        if os.path.isfile(out_fn):
            print('Already Exists: ' + str(out_fn))
        else:
            with open(out_fn, 'wb') as f:
                f.write(r.content)

    # let's merge the list of tiles together
    LandCoverFile = os.path.join(output_folder,"merged_ESA_LC.tif")
    
    # Merge rasters
    merged_raster = gdal.Warp(LandCoverFile, lc_list, options=gdal.WarpOptions(format='GTiff'))

    # Ensure data is written and file is closed properly
    if merged_raster:
        merged_raster.FlushCache()  # Save changes
        merged_raster = None  # Close dataset

    return LandCoverFile

def Write_Output_Raster(s_output_filename, raster_data, ncols, nrows, dem_geotransform, dem_projection, s_file_format, s_output_type):   
    """
    Creates a raster from the specified inputs using GDAL
       
    Parameters
    ----------
    s_output_filename: str
        The path and file name of the output raster
    raster_data: arr
        An array of data values that will be written to the output raster
    ncols: int
        The number of columns in the output raster
    nrows: int
        The number of rows in the output raster
    dem_geotransform: list
        A GDAL GetGeoTransform list that is passed to the output raster
    dem_projection: str
        A GDAL GetProjectionRef() string that contains the projection reference that is passed to the output raster
    s_file_format
        The string that specifies the type of raster that will be output (e.g., GeoTIFF = "GTiff")
    s_output_type
        The type of value that the output raster will be stored as (e.g., gdal.GDT_Int32)
    Returns
    -------
    None

    """
    o_driver = gdal.GetDriverByName(s_file_format)  #Typically will be a GeoTIFF "GTiff"
    #o_metadata = o_driver.GetMetadata()
    
    # Construct the file with the appropriate data shape
    o_output_file = o_driver.Create(s_output_filename, xsize=ncols, ysize=nrows, bands=1, eType=s_output_type)
    
    # Set the geotransform
    o_output_file.SetGeoTransform(dem_geotransform)
    
    # Set the spatial reference
    o_output_file.SetProjection(dem_projection)
    
    # Write the data to the file
    o_output_file.GetRasterBand(1).WriteArray(raster_data)
    
    # Once we're done, close properly the dataset
    o_output_file = None

    return

def Read_Raster_GDAL(InRAST_Name):
    """
    Retrieves the geograhic details of a raster using GDAL in a slightly different way than Get_Raster_Details()

    Parameters
    ----------
    InRAST_Name: str
        The file name and full path to the raster you are analyzing

    Returns
    -------
    RastArray: arr
        A numpy array of the values in the first band of the raster you are analyzing
    ncols: int
        The raster width in pixels
    nrows: int
        The raster height in pixels
    cellsize: float
        The pixel size of the raster longitudinally
    yll: float
        The lowest latitude of the raster
    yur: float
        The latitude of the top left corner of the top pixel of the raster
    xll: float
        The longitude of the top left corner of the top pixel of the raster
    xur: float
        The highest longitude of the raster
    lat: float
        The average of the yur and yll latitude values
    geoTransform: list
        A list of geotransform characteristics for the raster
    Rast_Projection:str
        The projection system reference for the raster
    """
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
    xll = geotransform[0]
    xur = xll + (ncols)*geotransform[1]
    lat = np.fabs((yll+yur)/2.0)
    Rast_Projection = dataset.GetProjectionRef()
    del dataset
    print(f'Spatial Data for Raster File:\n'
          f'\tncols = {ncols}\n'
          f'\tnrows = {nrows}\n'
          f'\tcellsize = {cellsize}\n'
          f'\tyll = {yll}\n'
          f'\tyur = {yur}\n'
          f'\txll = {xll}\n'
          f'\txur = {xur}')

    return RastArray, ncols, nrows, cellsize, yll, yur, xll, xur, lat, geotransform, Rast_Projection


def Get_Raster_Details(DEM_File):
    """
    Retrieves the geographic details of a raster using GDAL in a slightly different way than Read_Raster_GDAL()

    Parameters
    ----------
    DEM_File: str
        The file name and full path to the raster you are analyzing

    Returns
    -------
    minx: float
        The longitude of the top left corner of the top pixel of the raster
    miny: 
        The lowest latitude of the raster
    maxx: 
        The highest latitude of the raster
    maxy:
        The latitude of the top left corner of the top pixel of the raster
    dx: float
        The pixel size of the raster longitudinally
    dy: float
        The pixel size of the raster latitudinally 
    ncols: int
        The raster width in pixels
    nrows: int
        The raster height in pixels
    geoTransform: list
        A list of geotransform characteristics for the raster
    Rast_Projection:str
        The projection system reference for the raster
    """
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
    del data
    return minx, miny, maxx, maxy, dx, dy, ncols, nrows, geoTransform, Rast_Projection


def Create_AR_LandRaster(LandCoverFile, LAND_File, projWin_extents, ncols, nrows):
    """
    Creates a land cover raster that is clipped to a specified extent and cell size
    
   
    Parameters
    ----------
    LandCoverFile: str
        The path and file name of the source National Land Cover Database land-use/land-cover raster
    LAND_File: str
        The path and file name of the output land-use/land-cover dataset 
    projWin_extents: list
        A list of the minimum and maximum extents to which the LAND_File will be clipped, specified as [minimum longitude, maximum latitude, maximum longitude, minimum latitude]
    ncols: int
        The number of columns in the output LAND_File raster
    nrows: int
        The number of rows in the output LAND_File raster
    
    Returns
    -------
    None

    """
    ds = gdal.Open(LandCoverFile)
    ds = gdal.Translate(LAND_File, ds, projWin = projWin_extents, width=ncols, height = nrows)
    ds = None
    return


def Create_Water_Mask(lc_file, waterboundary_file, watervalue):
    (RastArray, ncols, nrows, cellsize, yll, yur, xll, xur, lat, geotransform, Rast_Projection) = Read_Raster_GDAL(lc_file)
    RastArray = np.where(RastArray==watervalue,1,0)   #Streams are identified with zeros
    Write_Output_Raster(waterboundary_file, RastArray, ncols, nrows, geotransform, Rast_Projection, "GTiff", gdal.GDT_Byte)
    return


def Get_Polygon_Geometry(lon_1, lat_1, lon_2, lat_2):
    return Polygon([[min(lon_1,lon_2),min(lat_1,lat_2)], [min(lon_1,lon_2),max(lat_1,lat_2)], [max(lon_1,lon_2),max(lat_1,lat_2)], [max(lon_1,lon_2),min(lat_1,lat_2)]])


if __name__ == "__main__":
    
    
    #Just leave blank if using Option 1 or 2 below
    if len(sys.argv) > 1:
        DEM_File = sys.argv[1]
        print('Input DEM File: ' + DEM_File)
    else:
        DEM_File = 'NED_n39w090_Clipped.tif'
        print('Did not input DEM, going with default: ' + DEM_File)
    
    
    year = 2021  # setting this to 2020 will download the v100 product instead
    output_folder = 'ESA_LC'  # use current directory or set a different one to store downloaded files
    if not os.path.exists(output_folder): 
        os.makedirs(output_folder)
    
    '''
    ###Option 1 - Get Geometry from a Shapefile
    #Get Geometry based on Country and Shapefile
    geom = Geom_Based_On_Country('Cyprus', 'ne_110m_admin_0_countries.shp')
    
    
    ###Option 2 - Get Geometry from Lat/Long Bounding Coordinates
    #Get Geometry based on Latitude and Longitude
    lat_1 = 42.5
    lat_2 = 43.0 
    lon_1 = -106.0 
    lon_2 = -106.5
    d = {'col1': ['name1'], 'geometry': LineString([[lon_1, lat_1], [lon_2, lat_2]])}
    geom = Get_Polygon_Geometry(lon_1, lat_1, lon_2, lat_2)
    '''
    
    
    ###Option 3 - Get Geometry from Raster File
    (lon_1, lat_1, lon_2, lat_2, dx, dy, ncols, nrows, geoTransform, Rast_Projection) = Get_Raster_Details(DEM_File)
    geom = Get_Polygon_Geometry(lon_1, lat_1, lon_2, lat_2)
    
    
    lc_list = Download_ESA_WorldLandCover(output_folder, geom, year)
    
    for lc_file in lc_list:
        lc_file_str = str(lc_file)
        
        if DEM_File != '':
            LAND_File_Clipped = lc_file_str.replace('.tif','_Clipped.tif')
            if os.path.isfile(LAND_File_Clipped):
                print('Already Exists: ' + str(LAND_File_Clipped))
            else:
                print('Creating: ' + str(LAND_File_Clipped))
                Create_AR_LandRaster(lc_file_str, LAND_File_Clipped, [lon_1, lat_2, lon_2, lat_1], ncols, nrows)
            lc_file_str = LAND_File_Clipped
            
        '''
        waterboundary_file = lc_file_str.replace('.tif','_wb.tif')
        if os.path.isfile(waterboundary_file):
            print('Already Exists: ' + str(waterboundary_file))
        else:
            print('Creating ' + str(waterboundary_file))
            Create_Water_Mask(lc_file_str, waterboundary_file, 80)
        '''
    