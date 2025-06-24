# RathCelon
A tool for extracting rating curves below a low-head dam. The current code is explicity set up to work with the [ECMWF GEOGLOWS Streamflow Service](https://geoglows.ecmwf.int/) datasets. 

## How to Set-up RathCelon
1a. Using Anaconda run:

```bash
conda create -n rathcelon_py310 python=3.10 pyarrow geopandas pandas netcdf4 dask fiona s3fs xarray zarr beautifulsoup4 dataretrieval geojson progress tqdm pygeos pillow=9.0.1 rasterio
```
1b. Alternatively, you can create an environment using conda and then get the libraries using pip

```bash
conda create -n rathcelon_py310 python=3.10
conda activate rathcelon_py310
pip install --upgrade pip
pip install pyarrow geopandas pandas netCDF4 dask fiona s3fs xarray zarr beautifulsoup4 dataretrieval geojson progress tqdm pygeos pillow==9.0.1 rasterio
pip install memory_profiler
conda install -c conda-forge gdal
```

2. Activate you conda environment with the command:
```bash
conda activate rathcelon_py310
```

3. Install the [Automated Rating Curve (ARC) Tool](https://github.com/MikeFHS/automated-rating-curve) in your conda environment by downloading the source code, navigating to the local instance on your machine, and typing the command:

```bash
pip install .
```

4. Now download this RathCelon repository and run this command from within RathCelon's local directory:
```bash
pip install .
```

5. Type the followng command into the prompt to ensure everything is working and for additional guidance
```bash
rathcelon -h
```

## Inputs and Runs
RathCelon can be ran from a command prompt in two ways. 

### Running with a JSON
The first is by creating an input JSON file like the one below. Note that if you are on a PC, the `\` in your path need to be changed to '/' in the JSON file:

```json
{
    "dams": [
        {
            "name": "152",
            "dam_csv":"C:/Users/jlgut/OneDrive/Desktop/rathcelon/examples/LHD_lat_long.csv",
            "dam_id_field": "BYU_SITE_ID",
            "dam_id": 152,
            "flowline": "C:/Users/jlgut/Downloads/streams_712.gpkg",
            "dem_dir": "C:/Users/jlgut/OneDrive/Desktop/rathcelon_test_data/DEM_3",
            "bathy_use_banks": false,
            "output_dir": "C:/Users/jlgut/OneDrive/Desktop/rathcelon_test_data/Results",
            "process_stream_network": true,
            "find_banks_based_on_landcover": true,
            "create_reach_average_curve_file": false,
            "known_baseflow": 5.0,
            "known_channel_forming_discharge": 10.0,
            "upstream_elevation_change_threshold": 0.5
        }
    ]
}
```

These inputs are:

    "name": (str) The name of the dam you're analyzing
    "dam_csv": (str) The path to the CSV containing the locations of the dams
    "dam_id_field": (str) The field in the dam_csv that stores the dam unique identifiers
    "dam_id": (int) The ID for this dam in the dam_id_field
    "flowline": (str) The path to your flowline dataset that you downloaded from [GEOGLOWS](http://geoglows-v2.s3-website-us-west-2.amazonaws.com/#streams/)
    "dem_dir": (str) The path to the folder containing your digital elevation models (DEMs) for the dam
    "bathy_use_banks": (true/false) Boolean indicating the option in ARC described [here](https://github.com/MikeFHS/automated-rating-curve/wiki/Estimating-Bathymetry-in-ARC#bank-or-water-surface-elevation)
    "output_dir": (str) Path to where the outputs of your RathCelon output will be stored
    "process_stream_network": (true/false) If true, RathCelon will take the flowline option and determine the local instance it will use to do its analysis. If false, the network has already been processed to find a local instance
    "find_banks_based_on_landcover": (true/false) If true, the banks of the stream will be found using [this logic](https://github.com/MikeFHS/automated-rating-curve/wiki/Estimating-Bathymetry-in-ARC#discovering-the-waterfront). If false, ARC will look to find the banks by assuming that the stream is a flat surface in the DEM. 
    "create_reach_average_curve_file": (true/false) If true, RathCelon with direct ARC to generate a curve file that is the same for all stream cells on a stream reach. If false, RathCelon will direct ARC to generate a curve file that has a local curve for each stream cell. 
    "known_baseflow": (float) Optional input expressed as a value in cubic meters per second. This is used in conjunction with "bathy_use_banks" = False. This will be substituted into the ARC's processing to estimate a bathymetry in each cross section. This value will also be used to calculate the top width used for distance calculations. 
    "known_channel_forming_discharge": (float) Optional input expressed as a value in cubic meters per second. This is used in conjunction with "bathy_use_banks" = True. This will be substituted into the ARC's processing to estimate a bathymetry in each cross section. 
    "upstream_elevation_change_threshold": (float) Required. The changed in upstream elevation used to identify the cross section upstream of the dam. The default is 1.0.

To run with the JSON, issue the following command in you command line window:

```bash
rathcelon json "path\to\your\input.json"
```

### Running from command line
You can also run RathCelon from the command line without the JSON.

This can be done by issueing the following command:

```bash
rathcelon cli name dam_csv dam_id_field dam_id flowline dem_dir output_dir --bathy_use_banks --process_stream_network --find_banks_based_on_landcover --create_reach_average_curve_file --known_baseflow 5.0 --known_channel_forming_discharge 10.0 --upstream_elevation_change_threshold 0.1
```

Issuing the commands `--bathy_use_banks`, `--process_stream_network`, `--find_banks_based_on_landcover`, and `--create_reach_average_curve_file` indicates that those options are set to True. 

## Outputs
In the output_dir you specify in your inputs, you will find a folder that is the name you specify in your inputs. In the name folder, will be several subfolders. 

The "VDT" folder contains the shapefiles {name}_Local_VDT_Database.shp and {name}_Local_CurveFile.shp. 

These will be the 4 nearest VDT database and curvefile stream cells that are 1, 2, 3, and 4x the top width distance of the stream downstream of the dam location and 1 stream cell that is located at the stream cell where the upstream elevation change is greater than or equal to your `upstream_elevation_change_threshold` value.

If you set "known_baseflow" in your inputs, this will be used by ARC to estimate bathymetry (if "bathy_use_banks" = False) and RathCelon to estimate the top width used in the distance calculations. Otherwise, top width was estimated using the 2yr discharge in the "output_dir/FLOW/{name}_Reanalysis.csv" file.

In either case, RathCelon finds the median top width for the stream reach on which the dam is located and uses this value to calculate the distance. 

If the original 1x top-width distance is <100 meters, RathCelon will default to using 100 meters for the top width distance. 

The "XS" folder contains the cross-section information for all stream cells. The shapefiles {name}_Local_XS_Lines.shp will contain the cross-section locations and information for the stream cells in the {name}_Local_CurveFile.shp. 

For a breakdown of the attributes in the curve file, see [this wiki](https://github.com/MikeFHS/automated-rating-curve/wiki/Running-ARC-and-Looking-at-ARC-Outputs). 



