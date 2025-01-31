# RathCelon
A tool for extracting rating curves below a low-head dam

## How to Set-up
1. Using Anaconda run:

```bash
conda create -n rathcelon_py310 python=3.10 pyarrow geopandas pandas netcdf4 dask fiona s3fs xarray zarr beautifulsoup4 dataretrieval geojson progress tqdm geoglows pygeos noise pillow=9.0.1 rasterio
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
The first is by creating an input JSON file like the one below:

```json
{
    "dams": [
        {
            "name": "272",
            "dam_csv":"C:/Users/jlgut/OneDrive/Desktop/rathcelon/examples/LHD_lat_long.csv",
            "dam_id_field": "BYU_SITE_ID",
            "dam_id": 272,
            "flowline": "C:/Users/jlgut/Downloads/streams_712.gpkg",
            "dem_dir": "C:/Users/jlgut/OneDrive/Desktop/rathcelon_test_data/DEM",
            "bathy_use_banks": false,
            "output_dir": "C:/Users/jlgut/OneDrive/Desktop/rathcelon_test_data/Results",
            "process_stream_network": false,
            "find_banks_based_on_landcover": false,
            "create_reach_average_curve_file": false
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
    "create_reach_average_curve_file": (true/false): If true, RathCelon with direct ARC to generate a curve file that is the same for all stream cells on a stream reach. If false, RathCelon will direct ARC to generate a curve file that has a local curve for each stream cell. 

To run with the JSON, issue the following command in you command line window:

```bash
rathcelon json "path\to\your\input.json"
```
### Running from command line
You can also run RathCelon from the command line without the JSON.

This can be done by issueing the following command:

```bash
rathcelon cli name dam_csv dam_id_field dam_id flowline dem_dir output_dir --bathy_use_banks --process_stream_network --find_banks_based_on_landcover --create_reach_average_curve_file
```
Issuing the commands `--bathy_use_banks`, `--process_stream_network`, `--find_banks_based_on_landcover`, and `--create_reach_average_curve_file` indicates that those options are set to True. 

## Outputs
In the output_dir you specify in your inputs, you will find a folder that is the name you specify in your inputs. In the name folder, will be several subfolders. The "VDT" folder contains the shapefiles {name}_Local_VDT_Database.shp and {name}_Local_CurveFile.shp. These will be the 3 nearest VDT database and curvefile stream cells that are 1, 2, and 3x the topwidth distance of the stream away from the dam. The top width was estimated using the 2yr discharge in the "FLOW\{name}_Reanalysis.csv" file and finding the maximum for the stream reach on which the dam is located. 

For a breakdown of the attributes in the curve file, see [this wiki](https://github.com/MikeFHS/automated-rating-curve/wiki/Running-ARC-and-Looking-at-ARC-Outputs). 



