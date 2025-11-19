# RathCelon
A tool for extracting rating curves below a low-head dam. The current code supports the [ECMWF GEOGLOWS Streamflow Service](https://geoglows.ecmwf.int/) and [National Water Model](https://registry.opendata.aws/nwm-archive/) datasets.

## How to Set-up RathCelon
1a. Using Anaconda run:

```bash
conda create -n rathcelon_py310 python=3.10 pyarrow geopandas pandas netcdf4 dask fiona s3fs xarray zarr beautifulsoup4 dataretrieval geojson progress tqdm pygeos pillow=9.0.1 rasterio
````

1b. Alternatively, you can create an environment using conda and then get the libraries using pip:

```bash
conda create -n rathcelon_py310 python=3.10
conda activate rathcelon_py310
pip install --upgrade pip
pip install pyarrow geopandas pandas netCDF4 dask fiona s3fs xarray zarr beautifulsoup4 dataretrieval geojson progress tqdm pygeos pillow==9.0.1 rasterio
pip install memory_profiler
conda install -c conda-forge gdal
```

2.  Activate you conda environment with the command:

<!-- end list -->

```bash
conda activate rathcelon_py310
```

3.  Install the [Automated Rating Curve (ARC) Tool](https://github.com/MikeFHS/automated-rating-curve) in your conda environment by downloading the source code, navigating to the local instance on your machine, and typing the command:

<!-- end list -->

```bash
pip install .
```

4.  Now download this RathCelon repository and run this command from within RathCelon's local directory:

<!-- end list -->

```bash
pip install .
```

5.  Type the following command into the prompt to ensure everything is working and for additional guidance:

<!-- end list -->

```bash
rathcelon -h
```

## Inputs and Runs

RathCelon can be run from a command prompt in two ways: using a JSON input file or via CLI arguments.

### Running with a JSON

The first method is by creating an input JSON file like the one below. Note that if you are on a PC, the `\` in your path need to be changed to `/` in the JSON file.

```json
{
    "dams": [
        {
            "name": "152",
            "dam_csv": "C:/Users/jlgut/OneDrive/Desktop/rathcelon/examples/LHD_lat_long.csv",
            "dam_id_field": "BYU_SITE_ID",
            "dam_id": 152,
            "flowline": "C:/Users/jlgut/Downloads/streams_712.gpkg",
            "dem_dir": "C:/Users/jlgut/OneDrive/Desktop/rathcelon_test_data/DEM_3",
            "output_dir": "C:/Users/jlgut/OneDrive/Desktop/rathcelon_test_data/Results",
            "bathy_use_banks": false,
            "process_stream_network": true,
            "find_banks_based_on_landcover": true,
            "create_reach_average_curve_file": false,
            "known_baseflow": 5.0,
            "known_channel_forming_discharge": null,
            "upstream_elevation_change_threshold": 0.5,
            "streamflow": null
        }
    ]
}
```

#### Input Parameters:

  * **`name`**: (str) The name of the dam you're analyzing.
  * **`dam_csv`**: (str) The path to the CSV containing the locations of the dams.
  * **`dam_id_field`**: (str) The field in the `dam_csv` that stores the unique identifiers for the dams.
  * **`dam_id`**: (int) The ID for this dam in the `dam_id_field`.
  * **`flowline`**: (str) The path to your flowline dataset (e.g., downloaded from [GEOGLOWS](http://geoglows-v2.s3-website-us-west-2.amazonaws.com/#streams/) or NHDPlus).
  * **`dem_dir`**: (str) The path to the folder containing your digital elevation models (DEMs) for the dam.
  * **`output_dir`**: (str) Path to where the outputs of your RathCelon run will be stored.
  * **`bathy_use_banks`**: (bool) Boolean indicating the option in ARC described [here](https://github.com/MikeFHS/automated-rating-curve/wiki/Estimating-Bathymetry-in-ARC#bank-or-water-surface-elevation).
  * **`process_stream_network`**: (bool) If `true`, RathCelon will take the `flowline` input and determine the local instance it will use for analysis. If `false`, it assumes the network has already been processed.
  * **`find_banks_based_on_landcover`**: (bool) If `true`, stream banks are found using [this logic](https://github.com/MikeFHS/automated-rating-curve/wiki/Estimating-Bathymetry-in-ARC#discovering-the-waterfront). If `false`, ARC assumes the stream is a flat surface in the DEM.
  * **`create_reach_average_curve_file`**: (bool) If `true`, ARC generates a curve file that is the same for all stream cells on a reach. If `false`, creates a local curve for each stream cell.
  * **`known_baseflow`**: (float or null) Optional input in cubic meters per second (cms). Used when `bathy_use_banks` is `false`. Substituted into ARC processing to estimate bathymetry and used to calculate top width for distance calculations.
  * **`known_channel_forming_discharge`**: (float or null) Optional input in cms. Used when `bathy_use_banks` is `true`. Substituted into ARC processing to estimate bathymetry.
  * **`upstream_elevation_change_threshold`**: (float) The change in upstream elevation (in meters) used to identify the cross-section upstream of the dam. Default is `1.0`.
  * **`streamflow`**: (str or null) Controls the source of streamflow data:
      * **`null`**: Use for **Pure GEOGLOWS**. RathCelon will fetch data directly from S3 for each segment.
      * **`"path/to/geoglows.gpkg"`**: Use for **Mixed Model**. Specify the GEOGLOWS `.gpkg` here if using NHDPlus flowlines to map NHD segments to GEOGLOWS IDs.
      * **`"path/to/data.parquet"`**: Use for **National Water Model**. Path to a local parquet file containing streamflow data.

To run with the JSON, issue the following command in your command line window:

```bash
rathcelon json "path/to/your/input.json"
```

### Running from command line

You can also run RathCelon from the command line without the JSON.

```bash
rathcelon cli name dam_csv dam_id_field dam_id flowline dem_dir output_dir --bathy_use_banks --process_stream_network --find_banks_based_on_landcover --create_reach_average_curve_file --known_baseflow 5.0 --upstream_elevation_change_threshold 0.1
```

Flags (`--bathy_use_banks`, etc.) set those options to `True`. Omit them to set them to `False`.

You can also specify the streamflow file if needed:

```bash
... --streamflow "path/to/file.gpkg"
```

## Outputs

In the `output_dir` you specify, you will find a folder named after your dam `name` containing several subfolders:

  * **VDT**: Contains `{name}_Local_VDT_Database.shp` and `{name}_Local_CurveFile.shp`. These include the nearest VDT and curvefile stream cells at 1x, 2x, 3x, and 4x the top width distance downstream, plus one upstream cell located where the elevation change meets your threshold.
  * **XS**: Contains `{name}_Local_XS_Lines.shp` with cross-section locations and information for the stream cells in the curve file.
  * **FLOW**: Contains the reanalysis CSV with discharge data used for the dam.

For a breakdown of the attributes in the curve file, see [this wiki](https://github.com/MikeFHS/automated-rating-curve/wiki/Running-ARC-and-Looking-at-ARC-Outputs).
