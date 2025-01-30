# RathCelon
A tool for extracting cross-sections below a low-head dam

## How to Set-up
1. Using Anaconda run:

`bash
conda create -n rathcelon_py310 python=3.10 pyarrow geopandas pandas netcdf4 dask fiona s3fs xarray zarr beautifulsoup4 dataretrieval geojson progress tqdm geoglows pygeos noise pillow=9.0.1 rasterio
`

2. Activate you conda environment with the command:
`bash
conda activate rathcelon_py310
`

3. Install the [Automated Rating Curve (ARC) Tool](https://github.com/MikeFHS/automated-rating-curve) in your conda environment by downloading the source code, navigating to the local instance on your machine, and typing the command:

`bash
pip install .
`

4. Now download this RathCelon repository and run this command from within RathCelon's local directory:
`bash
pip install .
`

5. Type the followng command into the prompt to ensure everything is working and for additional guidance
`bash
rathcelon -h
`

