from setuptools import setup, find_packages, Extension

setup(
    name='rathcelon',  # Replace with a unique name for your package
    version='0.0.0',
    description='A library for creating rating curves downstream of low-head dams',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Joseph Gutenson',
    author_email='joseph@follumhydro.com',
    url='https://github.com/jlgutenson/rathcelon',  # Update with your repo URL
    packages=find_packages(),
    install_requires=[
        'beautifulsoup4',
        'dask',
        'dataretrieval',
        'fiona',
        'gdal',
        'geoglows',
        'geojson',
        'geopandas',
        'netCDF4',
        'noise',
        'numpy',
        'pandas',
        'pillow==9.0.1',
        'progress',
        'pygeos',
        'pyproj',
        'rasterio',
        'requests',
        's3fs',
        'scipy',
        'shapely',
        'tqdm',
        'xarray'
    ],
    entry_points={
        "console_scripts": [
        "rathcelon=rathcelon.main:main",
        ],
    },
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD-3 License',
        'Operating System :: OS Independent',
    ],

)