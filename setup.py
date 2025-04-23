from setuptools import setup, find_packages, Extension

setup(
    name='rathcelon',  
    version='0.1.0',
    description='A library for creating rating curves downstream of low-head dams',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Joseph Gutenson',
    author_email='joseph@follumhydro.com',
    url='https://github.com/jlgutenson/rathcelon',  
    packages=find_packages(),
    install_requires=[
        'beautifulsoup4',
        'dask',
        'dataretrieval',
        'fiona',
        'gdal',
        'geojson',
        'geopandas',
        'netCDF4',
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