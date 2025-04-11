# Package Imports
import os # Interoperable file paths
from glob import glob
import pathlib 
import requests

# Data Tools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Spatial tools
import geopandas as gpd
import rioxarray as rxr
import rioxarray.merge as rxrm
import xarray as xr
import xrspatial

# GBIF packages
import pygbif.occurrences as occ
import pygbif.species as species
from getpass import getpass
import time
import zipfile

# Visualizing Tools
import holoviews as hv
import hvplot.pandas
import hvplot.xarray



def convert_longitude(longitude):
    """Convert longitude range from 0-360 to [-180:180]"""
    return (longitude - 360) if longitude > 180 else longitude



def derive_slope(site_gdf, elevation_da):
    # Estimate the UTM CRS based on the bounds of the site
    epsg_utm = site_gdf.estimate_utm_crs()
    # Reproject boundary and elevation data so that units are in meters
    elevation_proj_da = elevation_da.rio.reproject(epsg_utm)
    site_proj_gdf = site_gdf.to_crs(epsg_utm)
    
    # Calculate slope using xrspatial
    slope_full_da = xrspatial.slope(elevation_proj_da)
    slope_da = slope_full_da.rio.clip(site_proj_gdf.geometry)

    return slope_da



# Download .tif files
def download_tif(url, save_path):
    """Function to download and save .tif files"""
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(save_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        print(f"Downloaded: {save_path}")
    else: 
        print(f"Failed to download {url}. Status code: {response.status_code}")



# Download climate data
def MACAv2_data(variable, scenario, start_year, site_gdf):
    end_year = 2099 if start_year == 2096 else (start_year + 4)
    maca_url = (
        'http://thredds.northwestknowledge.net:8080/thredds/dodsC/MACAV2/BNU-ESM'
        f'/macav2metdata_{variable}_BNU-ESM_r1i1p1_'
        f'{scenario}_{start_year}_{end_year}_CONUS_monthly.nc')
    maca_da = xr.open_dataset(maca_url).squeeze().precipitation
    bounds = site_gdf.to_crs(maca_da.rio.crs).total_bounds
    # Reassign coordinates to [-180:180]
    maca_da = maca_da.assign_coords(
        lon=("lon", [convert_longitude(l) for l in maca_da.lon.values])
    )
    maca_da = maca_da.rio.set_spatial_dims(x_dim='lon', y_dim='lat')
    maca_da = maca_da.rio.clip_box(*bounds)


    return maca_da



# Plot boundary line of a gdf
def plot_site(gdf):
    gdf_plot = (
        gdf.dissolve().hvplot(
            geo=True, tiles='EsriImagery',
            title=gdf.FORESTNAME.values[0],
            fill_color=None, line_color='white', line_width=3,
            frame_width=600
        )
    )
    return gdf_plot



# Download soil data
def polaris_soil_data(variable, statistic, depth, gdf, save_folder):
    
    # Ensure the save directory exists
    polaris_dir = os.path.join(data_dir, 'polaris/'+save_folder)
    os.makedirs(polaris_dir, exist_ok=True)

    # Create list of soil data URLs
    soil_url_template = ("http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0"
                "/{variable}"
                "/{statistic}"
                "/{depth}"
                "/lat{min_lat}{max_lat}_lon{min_lon}{max_lon}.tif")

    bounds_min_lon, bounds_min_lat, bounds_max_lon, bounds_max_lat = (
        gdf.total_bounds)

    soil_url_list = []

    for lesser_lat in range(floor(bounds_min_lat), ceil(bounds_max_lat)):
        for lesser_lon in range(floor(bounds_min_lon), ceil(bounds_max_lon)):
            soil_url = soil_url_template.format(
                variable=variable,
                statistic=statistic,
                depth=depth,
                min_lat=lesser_lat, max_lat=lesser_lat+1, 
                min_lon=lesser_lon, max_lon=lesser_lon+1)
            soil_url_list.append(soil_url)

    # Download the files locally
    for url in soil_url_list:
        # Extract file name from url
        split_url = url.split("/")[-5:]
        filename = ''
        for i in split_url:
            filename+=i

        # Download the .tif file once
        save_path = os.path.join(polaris_dir, filename)
        if not os.path.exists(save_path):
            download_tif(url, save_path)



def process_image(file_pattern, boundary_gdf, buffer=0):
    """Load image, crop to study boundary, merge assays"""
    
    # Set a buffer to crop images to
    bounds_buffer = set_buffer(boundary_gdf, buffer)
    
    # Open and crop the images
    da_list = []
    for file_path in glob(file_pattern):
        tile_da = (
            rxr.open_rasterio(file_path, mask_and_scale=True)
            .squeeze())
        cropped_da = tile_da.rio.clip_box(*bounds_buffer)
        da_list.append(cropped_da)
    
    # Merge the list of cropped data arrays
    merged_da = rxrmerge.merge_arrays(da_list)
    
    # Returns a cropped and merged data array
    return merged_da



def set_buffer(boundary_gdf, buffer=0):
    """
    Increases the max bounds of a geo data frame by a set amount.
    Returns the max bounds as a tuple.
    """
    bounds = tuple(boundary_gdf.total_bounds)
    xmin, ymin, xmax, ymax = bounds
    bounds_buffer = (xmin-buffer, ymin-buffer, xmax+buffer, ymax+buffer)

    return bounds_buffer



# Download elevation data
def srtm_data_download(boundary_gdf, buffer=0.025):
    elevation_dir = os.path.join(data_dir, 
                                  'srtm',
                                  boundary_gdf.FORESTNAME.values[0])
    os.makedirs(elevation_dir, exist_ok=True)

    srtm_pattern = os.path.join(elevation_dir, '*.hgt.zip')

    bounds_buffer = set_buffer(boundary_gdf, buffer)

    if not glob(srtm_pattern):
        srtm_results = earthaccess.search_data(
            short_name='SRTMGL1',
            bounding_box=bounds_buffer
        )
        srtm_results = earthaccess.download(srtm_results, elevation_dir)
    
    srtm_da_list = []

    for srtm_path in glob(srtm_pattern):
        tile_da = rxr.open_rasterio(srtm_path, mask_and_scale=True).squeeze()
        cropped_da = tile_da.rio.clip_box(*bounds_buffer)
        srtm_da_list.append(cropped_da)
    
    srtm_da =rxrmerge.merge_arrays(srtm_da_list)

    return srtm_da