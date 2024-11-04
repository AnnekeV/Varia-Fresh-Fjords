#!/usr/bin/env python
# coding: utf-8
# In[501]:

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import geopandas as gpd
from shapely import wkt
import matplotlib.pyplot as plt
from geodatasets import get_path
import numpy as np
import xarray as xr
from scipy import interpolate

import sys
sys.path.append("..")
from scripts.dicts import *
from scripts.great_circle_distance import great_circle_distance
from scripts.paths import *

# Load datasets using paths from the second script
ds = xr.open_dataset(f"{pathIMAU01}/RACMO2.3p2/FGRN055_1940/Monthly/runoff_monthlyS_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202312.nc")
ds_runoff = ds.sel(time='2023')['runoff']
ds_lat = ds['lat']
ds_lon = ds['lon']

# File paths for adjusted section numbers and masks
fpath_adj_sect = f"{pathDataTemp}adjusted_section_numbers_slater.nc"
fpath_masks1k = f"{pathDataTemp}masks1k.nc"

# Load section numbers and masks
ds_adj_sect = xr.open_dataset(fpath_adj_sect)
ds_masks1k = xr.open_dataset(fpath_masks1k)

# Interpolate runoff data
lats_source = ds.lat.values
lons_source = ds.lon.values
lats_target = ds_masks1k.LAT.values
lons_target = ds_masks1k.LON.values
ds_runoff_tundra2023 = xr.Dataset()

for month in range(12):
    values_source = ds_runoff.isel(time=month).values
    time = ds_runoff.time.isel(time=month).values
    runoff_tundra = interpolate.griddata(
        (lats_source.flatten(), lons_source.flatten()), 
        values_source.flatten(),
        (lats_target, lons_target), 
        method='nearest'
    )
    runoff_tundra = np.reshape(runoff_tundra, lats_target.shape)
    da_runoff_tundra = xr.DataArray(runoff_tundra, dims=['y', 'x'], coords={'y': ds_adj_sect.y, 'x': ds_adj_sect.x}, name='Liquid Runoff Tundra').expand_dims(time=[time])
    ds_runoff_tundra2023 = xr.combine_by_coords([ds_runoff_tundra2023, da_runoff_tundra])

# Save output using the defined paths
ds_runoff_tundra2023.to_netcdf(f"{pathIMAU02}RACMO2.3p2/FGRN055/Downscaling_GR/Monthly/runoff_tundra.2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.tundra.nc")
ds_runoff_tundra2023.to_netcdf(f"{pathDataTemp}runoff_tundra.2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.tundra.nc")

# %%
