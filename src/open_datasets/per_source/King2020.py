#!/usr/bin/env python
# coding: utf-8

# In[521]:


import pandas as pd
import geopandas as gpd
from shapely import wkt
import matplotlib.pyplot as plt
from geodatasets import get_path
import numpy as np
import xarray as xr

import sys
sys.path.append('..')
from scripts.dicts import *
from scripts.great_circle_distance import great_circle_distance
from scripts.paths import *



# In[522]:


pathKing =  pathDataGithub +'raw/solid_discharge/King2020_doi_10_5061_dryad_qrfj6q5cb__v20200429/GrIS_D_1985-2018.xlsx'


# In[523]:


fpath_adj_sect  = pathDataGithub +'temp/adjusted_section_numbers_slater.nc'
fpath_masks1k = pathDataGithub + 'temp/masks1k.nc'
folder_MARRACMO1km = pathDataGithub + "raw/liquid/"

ds_adj_sect = xr.open_dataset(fpath_adj_sect)
ds_masks1k = xr.open_dataset(fpath_masks1k)


# In[524]:


def find_corresponding_section(ds_section, da_lat, da_lon, target_lat, target_lon):
    """Find the section in ds_section that is closest to the target location (target_lat, target_lon).
    
    Parameters
    ----------
    ds_section : xarray.Dataset
        Dataset with sections as values.
    da_lat : xarray.DataArray
        Latitude values of the sections.
    da_lon : xarray.DataArray
        Longitude values of the sections.
    target_lat : float
        Latitude of the target location.
    target_lon : float
        Longitude of the target location.
    
    Returns
    -------
    xarray.DataArray
        The section that is closest to the target location.
    """
    # Calculate the distance between the target location and all sections
    distances = great_circle_distance(da_lat, da_lon, target_lat, target_lon)
    indx_to_min = distances.argmin(dim=['x', 'y'])
    
    # Find the section with the smallest distance
    closest_section = int(ds_section[indx_to_min].values)

    # assert that closest distance is less than 20 km
    assert distances.min() < 20, f"Closest distance is {distances.min():.2f} km, which is more than 20 km"
    
    return closest_section


find_corresponding_section(ds_adj_sect['section_numbers_adjusted'], ds_masks1k.LAT, ds_masks1k.LON, 64.1, -50.5)


# # Import metadata s.a. longitude etcetera and match with basin_nr

# In[525]:


ds_King_meta = pd.read_excel(pathKing, header=[0, 1], nrows=1, index_col=0)
ds_King_meta =ds_King_meta.T.unstack(level=1)
ds_King_meta.columns = ds_King_meta.columns.droplevel()
def remove_quotes(x):
    return x.strip("'")
def remove_quotes_columns(ds):
    ds.columns = remove_quotes(ds.columns.str)
def remove_quotes_index(ds):
    ds.index = remove_quotes(ds.index.str)
remove_quotes_columns(ds_King_meta)
remove_quotes_index(ds_King_meta)
ds_King_meta['basin_nr'] = ds_King_meta.apply(lambda x: find_corresponding_section(ds_adj_sect['section_numbers_adjusted'], ds_masks1k.LAT, ds_masks1k.LON, x.latitude, x.longitude), axis=1)
ds_King_meta


# # Import Discharge Data

# In[526]:


import os
ds_D_king_xlsx = pd.read_excel(pathKing, header=[0,3], index_col=0, parse_dates=True)
glacier_names = remove_quotes(ds_D_king_xlsx.columns.levels[0].str)
ds_D_king_xlsx.columns = ds_D_king_xlsx.columns.set_levels(glacier_names, level=0)
ds_D_king_xlsx.columns = ds_D_king_xlsx.columns.set_levels(remove_quotes(ds_D_king_xlsx.columns.levels[1].str), level=1)
ds_D_king = ds_D_king_xlsx.unstack(level=0).reset_index().rename(columns={'level_0':'basin_nr', 'level_1':'Value_Uncertainty', 0:'Value'})
ds_D_king = pd.merge(ds_D_king, ds_King_meta['basin_nr'].reset_index(), on='Number_Name')
if not os.path.isfile(pathDataTemp + 'King2020_basins.csv'):
    ds_D_king.to_csv(pathDataTemp + 'King2020_basins.csv', index=False)


# # Correction on total values

# In[527]:


# Calculate the fraction of the total discharge that each glacier contributes between 2000 and 2005
df_D_2000_2005_mean = ds_D_king[ds_D_king.Value_Uncertainty=='Discharge (Gt/yr)'][(ds_D_king['Time'] > '2000')&(ds_D_king['Time'] < '2006')].groupby(['Number_Name'])['Value'].mean()
df_D_fraction_per_glacier = df_D_2000_2005_mean /df_D_2000_2005_mean.sum()

# Calculate the coverage of the King dataset
ds_D_king_coverage = ds_D_king[ds_D_king.Value_Uncertainty=='Discharge (Gt/yr)'].pivot(index='Time', columns='Number_Name', values='Value')
ds_D_king_coverage= (ds_D_king_coverage.notna()*df_D_fraction_per_glacier).sum(axis=1)
ds_D_king_coverage[ds_D_king_coverage.index > '2000'] = 1



# In[528]:


df_D_King_monthly_GrIS =  ds_D_king[ds_D_king.Value_Uncertainty=='Discharge (Gt/yr)'].pivot(index='Time', columns='Number_Name', values='Value').sum(axis=1)/(ds_D_king_coverage)
df_D_King_annual_GrIS = df_D_King_monthly_GrIS.resample('YS').mean()


# # Per basin

# In[501]:


ds_D_king_per_basin = ds_D_king.set_index(['Number_Name', 'Value_Uncertainty', 'Time'])
# only select discharge from multiindex
ds_D_king_per_basin = ds_D_king_per_basin.loc[(slice(None), "Discharge (Gt/yr)"), :]
ds_D_king_per_basin_sum_uncorrected = ds_D_king_per_basin.groupby([ 'Time', 'basin_nr']).sum() 
ds_D_king_per_basin_mean = ds_D_king_per_basin.groupby([ 'Time', 'basin_nr']).mean()


# In[453]:





# In[529]:


ds_King_meta['2000_2005_mean'] = df_D_2000_2005_mean
sum_per_basin = ds_King_meta.groupby('basin_nr')['2000_2005_mean'].sum().rename('Basin_sum_2000_2005')
df_D_fraction_of_basin = pd.merge(ds_King_meta.reset_index(), sum_per_basin, on='basin_nr')
df_D_fraction_of_basin['fraction_of_basin'] = df_D_fraction_of_basin['2000_2005_mean']/df_D_fraction_of_basin['Basin_sum_2000_2005']
df_D_fraction_of_basin = df_D_fraction_of_basin.set_index('Number_Name')['fraction_of_basin']

ds_D_king_basin_coverage = (ds_D_king[ds_D_king.Value_Uncertainty=='Discharge (Gt/yr)'].pivot(index='Time', columns='Number_Name', values='Value').notna()*df_D_fraction_of_basin).T.groupby(ds_King_meta['basin_nr'] ).sum().T
ds_D_king_basin_coverage[ds_D_king_basin_coverage.index > '2000'] = 1

ds_D_king_per_basin_sum = (ds_D_king_per_basin_sum_uncorrected.reset_index().pivot(index='Time', columns='basin_nr', values='Value').T/ds_D_king_coverage).T




# In[544]:


ds_D_king_basin_coverage[ds_D_king_basin_coverage>0]['1985': '1995'].mean()





# In[264]:




# In[557]:


ds_D_king_per_basin_sum_annual =  ds_D_king_per_basin_sum.resample('YS').mean()
ds_D_king_sum_annual_total = ds_D_king_per_basin_sum_annual.sum(axis=1)
ds_D_king_sum_annual_total[ds_D_king_sum_annual_total == 0] = np.nan




# In[261]:


nr_per_basin_King = ds_King_meta.groupby('basin_nr').count()
nr_per_basin_King.drop(columns=['latitude']).rename(columns={'longitude':'nr_glaciers'})


