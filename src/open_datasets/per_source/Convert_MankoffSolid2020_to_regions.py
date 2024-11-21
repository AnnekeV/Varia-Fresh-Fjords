#!/usr/bin/env python
# coding: utf-8

# # preamble

# In[4]:


import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import geopandas as gpd
from shapely import wkt
import matplotlib.pyplot as plt
from geodatasets import get_path
import numpy as np
import xarray as xr

import sys
sys.path.append('..')
from scripts.dicts import *
from scripts.paths import *
from tqdm import tqdm


# In[5]:


dfRegionD = pd.read_csv(path_Mankoff2020Solid + "region_D.csv", index_col=0)
dfErrorRegionD = pd.read_csv(path_Mankoff2020Solid + "region_err.csv", index_col=0)


# # GATE

# In[6]:


dfgatemeta = pd.read_csv(path_Mankoff2020Solid + "gate_meta.csv", index_col=0)


# In[7]:


dfgadeD = pd.read_csv(path_Mankoff2020Solid + "gate_D.csv", index_col=0)


# In[9]:


folder_base = "/Volumes/imau02/rapid/Anneke/"

fpath_adj_sect  = f'{pathDataTemp}adjusted_section_numbers_slater.nc'
ds_mask_sections = xr.open_dataset(fpath_adj_sect)

fpath_masks1k = pathDataGithub + 'temp/masks1k.nc'
mask1k = xr.open_dataset(fpath_masks1k)



# In[23]:


ds_mask_sections = ds_mask_sections.assign_coords(x=ds_mask_sections.x.astype(int).values)
ds_mask_sections = ds_mask_sections.assign_coords(y=ds_mask_sections.y.astype(int).values)
# assign lat and lon
ds_mask_sections['LON']= ds_mask_sections['section_numbers_adjusted'].copy()
ds_mask_sections['LAT']= ds_mask_sections['section_numbers_adjusted'].copy()
ds_mask_sections['LON'].values = mask1k['LON'].values
ds_mask_sections['LAT'].values = mask1k['LAT'].values


# In[25]:



ds = ds_mask_sections


# Load your dataset (assuming it's already loaded as 'ds')
# ds = xr.open_dataset('your_dataset.nc')

# Define the haversine function
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points
    on the earth specified in decimal degrees.
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of Earth in kilometers
    return c * r


# Target coordinates
target_lon = -50
target_lat = 80

def give_section(lat, lon, ds):
    ''' Function to give the section number closest to the target coordinates
    ds = xarray dataset with the section numbers and coordinates
    '''
    # Apply the haversine function to each grid point in the dataset
    distances =  haversine(
        ds.LON.values.flatten(),
        ds.LAT.values.flatten(),
        lon,
        lat,
        )
    arg = np.argmin(distances)
    # print(ds.LON.values.flatten()[arg], ds.LAT.values.flatten()[arg] )
    # print(f"{ds['distance_to_target'].values.flatten()[arg]:.0f} m")
    return ds['section_numbers_adjusted'].values.flatten()[arg]

give_section(target_lat, target_lon, ds)


# In[57]:


section_number = []
for i in tqdm(range(len(dfgatemeta))):
   nr = give_section(dfgatemeta['lat'].iloc[i], dfgatemeta['lon'].iloc[i], ds)
   section_number.append(nr)
dfgatemeta['section_number'] = section_number


# In[58]:


dfgatemeta.to_csv(path_Mankoff2020Solid_adjusted + "gate_meta_section_nr_adjusted.csv")


# In[ ]:


dfErrorGate = pd.read_csv(path_Mankoff2020Solid + "gate_err.csv", index_col=0)
dfCoverageGate = pd.read_csv(path_Mankoff2020Solid + "gate_coverage.csv", index_col=0)


# In[60]:


dfSectionD = pd.DataFrame(index=dfgadeD.index, columns=dfRegionD.columns)
dfSectionErr = pd.DataFrame(index=dfgadeD.index, columns=dfRegionD.columns)
dfSectionCoverage = pd.DataFrame(index=dfgadeD.index, columns=dfRegionD.columns)    

for section in range(1,7+1):
    sectionName = dict_sections[section]
    gates = dfgatemeta[dfgatemeta['section_number']==section].index
    dfgadeD.columns = dfgadeD.columns.astype(int) 
    dfErrorGate.columns = dfErrorGate.columns.astype(int)
    dfCoverageGate.columns = dfCoverageGate.columns.astype(int)
    dfSectionCoverage[sectionName] = dfCoverageGate.loc[:, gates.values].sum(axis=1)
    sectionName = dict_sections[section]
    dfSectionD[sectionName] = dfgadeD.loc[:, gates.values].sum(axis=1)
    dfSectionErr[sectionName] = dfErrorGate.loc[:, gates.values].sum(axis=1)
    dfSectionCoverage[sectionName] = dfCoverageGate.loc[:, gates.values].sum(axis=1)



    


# In[62]:


dfSectionD.to_csv(path_Mankoff2020Solid_adjusted + "section_D.csv")
dfSectionErr.to_csv(path_Mankoff2020Solid_adjusted + "section_err.csv")
dfSectionCoverage.to_csv(path_Mankoff2020Solid_adjusted + "section_coverage.csv")

