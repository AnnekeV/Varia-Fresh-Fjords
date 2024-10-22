#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import numpy as np
import glob as glob
import xarray as xr
import numpy as np
from scipy import interpolate
import pandas as pd
import numpy as np
import xarray as xr

import sys

sys.path.append("..")
from scripts.dicts import *
from scripts.paths import *


# In[2]:


path_BasalMelt = "/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/FWclean/data/raw/liquid/Karlsson2021_basalmelt/basalmelt.nc"


# In[3]:


fpath_adj_sect = pathDataGithub + "temp/adjusted_section_numbers_slater.nc"
fpath_masks1k = pathDataGithub + "temp/masks1k.nc"
folder_MARRACMO1km = pathDataGithub + "raw/liquid/"

ds_adj_sect = xr.open_dataset(fpath_adj_sect)
ds_masks1k = xr.open_dataset(fpath_masks1k)


# In[4]:


ds_Basal = xr.open_dataset(path_BasalMelt)


# In[5]:


lats_target = ds_Basal["latitude"].values
lons_target = ds_Basal["longitude"].values

lats_source = ds_masks1k.LAT.values
lons_source = ds_masks1k.LON.values

values_source = ds_adj_sect["section_numbers_adjusted"].values

lats_target.shape, lons_target.shape, lats_source.shape, lons_source.shape, values_source.shape

adj_sect_basal = interpolate.griddata(
    (lats_source.flatten(), lons_source.flatten()),
    values_source.flatten(),
    (lats_target, lons_target),
    method="nearest",
)


# In[6]:

ds_Basal["adj_sect"] = xr.DataArray(
    adj_sect_basal, dims=["y", "x"], coords={"y": ds_Basal.y, "x": ds_Basal.x}
)
ds_Basal["adj_sect"].plot()


# In[7]:
df_total_melt_per_basin = (
    ds_Basal["totalmelt"].groupby(ds_Basal["adj_sect"]).sum().to_dataframe() / 1e3
)
df_total_melt_per_basin.plot(kind="bar")

df_total_melt_per_basin.round(1)

df_melt_per_basin = df_total_melt_per_basin.copy()


# Also for individual components

for melt_type in ["gfmelt", "fricmelt", "vhdmelt"]:
    df_melt_per_basin[melt_type] = (
        ds_Basal[melt_type].groupby(ds_Basal["adj_sect"]).sum().to_dataframe() / 1e3
    )
df_melt_per_basin


# # Other way around



lats_source = ds_Basal["latitude"].values
lons_source = ds_Basal["longitude"].values
values_source = ds_Basal["totalmelt"].values

lats_target = ds_masks1k.LAT.values
lons_target = ds_masks1k.LON.values

melt_basal = interpolate.griddata(
    (lats_source.flatten(), lons_source.flatten()),
    values_source.flatten(),
    (lats_target, lons_target),
    method="nearest",
)
melt_basal = np.reshape(melt_basal, lats_target.shape)
melt_basal = xr.DataArray(
    melt_basal, dims=["y", "x"], coords={"y": ds_adj_sect.y, "x": ds_adj_sect.x}
)




# In[12]:


df_melt_basal_per_basin_racmo_grid = (
    melt_basal.groupby(ds_adj_sect["section_numbers_adjusted"])
    .sum()
    .to_dataframe(name="Total melt (m/yr)")
    / 1e3
)
print(f"Sum is {df_melt_basal_per_basin_racmo_grid.sum().values[0]:.1f} km3/yr")
df_melt_basal_per_basin_racmo_grid.round(1)




# In[14]:
df_monthly_fluxes = pd.read_csv(
    "/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/FWclean/data/processed/Seasonal_Greenland_2009_2022.csv",
    index_col=0,
)
df_monthly_fraction_runoff_of_total_runoff = (
    df_monthly_fluxes / df_monthly_fluxes.sum(axis=0)
)["Liquid Runoff Ice Sheet"]


# In[15]:
folder = "/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/FWclean/data/processed/Seasonal cycle per sector/"
files = glob.glob(folder + "*.csv")

df_basin_monthly_runoff = pd.DataFrame()
for f in files:
    df = pd.read_csv(f, index_col=0)["Liquid Runoff Ice Sheet"].rename(
        f.split("/")[-1].split("_")[0]
    )
    df_basin_monthly_runoff = pd.concat([df_basin_monthly_runoff, df], axis=1)
df_basin_monthly_runoff_fraction = (
    df_basin_monthly_runoff / df_basin_monthly_runoff.sum(axis=0)
)
df_basin_monthly_runoff_fraction



df_vhd_monthly_basin = (
    df_melt_per_basin.rename(index=dict_sections)["vhdmelt"]
    * df_basin_monthly_runoff_fraction
)
monthly_equal = pd.DataFrame(
    index=df_basin_monthly_runoff_fraction.index,
    columns=df_basin_monthly_runoff_fraction.columns,
    data=1 / 12,
)
df_gf_monthly_basin = (
    df_melt_per_basin.rename(index=dict_sections)["gfmelt"] * monthly_equal
)
df_fric_monthly_basin = (
    df_melt_per_basin.rename(index=dict_sections)["fricmelt"] * monthly_equal
)
df_Basal_basin_monthly = df_vhd_monthly_basin + df_gf_monthly_basin + df_fric_monthly_basin
df_Basal_GrIS_monthly = df_Basal_basin_monthly.sum(axis=1)
