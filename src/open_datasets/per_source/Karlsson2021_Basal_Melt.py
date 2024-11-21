#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import glob
import xarray as xr
from scipy import interpolate
import sys

sys.path.append("..")
from scripts.dicts import *
from scripts.paths import *






# Define paths
# Basal melt from : Karlsson, N.B., Solgaard, A.M., Mankoff, K.D. et al. A first constraint on basal melt-water production of the Greenland ice sheet. Nat Commun 12, 3461 (2021). https://doi-org.utrechtuniversity.idm.oclc.org/10.1038/s41467-021-23739-z
# mask1k from racmo mask downscaled to 1km for Promice mask
 
path_BasalMelt = pathDataRaw + "liquid/Karlsson2021_basalmelt/basalmelt.nc"
fpath_adj_sect = pathDataTemp + "adjusted_section_numbers_slater.nc"
fpath_masks1k = pathDataTemp + "masks1k.nc"

# Load datasets
ds_adj_sect = xr.open_dataset(fpath_adj_sect)
ds_masks1k = xr.open_dataset(fpath_masks1k)
ds_Basal = xr.open_dataset(path_BasalMelt)

# Interpolate section numbers to Basal melt grid
lats_target = ds_Basal["latitude"].values
lons_target = ds_Basal["longitude"].values
lats_source = ds_masks1k.LAT.values
lons_source = ds_masks1k.LON.values
values_source = ds_adj_sect["section_numbers_adjusted"].values

adj_sect_basal = interpolate.griddata(
    (lats_source.flatten(), lons_source.flatten()),
    values_source.flatten(),
    (lats_target, lons_target),
    method="nearest",
)

ds_Basal["adj_sect"] = xr.DataArray(
    adj_sect_basal, dims=["y", "x"], coords={"y": ds_Basal.y, "x": ds_Basal.x}
)

# Calculate total melt per basin
df_total_melt_per_basin = (
    ds_Basal["totalmelt"].groupby(ds_Basal["adj_sect"]).sum().to_dataframe() / 1e3
)
df_melt_per_basin = df_total_melt_per_basin.copy()

# Calculate individual components
for melt_type in ["gfmelt", "fricmelt", "vhdmelt"]:
    df_melt_per_basin[melt_type] = (
        ds_Basal[melt_type].groupby(ds_Basal["adj_sect"]).sum().to_dataframe() / 1e3
    )

# Interpolate total melt to RACMO grid
values_source = ds_Basal["totalmelt"].values
melt_basal = interpolate.griddata(
    (lats_target.flatten(), lons_target.flatten()),
    values_source.flatten(),
    (lats_source, lons_source),
    method="nearest",
)
melt_basal = np.reshape(melt_basal, lats_source.shape)
melt_basal = xr.DataArray(
    melt_basal, dims=["y", "x"], coords={"y": ds_adj_sect.y, "x": ds_adj_sect.x}
)

df_melt_basal_per_basin_racmo_grid = (
    melt_basal.groupby(ds_adj_sect["section_numbers_adjusted"])
    .sum()
    .to_dataframe(name="Total melt (m/yr)")
    / 1e3
)
print(f"Sum is {df_melt_basal_per_basin_racmo_grid.sum().values[0]:.1f} km3/yr")


# Calculate monthly basal melt per basin, based on monthly runoff
# Load monthly fluxes
df_monthly_fluxes = pd.read_csv(
    pathDataProcessed + "Seasonal_Greenland_2009_2022.csv", index_col=0
)
df_monthly_fraction_runoff_of_total_runoff = (
    df_monthly_fluxes / df_monthly_fluxes.sum(axis=0)
)["Liquid Runoff Ice Sheet"]

# Load basin monthly runoff
folder = pathDataProcessed + "Archive/Seasonal cycle per sector without basal/"
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

# Calculate monthly basal melt per basin
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
df_Basal_GrIS_annual = df_Basal_GrIS_monthly.sum()

