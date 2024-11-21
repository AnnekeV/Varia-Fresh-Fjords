#!/usr/bin/env python
# coding: utf-8

### Introduction

# This script processes the MankoffSolid2020 dataset to assign section numbers used in this study to gates based on their geographical coordinates. It calculates the total discharge, error, and coverage for each section by summing the values from the gates within each section. The results are saved to CSV files for further analysis.
# Anneke Vries Nov 2024



import pandas as pd
import numpy as np
import xarray as xr
from tqdm import tqdm
import sys

sys.path.append('..')
from scripts.dicts import *
from scripts.paths import *

# Load datasets
dfRegionD = pd.read_csv(path_Mankoff2020Solid + "region_D.csv", index_col=0)
dfErrorRegionD = pd.read_csv(path_Mankoff2020Solid + "region_err.csv", index_col=0)
dfgatemeta = pd.read_csv(path_Mankoff2020Solid + "gate_meta.csv", index_col=0)
dfgadeD = pd.read_csv(path_Mankoff2020Solid + "gate_D.csv", index_col=0)
dfErrorGate = pd.read_csv(path_Mankoff2020Solid + "gate_err.csv", index_col=0)
dfCoverageGate = pd.read_csv(path_Mankoff2020Solid + "gate_coverage.csv", index_col=0)

fpath_adj_sect = f'{pathDataTemp}adjusted_section_numbers_slater.nc'
ds_mask_sections = xr.open_dataset(fpath_adj_sect)

fpath_masks1k = pathDataGithub + 'temp/masks1k.nc'
mask1k = xr.open_dataset(fpath_masks1k)

# Assign coordinates
ds_mask_sections = ds_mask_sections.assign_coords(x=ds_mask_sections.x.astype(int).values)
ds_mask_sections = ds_mask_sections.assign_coords(y=ds_mask_sections.y.astype(int).values)
ds_mask_sections['LON'] = mask1k['LON']
ds_mask_sections['LAT'] = mask1k['LAT']

# Define the haversine function
def haversine(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of Earth in kilometers
    return c * r

# Function to find the closest section number to the target coordinates
def give_section(lat, lon, ds):
    distances = haversine(ds.LON.values.flatten(), ds.LAT.values.flatten(), lon, lat)
    arg = np.argmin(distances)
    return ds['section_numbers_adjusted'].values.flatten()[arg]

# Assign section numbers to gates
section_number = []
for i in tqdm(range(len(dfgatemeta))):
    nr = give_section(dfgatemeta['lat'].iloc[i], dfgatemeta['lon'].iloc[i], ds_mask_sections)
    section_number.append(nr)
dfgatemeta['section_number'] = section_number

# Save the updated gate metadata
dfgatemeta.to_csv(path_Mankoff2020Solid_adjusted + "gate_meta_section_nr_adjusted.csv")

# Calculate total discharge, error, and coverage for each section
dfSectionD = pd.DataFrame(index=dfgadeD.index, columns=dfRegionD.columns)
dfSectionErr = pd.DataFrame(index=dfgadeD.index, columns=dfRegionD.columns)
dfSectionCoverage = pd.DataFrame(index=dfgadeD.index, columns=dfRegionD.columns)

for section in range(1, 8):
    sectionName = dict_sections[section]
    gates = dfgatemeta[dfgatemeta['section_number'] == section].index
    dfgadeD.columns = dfgadeD.columns.astype(int)
    dfErrorGate.columns = dfErrorGate.columns.astype(int)
    dfCoverageGate.columns = dfCoverageGate.columns.astype(int)
    dfSectionCoverage[sectionName] = dfCoverageGate.loc[:, gates.values].sum(axis=1)
    dfSectionD[sectionName] = dfgadeD.loc[:, gates.values].sum(axis=1)
    dfSectionErr[sectionName] = dfErrorGate.loc[:, gates.values].sum(axis=1)

# Save the results to CSV files
dfSectionD.to_csv(path_Mankoff2020Solid_adjusted + "section_D.csv")
dfSectionErr.to_csv(path_Mankoff2020Solid_adjusted + "section_err.csv")
dfSectionCoverage.to_csv(path_Mankoff2020Solid_adjusted + "section_coverage.csv")
```
