
from open_preprocess_racmo import *
from paths import *
import os

import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import matplotlib.colors as mcolors

import numpy as np


import pandas as pd

time_resolution = 'Monthly'
spatial_resolution = '1k'

years= ['2011', '2022']


import time

print("\nReading RACMO data")
start_time = time.time()
dsRunoff = read_RACMO(spatial_resolution, time_resolution, years=years, variable="runoff").rename({'runoffcorr': 'runoff'})
print(f"{ time.time() - start_time:.0f}s")

print("\nConverting to volume")
dsRunoffIce = volume(mask_data(dsRunoff['runoff'], "PROMICE_Ice_caps", spatial_resolution), spatial_resolution)
end_time = time.time()
print(f"{ time.time() - start_time:.0f}s")


def mask_MougBasins_ice(ds, IceOrTundra):
    '''Mask the RACMO data with the Mouginot basins and sum the values for each basin
    IceOrTundra: "Ice" or "Tundra" to select the mask to use'''
    if spatial_resolution == '5_5k':

        mask55 = open_mask_5_5k(spatial_resolution)
        mask55['rlat'] = ds['rlat']
        mask55['rlon'] = ds['rlon']

        if IceOrTundra == 'Ice':
            ds['Mouginot_basins'] = mask55['Mouginot_basins']
            dsMougTime = ds.groupby('Mouginot_basins').sum().isel(height=0)
        elif IceOrTundra == 'Tundra':
            ds['Mouginot_basins'] = mask55['Mouginot_Tundra']
            dsMougTime = ds.groupby('Mouginot_basins').sum().isel(height=0)
        ds.attrs['Description'] = f"Sum per sector of the Mouginot basins for {spatial_resolution} RACMO2.3p2, for the {IceOrTundra} "
        return dsMougTime
        
    elif spatial_resolution == '1k':
        if not 'masks1k' in globals():
            mask1k = open_mask_1k()
        path_mask_1k_with_tundra = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/GrIS_topo_icemask_lsm_tundra_basins_lon_lat_1km.nc"
        mask_1k_with_tundra = xr.open_dataset(path_mask_1k_with_tundra)
        if IceOrTundra == 'Ice':
            ds['Mouginot_basins'] = mask_1k_with_tundra['Basins_All_Greenland']
            dsMougTime = ds.groupby('Mouginot_basins').sum()
        elif IceOrTundra == 'Tundra':
            ds['Mouginot_basins'] = mask_1k_with_tundra['Tundra_basins']
            dsMougTime = ds.groupby('Mouginot_basins').sum()
        # give a attribute
        ds.attrs['Description'] = f"Sum per sector of the Mouginot basins for {spatial_resolution} RACMO2.3p2, for the {IceOrTundra} "
        return dsMougTime




dict_Moug = {0:"Not GRIS", 1: "NO", 2: "NE", 3: "CE", 4: "SE", 5: "SW", 6: "CW", 7: "NW"}

dsRunoffIceSector = mask_MougBasins_ice(dsRunoffIce, "Ice")


# now write to pathdatatemp

filename = f"runoff_monthly_FGRN055_Downscaled_RACMO2.3p2_RACMO2.3p2_{spatial_resolution}_runoff_{time_resolution}_sector_{years[0]}_{years[-1]}_Ice_Caps.nc"
path_file = os.path.join(pathAnnekeFolderIMAU02, "Downscaling_GR/Monthly", filename)
print("Writing to ", path_file)


dsRunoffIceSector.to_netcdf(path_file)

print("\nDone, in total ", f"{ time.time() - start_time:.0f}s")