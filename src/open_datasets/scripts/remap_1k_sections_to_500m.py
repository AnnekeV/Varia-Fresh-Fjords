#%%
import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata


folder_base = "/Volumes/imau02/rapid/Anneke/"
# folder_base = '/science/projects/imau02/rapid/Anneke/
folder_CARRA = folder_base + "CARRA/"
folder_RACMO = folder_base + "RACMO2.3p2/"
folder_out = folder_CARRA + "Monthly/on_RACMO_grid"

#%%
print("Opening datasets...")
dsOldGrid = xr.open_dataset( folder_RACMO+ "FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_with_lon_lat_13_June.nc")
# dsNewGrid = xr.open_dataset(folder_RACMO+ "FGRN055/Downscaling_GR_500m/GrIS_topo_icemask_lsm_lon_lat_0.5km.nc")
dsNewGrid = xr.open_mfdataset(folder_RACMO + "FGRN055/Downscaling_GR_500m/GrIS_topo_icemask_lsm_lon_lat_0.5km.nc")


#%% Define data that will be analyzed
# List of years
variablesOld = ['adjusted_sections']
# Map of variable name for every dataset
variable_mapping = {
     'adjusted_sections': {'c': 'section_numbers_adjusted', 'r': 'section_numbers_adjusted'},
}

variable_location = {
    'adjusted_sections' :  {'in':   folder_RACMO+ "FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_with_lon_lat_13_June.nc" , 
    # 'out': folder_RACMO + "Downscaling_GR_500m/Mask_adjusted_section_numbers_slater_0.5km_13_June" }
    'out': "Mask_adjusted_section_numbers_slater_0.5km_13_June.nc" }
    }




#%% Open CARRA data and store in CARRA dictionary with year, month, variable as key


# Initialize an empty dictionary to store all CARRA data
carra_data_dict = {}

for variable in variablesOld:
    # Open CARRA data for the current variable
    carra_data_dict[f'{variable}'] = xr.open_dataset(variable_location[f'{variable}']['in'])
    print('Opened Old ------>', variable)




#%% Plot difference between CARRA and interpolated RACMO
   
for variable in variablesOld:
    ofile = variable_location[variable]['out']
    print('----->' + ofile)
    # check if data is not already available
    if os.path.isfile(ofile): 
        print('file already exists, moving on...')
        continue

    ### Get data ###
    # Extract the carra and racmo data
    dataNew = dsNewGrid
    dataOld = dsOldGrid

    daNew = dsNewGrid["LSM"]
    daOldInterp = xr.DataArray(
                        coords={
                            'y': daNew['y'],
                            'x': daNew['x']
                        },
                        dims=['time', 'y', 'x']
                    )

    print("Preparing grid points...")
    # Prepare grid points        
    points_carra = (dsOldGrid['LON'].values.ravel(), dsOldGrid['LAT'].values.ravel())

    # for t in range(3)
    values_carra = getattr(carra_data_dict[f'{variable}'],variable_mapping[f'{variable}']['c']).values.ravel()

    print("Interpolating...")
    # Interpolate CARRA data to downscaled RACMO grid
    interp_carra_values = griddata(points_carra, values_carra, (dataNew.lon, dataNew.lat), method='nearest') #method='linear')
    daOldInterp.values[:, :] = interp_carra_values[ :, :]

    # set all attributes of the interpolated data array
    daOldInterp.attrs = getattr(carra_data_dict[variable],variable_mapping[variable]['c']).attrs
    ds_carra_interp = daOldInterp.to_dataset(name=variable_mapping[variable]['r'])
    ds_carra_interp.attrs = carra_data_dict[variable].attrs
    ds_carra_interp.to_netcdf(ofile)





# %%
