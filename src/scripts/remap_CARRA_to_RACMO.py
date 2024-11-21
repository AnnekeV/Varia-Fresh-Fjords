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
LSM_CARRA = xr.open_dataset(folder_CARRA+ "Land_sea_mask_CARRA.nc")
# LSM_RACMO = xr.open_dataset(folder_RACMO+ "FGRN055/Downscaling_GR_500m/GrIS_topo_icemask_lsm_lon_lat_0.5km.nc")
LSM_RACMO = xr.open_mfdataset("/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR_500m/GrIS_topo_icemask_lsm_lon_lat_0.5km.nc")


#%% Define data that will be analyzed
# List of years
years = [2023]

# List of months
#months = [1,2,3,4,5,6,7,8,9,10,11,12]


variables_R = ['LSM']
# List of CARRA variables
variables_C = ['precip2009_2023_yearly']

# Map of variable name for every dataset
variable_mapping = {
    #  'LWD': {'c': 'strd', 'r': 'lwsd'},
    #  'lat': {'c': 'latitude', 'r': 'lat'},
    #  'lon': {'c': 'longitude', 'r': 'lon'},
     'precip2009_2023_yearly': {'c': 'tp', 'r': 'precip'},
    #  'lsm': {'c': 'lsm', 'r': 'LSM'},
}

variable_location = {
    'precip2009_2023_yearly' :  {'in':    folder_CARRA + "Yearly/total_precipitation.CARRA.west_domain.forecast.Yearly.2009-2023.nc" , 
    'out': folder_CARRA + "Yearly/RACMOgrid/total_precipitation.CARRA.west_domain.2009-2023.0.5km.YY.nc" }
    }


# Open a dataset with the RACMO GRID
racmo_data_dict = {}
racmo_data_dict[f'{variable}'] = LSM_RACMO
print('Opened RACMO ------>', variable)


#%% Open CARRA data and store in CARRA dictionary with year, month, variable as key


# Initialize an empty dictionary to store all CARRA data
carra_data_dict = {}

for variable in variables_C:
    # Open CARRA data for the current variable
    carra_data_dict[f'{variable}'] = xr.open_dataset(variable_location[f'{variable}']['in'])
    print('Opened CARRA ------>', variable)




#%% Plot difference between CARRA and interpolated RACMO
   
for variable in variables_C:
    ofile = variable_location[f'{variable}']['out']
    print('----->' + ofile)
    # check if data is not already available
    if os.path.isfile(ofile): 
        print('file already exists, moving on...')
        continue

    ### Get data ###
    # Extract the carra and racmo data
    data_racmo = LSM_RACMO
    data_carra = carra_data_dict[f'{variable}']

    da_racmo = LSM_RACMO["LSM"]
    da_carra_interp = xr.DataArray(
                        coords={
                            'time': data_carra['time'],
                            'y': da_racmo['y'],
                            'x': da_racmo['x']
                        },
                        dims=['time', 'y', 'x']
                    )

    print("Preparing grid points...")
    # Prepare grid points        
    points_carra = (LSM_CARRA['longitude'].values.ravel()-360, LSM_CARRA['latitude'].values.ravel())

    for t in tqdm(range(len(data_carra['time']))):
    # for t in range(3):
        values_carra = getattr(carra_data_dict[f'{variable}'],variable_mapping[f'{variable}']['c']).isel(time=t).values.ravel()#*100
    
        print("Interpolating...")
        # Interpolate CARRA data to downscaled RACMO grid
        interp_carra_values = griddata(points_carra, values_carra, (data_racmo.lon, data_racmo.lat), method='nearest') #method='linear')
        print(t)
        da_carra_interp.values[t, :, :] = interp_carra_values[ :, :]
        print( da_carra_interp.values[t, :, :])

    # set all attributes of the interpolated data array
    da_carra_interp.attrs = getattr(carra_data_dict[variable],variable_mapping[variable]['c']).attrs
    ds_carra_interp = da_carra_interp.to_dataset(name=variable_mapping[variable]['r'])
    ds_carra_interp.attrs = carra_data_dict[variable].attrs
    ds_carra_interp.to_netcdf(ofile)


    
print("Plotting...")
plt.imshow(interp_carra_values)



# %%
