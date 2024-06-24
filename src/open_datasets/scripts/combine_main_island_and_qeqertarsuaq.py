import xarray as xr
from paths import *

base_folder = "/Volumes/imau02/rapid/Anneke/"
# base_folder = "/science/projects/imau02/rapid/Anneke/"
ofolder = pathDataTemp

# Load the datasets for qeqertarsuaq and fjords of main island
ds_qqt = xr.open_dataset(base_folder + "RACMO2.3p2/FGRN055/Downscaling_GR/fjords_Qeqertarsuaq.nc")
masks1k = xr.open_dataset(base_folder + "RACMO2.3p2/FGRN055/Downscaling_GR/masks1k_copy.nc")

# Create a copy of the 'fjords' variable
masks1k['fjords_no_qqt'] = masks1k['fjords'].copy(deep=True)

# Add the 'fjords' variable from masks1k with 'fjords_Qeqertarsuaq' from ds_qqt
masks1k['fjords'] = masks1k['fjords'] + ds_qqt['fjords_Qeqertarsuaq']

# Ensure no values are above 1
masks1k['fjords'] = masks1k['fjords'].where(masks1k['fjords'] <= 1, 1)


# Iterate over all data variables in the dataset
for var in masks1k.data_vars:
    # Check if both _FillValue and missing_value attributes exist
    if '_FillValue' in masks1k[var].encoding and 'missing_value' in masks1k[var].encoding:
        # Set the _FillValue and missing_value attributes to the same value
        masks1k[var].encoding['_FillValue'] = masks1k[var].encoding['missing_value']

# Save the dataset to a netCDF file
masks1k.to_netcdf(pathDataTemp + "masks1k.nc")

