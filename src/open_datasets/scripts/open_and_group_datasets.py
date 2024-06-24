
#%%
import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt


# Define paths and year range

base_imau1 = "/Volumes/imau01/rapid/"
base_imau2 = "/Volumes/imau02/rapid/"
base_runoff_path =base_imau1 + "RACMO2.3p2/FGRN055/Downscaling_GR/Monthly"
runoff_filename_template = "runoff.{year}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.nc.gz"
mask_path = base_imau1 +"RACMO2.3p2/FGRN055/Downscaling_GR/Icemask_Topo_Iceclasses_lon_lat_average_1km.nc"
adjusted_sections_path = base_imau2+ "Anneke/RACMO2.3p2/FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_may24_copy_copy.nc"
ofolder = base_imau2 + "Anneke/RACMO2.3p2/FGRN055/Downscaling_GR/Monthly"

years = range(1990, 2023 + 1)  # Including 2023

ofile_gic = os.path.join(ofolder, f"runoff_GIC.{years[0]}-{years[-1]}.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.GIC.section_sum.nc")
ofile_gris = os.path.join(ofolder, f"runoff_GrIS.{years[0]}-{years[-1]}.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.GrIS.section_sum.nc")




# List to hold individual year datasets
datasets_gic = []
datasets_gris = []

# Load mask datasets
print("Loading mask datasets...")
mask = xr.open_dataset(mask_path)
mask_sections = xr.open_dataset(adjusted_sections_path)
print("Mask datasets loaded.")

# Load and concatenate runoff datasets for all years
print("Loading and concatenating runoff datasets...")
for year in years:
    runoff_path = os.path.join(base_runoff_path, str(year), runoff_filename_template.format(year=year))
    print(f"Loading dataset for year {year}:\n ====> {runoff_path}")
    dsrunoff = xr.open_dataset(runoff_path, decode_times=True)


    dsrunoff['GIC'] = dsrunoff.LAT
    dsrunoff['GIC'].values = mask['GIC'].values
    dsrunoff['GrIS'] = dsrunoff.LAT
    dsrunoff['GrIS'].values = mask['GrIS'].values
    dsrunoff['section_numbers_adjusted'] = dsrunoff.LAT
    dsrunoff['section_numbers_adjusted'].values = mask_sections['section_numbers_adjusted'].values

    # Groupby and sum for GIC
    print("Processing GIC runoff data...")
    gic_grouped = dsrunoff.where(dsrunoff['GIC']).groupby(dsrunoff['section_numbers_adjusted']).sum()
    datasets_gic.append(gic_grouped)
    print("GIC runoff data processed.")

    # Groupby and sum for GrIS
    print("Processing GrIS runoff data...")
    gris_grouped = dsrunoff.where(dsrunoff['GrIS']).groupby(dsrunoff['section_numbers_adjusted']).sum()
    datasets_gris.append(gris_grouped)
    print("GrIS runoff data processed.")

# Concatenate all datasets along the time dimension
dsrunoffGIC = xr.concat(datasets_gic, dim='time')
dsrunoffGrIS = xr.concat(datasets_gris, dim='time')

print("All runoff datasets concatenated.")
# Save GIC runoff data to netCDF
print(f"Saving GIC runoff data to {ofile_gic}...")
dsrunoffGIC.to_netcdf(ofile_gic, format='NETCDF4')
print(f"GIC runoff data saved to {ofile_gic}")


# Save GrIS runoff data to netCDF
print(f"Saving GrIS runoff data to {ofile_gris}...")
dsrunoffGrIS.to_netcdf(ofile_gris, format='NETCDF4')
print(f"GrIS runoff data saved to {ofile_gris}")



# %%
