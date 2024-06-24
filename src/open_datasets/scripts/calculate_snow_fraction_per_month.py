import os
import xarray as xr
import numpy as np
from tqdm import tqdm
import glob


pathrapid = "/Volumes/imau01/rapid/RACMO2.3p2/FGRN055"
variable = "snowfrac"
time_resolution = "Daily"
downscaling_type = "Downscaling_GR"
path_daily1k = os.path.join(pathrapid, downscaling_type, time_resolution)
pattern = os.path.join(path_daily1k, "**", "*.[ng]z")
files = glob.glob(pattern, recursive=True)


pathMonthlylocal = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/Monthly"
pathYearlyLocal = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/Annual/"
ds_yearly = []
all_years = []

for year_int in tqdm(np.arange(2001,2022)):
    year = str(year_int)
    all_years.append(year_int)
    print(f"========={year}=========")

    file_export_name = os.path.join(pathMonthlylocal, year, f"snowfrac.{year}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.nc")
    if os.path.exists(file_export_name):
        print(f"file {file_export_name} already exists")
        ds_yearly.append(xr.open_dataset(file_export_name))
        continue
    if year == "2001":
        print("skipping 2001")
        continue

    # open all files of a specific year and variable
    files_snowfrac_year = [file for file in files if "snowfrac" in file and year in file]
    ds_snowfrac = xr.open_mfdataset(files_snowfrac_year, chunks={"snowfraccorr": "auto"})
    print("loaded snowfrac", year)

    files_precip_year = [file for file in files if "precip" in file and year in file]
    ds_precip = xr.open_mfdataset(files_precip_year, chunks={"precip": "auto"})
    print("loaded precip", year)

    # calculate monthly mean snowfrac per month and precip_snow
    ds_precip_snow = (ds_snowfrac["snowfraccorr"] * ds_precip["precipcorr"]).groupby("time.month").sum(dim="time").to_dataset(name="precip_snow")  # sum over snow precipitation
    ds_precip_snow["snowfrac"] = ds_precip_snow["precip_snow"] / ds_precip["precipcorr"].groupby("time.month").sum(dim="time")
    time_monthly_date = ds_precip.time.groupby("time.month").mean(dim="time")
    ds_precip_snow["month"] = time_monthly_date
    ds_precip_snow = ds_precip_snow.rename({"month": "time"})

    # assign attributes
    ds_precip_snow.attrs = ds_snowfrac.attrs
    ds_precip_snow.attrs["title"] = "Monthly mean snowfrac field at 1km (RACMO2.3p2 FGRN055 ERA)"
    ds_precip_snow["snowfrac"].attrs = ds_snowfrac["snowfraccorr"].attrs
    ds_precip_snow["snowfrac"].attrs["units"] = "fraction"
    ds_precip_snow["precip_snow"].attrs = ds_precip["precipcorr"].attrs
    ds_precip_snow["precip_snow"].attrs["long_name"] = "1km Topography precipitation as snow"
    ds_precip_snow["precip_snow"].attrs["standard_name"] = "1km_Topography_precip_snow"
    ds_precip_snow["precip_snow"].attrs["units"] = "mm w.e. per month"

    print("computed monthly mean snowfrac", year)

    # yearly 
    ds_precip_snow_yearly = ds_precip_snow.groupby('time.year').sum(dim='time').copy(deep=True)
    ds_precip_snow_yearly['snowfrac'] = ds_precip_snow_yearly['precip_snow'] / ds_precip['precipcorr'].groupby('time.year').sum(dim='time')
    ds_precip_snow_yearly = ds_precip_snow_yearly.rename({"year": "time"})
    # assign attributes
    ds_precip_snow_yearly.attrs = ds_precip_snow.attrs
    ds_precip_snow_yearly.attrs["title"] = "Yearly mean snowfrac and snow accumulation field at 1km (RACMO2.3p2 FGRN055 ERA)"
    ds_precip_snow_yearly["snowfrac"].attrs = ds_precip_snow["snowfrac"].attrs
    ds_precip_snow_yearly["precip_snow"].attrs = ds_precip_snow["precip_snow"].attrs
    ds_yearly.append(ds_precip_snow_yearly)


    # export to netcdf
    ds_precip_snow = ds_precip_snow.drop_vars("precip_snow")
    encoding = {var: {'zlib': True, 'complevel': 8} for var in ds_precip_snow.data_vars}
    ds_precip_snow.to_netcdf(file_export_name,  format='NETCDF4', encoding=encoding)
    print("exported monthly mean snowfrac", year)

    # close all ds to free memory
    ds_snowfrac.close()
    ds_precip.close()
    ds_precip_snow.close()

print("now exporting annual")
ds_yearly = xr.concat(ds_yearly, dim="time")
print(ds_yearly)
encoding = {var: {'zlib': True, 'complevel': 8} for var in ds_yearly.data_vars}

path_ds_yearly = os.path.join(pathYearlyLocal, f"snowfrac.{all_years[0]}-{all_years[-1]}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.YY.nc")
ds_yearly.to_netcdf(path_ds_yearly, format='NETCDF4')#, encoding=encoding)
                              
