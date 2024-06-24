import os
import xarray as xr
import numpy as np
from tqdm import tqdm
import glob


pathrapid = "/Volumes/imau01/rapid/RACMO2.3p2/FGRN055"
variable = "snowfrac"
time_resolution = "Monthly"
downscaling_type = "Downscaling_GR"
pathmonthly = os.path.join(pathrapid, downscaling_type, time_resolution)
pattern = os.path.join(pathmonthly, "**", "*.[ng]z")
files = glob.glob(pattern, recursive=True)


pathMonthlylocal = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/Monthly"
pathYearlyLocal = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/Annual/"
ds_yearly = []
ds_precip_month = []
all_years = []
fnames = []

year_int = 2000
while year_int  < 2022:

    year = str(year_int)
    all_years.append(year_int)
    print(f"========={year}=========")

    file_export_name = os.path.join(pathMonthlylocal, year, f"snowfrac.{year}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.nc")
    if os.path.exists(file_export_name):
        ds_yearly.append(xr.open_dataset(file_export_name))
    else:
        print(f"file does not exist: {file_export_name} ")
        year_int += 1
        continue
    print("loaded snowfrac", year)

     # open all files of a specific year and variabl
    files_precip_year = [file for file in files if "precip" in file and year in file]
    ds_precip_month.append(xr.open_dataset(files_precip_year[0]))
    print("loaded precip", year)

    if len(ds_yearly) == 2:
        ds_snowfrac = xr.concat(ds_yearly, dim="time")
        ds_precip = xr.concat(ds_precip_month, dim="time")
        ds_precip_snow_month = (ds_snowfrac["snowfrac"] * ds_precip["precipcorr"])
        ds_precip_snow_year = ds_precip_snow_month.groupby("time.year").sum(dim="time")
       
        ds_precip_year = ds_precip["precipcorr"].groupby("time.year").sum(dim="time")
        ds_snowfrac_year = (ds_precip_snow_year / ds_precip_year).to_dataset(name="snowfrac")
        fname = f"subwet_snowfrac.{year_int-1}-{year_int}.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.AN.nc"
        fnames.append(fname)
        ds_snowfrac_year.to_netcdf(fname)
        ds_yearly = []
        ds_precip_month = []
        print(f"saved {year_int-1}-{year_int}")

    year_int += 1
    # close all files


ds_annual = []
for file in fnames:
    print(file)
    ds = xr.open_dataset(file)
    ds_annual.append(ds)
    ds.close()
ds_annual = xr.concat(ds_annual, dim="year")
# rename year to time
ds_annual = ds_annual.rename({"year": "time"})

# assign attributes
ds_annual.attrs["units"] = "fraction"
ds_annual.attrs["long_name"] = "Annual 1km Topography corrected snowfrac"
ds_annual.attrs["title"] = "Annual mean snowfrac field at 1km (RACMO2.3p2 FGRN055 ERA)"
ds_annual.attrs["actual_range"] = "0, 1"

ds_annual["snowfrac"].attrs["units"] = "mm w.e. per year"
ds_annual["snowfrac"].attrs["long_name"] = "1km Topography corrected snow fraction"
ds_annual["snowfrac"].attrs["standard_name"] = "1km_Topography_snowfrac"
ds_annual["snowfrac"].attrs["actual_range"] = "[0, 1]"



ds_annual.to_netcdf(pathYearlyLocal+"/snowfrac.2000-2021.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.AN.nc")
print("saved snowfrac.2000-2021.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.AN.nc")

# delete all files in fnames
for file in fnames:
    os.remove(file)

print("deleted all files in fnames")