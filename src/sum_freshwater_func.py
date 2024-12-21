# %%
from scripts.open_preprocess_racmo import *
from scripts.paths import *
from scripts.dicts import *
from scripts.masks import mask_sections

import numpy as np
import gzip
import shutil
import os
import xarray as xr

import pandas as pd

# %%


# =====================================
# open all mask files
# =====================================
print("Import and make masks")


pathDataRaw = pathGithubFolder + "data/raw/"
fpath_adj_sect = f"{pathDataTemp}adjusted_section_numbers_slater.nc"
fpath_masks1k = f"{pathDataTemp}masks1k.nc"
folder_MARRACMO1km = f"{pathDataRaw}liquid/MAR_RACMO_Annual/"

ds_masks1k = xr.open_dataset(fpath_masks1k)
ds_adj_sect = xr.open_dataset(fpath_adj_sect)

dsmask_sections = xr.open_mfdataset(mask_sections)
# masks1k = xr.open_dataset("/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/masks1k.nc")
# fname = (
#     pathIMAU02
#     + "RACMO2.3p2/FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_may24_copy.nc"
# )
# mask = xr.open_dataset(fname)
masks1k = ds_masks1k
mask_tundra = masks1k["LSM"].where(masks1k["Promicemask"] == 0)
dsmask_sections["Tundra_sections"] = dsmask_sections["section_numbers_adjusted"].where(mask_tundra.values == 1)
mask_tundra_sections = dsmask_sections["Tundra_sections"]
mask_tundra_sections = mask_tundra_sections.rename({"x": "rlon", "y": "rlat"})


#  set basic parameters
dict_coords = {}
dict_vars = {}
time_resolution = "Annual"
years = ["1940", "2023"]
spatial_resolution = "1k"


# ===================================================
# Import all FW input source variables for ANNUAL plot
# ====================================================
print("Import all FW input source variables for ANNUAL plot")

dsRunoffTundra_mmyear = read_RACMO("1k_reproject", time_resolution, years=years, variable="runoff")
dsRunoffTundra_mmyear = dsRunoffTundra_mmyear.rename({"rlon": "x", "rlat": "y"})
dsRunoffTundra = volume(
    mask_data(dsRunoffTundra_mmyear["runoff"], "LSM Tundra", spatial_resolution),
    spatial_resolution,
)

dict_coords["Tundra"] = dict({"x": "x", "y": "y"})
dsRunoffTundraSum = dsRunoffTundra.sum(dim=["x", "y"]).resample(time="YS").sum()


# dsRunoff = read_RACMO(spatial_resolution, time_resolution, years=years, variable="runoff")
# # dsRunoff['time'] =  pd.to_datetime("1958-01-01") + pd.to_timedelta(dsRunoff.time.values*365, 'D')

# dsRunoff_RACMO_1k_YY_GrIS = (
#     mask_data(dsRunoff["runoffcorr"], "PROMICE_GrIS", spatial_resolution).resample(time="YS").sum()
# )
# dsRunoff_RACMO_1k_YY_GIC = (
#     mask_data(dsRunoff["runoffcorr"], "PROMICE_Ice_caps", spatial_resolution).resample(time="YS").sum()
# )
# dsRunoff_RACMO_1k_YY_GrIS_sum = dsRunoff_RACMO_1k_YY_GrIS.sum(dim=["x", "y"]) / 1e6
# dsRunoff_RACMO_1k_YY_GIC_sum = dsRunoff_RACMO_1k_YY_GIC.sum(dim=["x", "y"]) / 1e6

# dsRunoff_RACMO_1k_YY_GrIS_sections = (
#     dsRunoff_RACMO_1k_YY_GrIS.groupby(dsmask_sections["section_numbers_adjusted"]).sum() / 1e6
# )
# dsRunoff_RACMO_1k_YY_GIC_sections = (
#     dsRunoff_RACMO_1k_YY_GIC.groupby(dsmask_sections["section_numbers_adjusted"]).sum() / 1e6
# )


# =====================================
# Import SOLID ice DISCHARGE
# =====================================


# %% SOLID DISCHARGE
dfGISDMankoff = pd.read_csv(path_Mankoff2020Solid + "GIS_D.csv", index_col=0, parse_dates=True)
dfErrorGISDMankoff = pd.read_csv(path_Mankoff2020Solid + "GIS_err.csv", index_col=0, parse_dates=True)
dfGIScoverageMankoff = pd.read_csv(path_Mankoff2020Solid + "GIS_coverage.csv", index_col=0, parse_dates=True)
filterCov = (dfGIScoverageMankoff < 0.5).values
dfGISDMankoff = dfGISDMankoff.groupby(dfGISDMankoff.index.year).mean()
dfErrorGISDMankoff = dfErrorGISDMankoff.groupby(dfErrorGISDMankoff.index.year).mean()


# =====================================
# Import precipitation
# =====================================

# RACMO precipitation

file_annual_fjord = f"{pathDataTemp}RACMO2.3p2_1km_precip_fjords_Annual_1958_2023.nc"
dsPrecipFjordsVol = xr.open_dataset(file_annual_fjord)
dsPrecipFjordsVolSum = dsPrecipFjordsVol.sum(dim=["x", "y"]).resample(time="YS").sum()

# CARRA precipitation
# ds_precip_carra_1991_2008_sum = xr.open_mfdataset(
#     pathIMAU02
#     + "CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2008.1km.YY.fjords_only.sum_per_basin.nc"
# )
# ds_precip_carra_2009_2023_sum = xr.open_mfdataset(
#     pathIMAU02
#     + "CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.2009-2023.1km.YY.fjords_only.sum_per_basin.nc"
# )
# # concatenate
# ds_precip_carra_1991_2023_sum = xr.concat(
#     [ds_precip_carra_1991_2008_sum, ds_precip_carra_2009_2023_sum], dim="time"
# )
# dsPrecipFjordsCARRA_Annual_Sum = (
#     (ds_precip_carra_1991_2023_sum.sum(dim="section_numbers_adjusted") / 1e6)
#     .resample(time="YS")
#     .sum()
# )


# CARRA precipitation
ds_precip_carra_1991_2023_sum = (
    xr.open_mfdataset(
        [
            pathDataTemp
            + "Summed_grouped/"
            + "total_precipitation.CARRA.west_domain.1991-2008.1km.YY.fjords_only.sum_per_basin.nc",
            pathDataTemp
            + "Summed_grouped/"
            + "total_precipitation.CARRA.west_domain.2009-2023.1km.YY.fjords_only.sum_per_basin.nc",
        ]
    )
    / 1e6
)
ds_precip_carra_1991_2023_sum.resample(time="YS").sum()

# =====================================
# Prepare files for combining in one dataframe
# =====================================

dfRunoffTundra = dsRunoffTundraSum.squeeze().to_dataframe(name="Runoff Tundra").drop(columns="height")

dfPrecipFjordsVol = dsPrecipFjordsVolSum["precipcorr"].to_dataframe(name="Precipitation Fjords")
dfSolidMankoff = dfGISDMankoff.copy(deep=True)

dfRunoffTundra.index = dfRunoffTundra.index.year
dfPrecipFjordsVol = dfPrecipFjordsVol.groupby(dfPrecipFjordsVol.index.year).mean()

# dfRunoff_RACMO_1k_YY_GIC_sum = dsRunoff_RACMO_1k_YY_GIC_sum.squeeze().to_dataframe(name="Runoff Ice Caps (RACMO 1km)")
# dfRunoff_RACMO_1k_YY_GrIS_sum = dsRunoff_RACMO_1k_YY_GrIS_sum.squeeze().to_dataframe(name="Runoff GrIS (RACMO 1km)")
dfPrecipFjordsCARRA = (
    dsPrecipFjordsCARRA_Annual_Sum["precip"].squeeze().to_dataframe(name="Precipitation Fjords (CARRA)")
)
dfPrecipFjordsCARRA.index = dfPrecipFjordsCARRA.index.year
# dfRunoff_RACMO_1k_YY_GIC_sum.index = dfRunoff_RACMO_1k_YY_GIC_sum.index.year
# dfRunoff_RACMO_1k_YY_GrIS_sum.index = dfRunoff_RACMO_1k_YY_GrIS_sum.index.year


# df_sum_GIS_55 = pd.concat(
#     [
#         dfSolidMankoff,
#         dfRunoff_RACMO_1k_YY_GIC_sum,
#         dfRunoff_RACMO_1k_YY_GrIS_sum,
#         dfRunoffTundra,
#         dfPrecipFjordsVol,
#         dfPrecipFjordsCARRA,
#         # dfRunoffRACMOGrIS,
#         # dfRunoffRACMOGIC,
#     ],
#     axis=1,
# ).sort_index()

# df_sum_GIS_55.index = pd.to_datetime(df_sum_GIS_55.index, format="%Y")

# filterNanMonths = df_sum_GIS_55.isna().sum(axis=1) > 0  # find months with NaN values in any of the variables
# df_sum_GIS_55Relative = (df_sum_GIS_55.iloc[:,:5].mask(filterNanMonths).T / df_sum_GIS_55.iloc[:,:5].sum(axis=1)).T


# %% =====================
# PER SECTOR and per month
# =======================
def kgperm2_to_Gt(ds):
    ds = ds / 1e6
    ds.attrs["units"] = "Gt"
    return ds


# %% loading datasets

# ======
# Load GrIS
# ======
dsRunoffIceSectormm = xr.open_dataset(
    pathDataTemp + "Summed_grouped/" + "runoff_GrIS.1990-2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.GrIS.section_sum.nc"
)["runoffcorr"]
dsRunoffIceSector = kgperm2_to_Gt(dsRunoffIceSectormm)
dfRunoffIceSector = (
    dsRunoffIceSector.to_dataframe(name="Liquid Runoff Ice Sheet")
    .reset_index()
    .set_index("time")
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)

dfRunoffIceSector.columns = dfRunoffIceSector.columns.get_level_values(1)
dsRunoffIceCapSectormm = xr.open_dataset(
    pathDataTemp + "Summed_grouped/" + "runoff_GIC.1990-2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.GIC.section_sum.nc"
)["runoffcorr"]

dsRunoffIceCapSector = kgperm2_to_Gt(dsRunoffIceCapSectormm)
dfRunoffIceCapSector = (
    dsRunoffIceCapSector.to_dataframe(name="Liquid Runoff Ice Sheet")
    .reset_index()
    .set_index("time")
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)
dfRunoffIceCapSector.columns = dfRunoffIceCapSector.columns.get_level_values(1)

# ======
# Load Tundra
# ======
time_resolution = "Monthly"
years = list(range(2012, 2023))
dsRunoffTundra_mmyear = read_RACMO("1k_reproject", time_resolution, years=years, variable="runoff")

# rename x and y to rlon and rlat
dsRunoffTundra_mmyear["Tundra_basins"] = mask_tundra_sections
path_monthly_racmo = pathAnnekeFolderIMAU02 + "/Downscaling_GR/Monthly/"
fname_runofftundraMM = f"runoff_tundra.{years[0]}-{years[1]}.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.tundra.section_sum.nc"

if os.path.isfile(path_monthly_racmo + fname_runofftundraMM):
    dsRunoffTundra_mmyear_SECTOR = xr.open_dataset(path_monthly_racmo + fname_runofftundraMM)
else:
    dsRunoffTundra_mmyear_SECTOR = (
        dsRunoffTundra_mmyear["runoff"].groupby(dsRunoffTundra_mmyear["Tundra_basins"]).sum().drop("height")
    )
    dsRunoffTundra_mmyear_SECTOR.to_netcdf(path_monthly_racmo + fname_runofftundraMM)
# add 2023
ds_runoff_tundra2023 = xr.open_dataset(
    f"{pathDataTemp}runoff_tundra.2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.tundra.nc", engine="netcdf4"
)
ds_runoff_tundra2023["Tundra_basins"] = ds_masks1k["LSM"]
ds_runoff_tundra2023["Tundra_basins"].values = ds_adj_sect["section_numbers_adjusted"].values
ds_runoff_tundra2023_Basins = (
    ds_runoff_tundra2023.groupby("Tundra_basins").sum().rename({"Liquid Runoff Tundra": "runoff"})
)
dsRunoffTundra_mmyear_SECTOR = (
    xr.concat([ds_runoff_tundra2023_Basins, dsRunoffTundra_mmyear_SECTOR], dim="time")
    .sortby("time")
    .resample(time="MS")
    .mean()
)

dfRunoffTundraSector = (
    kgperm2_to_Gt(dsRunoffTundra_mmyear_SECTOR["runoff"])
    .to_dataframe()
    .reset_index()
    .drop(columns="height")
    .set_index("time")
    .pivot(columns="Tundra_basins")
    .rename(columns=dict_sections)
)

dfRunoffTundraSector.columns = dfRunoffTundraSector.columns.get_level_values(1)


# =====
# Load Solid ice discharge
# =====

#  from https://essd.copernicus.org/articles/12/1367/2020/
dfSectionD = pd.read_csv(path_Mankoff2020Solid_adjusted + "section_D.csv", index_col=0, parse_dates=True)
dfErrorSectionDMankoff = pd.read_csv(path_Mankoff2020Solid_adjusted + "section_err.csv", index_col=0, parse_dates=True)
dfSectionCovManoff = pd.read_csv(path_Mankoff2020Solid_adjusted + "section_coverage.csv", index_col=0, parse_dates=True)
# filterCovReg =  (dfSectionCovManoff < 0.5).values

dfSectionD = dfSectionD / 12
dfErrorSectionDMankoff = dfErrorSectionDMankoff / 12
dfSectionD.columns.name = "Basins"


dsPrecipFjords = xr.open_mfdataset(pathDataTemp + "RACMO2.3p2_1km_precip_fjords_Annual_1958_2023.nc")
dsPrecipFjordsAnnualRACMO_58_23 = dsPrecipFjords.copy(deep=True)

dsmask_sections = xr.open_mfdataset(mask_sections)
dsPrecipFjordsSectormm_sum = dsPrecipFjords.groupby(dsmask_sections["section_numbers_adjusted"]).sum()
dsPrecipFjordsSectormm_sum = dsPrecipFjordsSectormm_sum.resample(time="YS").sum()
# %% Open precipitation monthly
# monthly_precip_cara1 = (
#     pathIMAU02
#     + "CARRA/Monthly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.2009-2023.1km.MM.fjords_only.sum_per_basin.nc"
# )
# ds_precip_carra_2009_2023_month = xr.open_dataset(monthly_precip_cara1) / 1e6

# ds_precip_carra_1991_2008_month = (
#     xr.open_dataset(
#         pathIMAU02
#         + "CARRA/Monthly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2008.1km.MM.fjords_only.sum_per_basin.nc"
#     )
#     / 1e6
# )
# ds_precip_carra_1991_2023_month = xr.concat(
#     [ds_precip_carra_1991_2008_month, ds_precip_carra_2009_2023_month], dim="time"
# )

# path_monthly_precip_carra_1991_2008 = (
#     pathIMAU02 + "CARRA/Monthly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2008.1km.MM.fjords_only.sum_per_basin.nc")
# path_monthly_precip_carra_2009_2023 = (pathIMAU02 + "CARRA/Monthly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.2009-2023.1km.MM.fjords_only.sum_per_basin.nc")

# ds_precip_carra_1991_2008_month = xr.open_dataset(path_monthly_precip_carra_1991_2008)
# ds_precip_carra_2009_2023_month = xr.open_dataset(path_monthly_precip_carra_2009_2023)
# ds_precip_carra_1991_2023_month = xr.concat([ds_precip_carra_1991_2008_month, ds_precip_carra_2009_2023_month], dim="time") /1e6

ds_precip_racmo_1990_2023_sum_monthly = xr.open_mfdataset(
    pathDataTemp
    + "Summed_grouped/"
    + "precip.1990-2023.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.fjords_only.sum_per_basin.nc",
    decode_times=False,
)
ds_precip_racmo_1990_2023_sum_monthly["time"] = pd.to_datetime("1990-01-01") + pd.to_timedelta(
    ds_precip_racmo_1990_2023_sum_monthly.time, "D"
)

dfPrecipFjordsSector_RACMO_year = (
    dsPrecipFjordsSectormm_sum.to_dataframe()
    .reset_index()
    .set_index("time")
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)
dfPrecipFjordsSector_RACMO_year.columns = dfPrecipFjordsSector_RACMO_year.columns.get_level_values(1)

# dfPrecipFjordsSector_CARRA_year = (
#     ds_precip_carra_1991_2023_sum.to_dataframe()
#     .reset_index()
#     .set_index("time")
#     .pivot(columns="section_numbers_adjusted")
#     .rename(columns=dict_sections)
# )
# dfPrecipFjordsSector_CARRA_year.columns = dfPrecipFjordsSector_CARRA_year.columns.get_level_values(
#     1
# )

dfPrecipFjordsSector_RACMO_month = (
    ds_precip_racmo_1990_2023_sum_monthly.to_dataframe()
    .reset_index()
    .set_index(["time", "section_numbers_adjusted"])[["precipcorr"]]
    .reset_index()
    .set_index("time")
    .drop_duplicates()
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)
dfPrecipFjordsSector_RACMO_month.columns = dfPrecipFjordsSector_RACMO_month.columns.get_level_values(1)

dfPrecipFjordsSector_CARRA_month = (
    ds_precip_carra_1991_2023_month.to_dataframe()
    .reset_index()
    .set_index("time")
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)
dfPrecipFjordsSector_CARRA_month.columns = dfPrecipFjordsSector_CARRA_month.columns.get_level_values(1)
# %% =====
# Load Basal melt
# =====

df_Basal_Basin_monthly = pd.read_csv(pathDataTemp + "Basal_melt/Basal_basin_monthly.csv", index_col=0)
start_year = 1980
end_year = 2023
date_range = pd.date_range(start=f"{start_year}-01-01", end=f"{end_year}-12-31", freq="MS")
repeated_cycle = np.tile(np.squeeze(df_Basal_Basin_monthly.values), (end_year - start_year + 1, 1))
df_time_series_Basal = pd.DataFrame(data=repeated_cycle, index=date_range, columns=df_Basal_Basin_monthly.columns)
ds_time_series_Basal = (
    df_time_series_Basal.stack()
    .to_xarray()
    .to_dataset(name="Basal melt")
    .rename({"level_0": "time", "level_1": "Basins"})
)

# %%  =========================================
# Combine all monthly sector data in one dataset
# =========================================
open_previous = False

if open_previous:
    # dsSectorSum = xr.open_dataset(pathDataTemp + "RACMO2.3p2_1k_sector_sum_2024_06_12.nc")
    dsSectorSum = xr.open_dataset(pathDataTemp + "Sector_sum_2024_10_23.nc")

else:

    dsSectorSum = xr.merge(
        [
            dfRunoffIceSector.stack()
            .to_xarray()
            .to_dataset(name="Liquid Runoff Ice Sheet")
            .rename({"section_numbers_adjusted": "Basins"}),
            dfRunoffIceCapSector.stack()
            .to_xarray()
            .to_dataset(name="Liquid Runoff Ice Caps")
            .rename({"section_numbers_adjusted": "Basins"}),
            dfRunoffTundraSector.stack()
            .to_xarray()
            .to_dataset(name="Liquid Runoff Tundra")
            .rename({"Tundra_basins": "Basins"}),
            dfPrecipFjordsSector_RACMO_month.stack()
            .to_xarray()
            .to_dataset(name="Precipitation Fjords")
            .rename({"section_numbers_adjusted": "Basins"})
            / 1e6,
            ds_time_series_Basal,
        ]
    )

    def compute_weighted_monthly_mean(data):
        # Calculate the weights as the difference between each date and the previous date

        diff_next = data.Date.diff("Date").shift(Date=-1).dt.days
        diff_prev = data.Date.diff("Date").dt.days
        diff = (diff_next + diff_prev) / 2
        weights = diff
        total_weights = weights.resample(Date="1ME").sum()

        # Compute the weighted mean
        weighted_data = data * weights
        weighted_sum = weighted_data.resample(Date="1ME").sum()
        weighted_mean = weighted_sum / total_weights
        weighted_mean["Date"] = weighted_mean.Date - pd.Timedelta(15, unit="D")
        return weighted_mean

    dsSectorSum["Solid Ice discharge (weighted mean)"] = (
        compute_weighted_monthly_mean(dfSectionD.stack().to_xarray())
        .interp(Date=dsSectorSum.time, method="nearest")
        .drop("Date")
    )

    #
    # export dsSectorSum to csv
    start = dsSectorSum.time[0].dt.strftime("%Y").values
    end = dsSectorSum.time[-1].dt.strftime("%Y").values
    dsSectorSum.to_netcdf(pathDataTemp + "RACMO2.3p2_1k_sector_sum_2024_10_22.nc")

dsSectorSum = dsSectorSum.where(dsSectorSum != 0)
# (dsSectorSum.resample(time="YS").mean() * 12).sum(
#     dim="Basins"
# ).to_dataframe().reset_index().set_index("time").to_csv(
#     pathDataTemp + f"Greenland_Sum_2024_06_12.csv"
# )

# Sort dataset based on mean
means = {var: dsSectorSum[var].mean().item() for var in dsSectorSum.data_vars}
variables_sorted = (
    pd.DataFrame(list(means.items()), columns=["Variable", "Mean"])
    .sort_values(by="Mean", ascending=False)["Variable"]
    .values
)
dsSectorSum = dsSectorSum[variables_sorted]


try:
    dsMonthlyGr = dsSectorSum.copy(deep=True).drop("winter_year").drop("summer_year")
except:
    dsMonthlyGr = dsSectorSum.copy(deep=True)

dsMonthlyGr = dsMonthlyGr.where(dsMonthlyGr > 0)
dsMonthlyGr = dsMonthlyGr.rename(dict_shorter_name)
dsMonthlyGr = dsMonthlyGr[col_order_rel]
dsMonthlyGr = dsMonthlyGr.resample(time="MS").mean()


# Calculate the mean and standard deviation for each month

period1 = {"start": "1990", "end": "2004"}
period2 = {"start": "2005", "end": "2023"}


def calculate_seasonal_statistics(ds, period, dim="Basins"):
    """
    Calculate the mean and standard deviation of seasonal data for a given period.

    Parameters:
    ds (xarray.Dataset): The dataset to calculate statistics on.
    period (dict): A dictionary with 'start' and 'end' keys defining the time period.
    dim (str): The dimension to sum over before calculating statistics.

    Returns:
    tuple: A tuple containing the mean and standard deviation of the seasonal data.
    """
    data_mean_seasonal = ds.sel(time=slice(period["start"], period["end"])).sum(dim=dim).groupby("time.month").mean()
    data_std_seasonal = ds.sel(time=slice(period["start"], period["end"])).sum(dim=dim).groupby("time.month").std()
    return data_mean_seasonal, data_std_seasonal


data_mean_seasonal_period2, data_std_seasonal_period2 = calculate_seasonal_statistics(dsMonthlyGr, period2)
data_mean_seasonal_period1, data_std_seasonal_period1 = calculate_seasonal_statistics(dsMonthlyGr, period1)


#  Recalculate solid ice discharge for different periods
dfGIS_period1 = pd.read_csv(path_Mankoff2020Solid + "GIS_D.csv", index_col=0, parse_dates=True)[
    period1["start"]: period1["end"]
]
dfGIS_period2 = pd.read_csv(path_Mankoff2020Solid + "GIS_D.csv", index_col=0, parse_dates=True)[
    period2["start"]: period2["end"]
]
data_mean_seasonal_period1["Solid Ice Discharge"].values = (
    np.squeeze(dfGIS_period1.groupby(dfGIS_period1.index.month).mean().values) / 12
)
data_mean_seasonal_period2["Solid Ice Discharge"].values = (
    np.squeeze(dfGIS_period2.groupby(dfGIS_period2.index.month).mean().values) / 12
)
data_std_seasonal_period2["Solid Ice Discharge"].values = (
    np.squeeze(dfGIS_period2.groupby(dfGIS_period2.index.month).std().values) / 12
)
data_std_seasonal_period2["Solid Ice Discharge"].values = (
    np.squeeze(dfGIS_period2.groupby(dfGIS_period2.index.month).std().values) / 12
)


data_mean_seasonal_good_data, data_std_seasonal_good_data = calculate_seasonal_statistics(
    dsMonthlyGr, {"start": "2009", "end": "2023"}
)


# %% MAR 1 KM RUNOFF


def convert_months_to_date(data_array, start_date):
    """
    This function converts the months since a start date to a date. For format dat_array as ds.time, start_date should 'YYYY-MM-DD'
    """
    array_dates = np.empty(len(data_array)).astype("datetime64[ns]")
    for i in range(len(data_array)):
        array_dates[i] = pd.to_datetime(start_date) + pd.DateOffset(
            months=np.floor(data_array[i]), days=((data_array[i] % 1) * 30)
        )
    return array_dates


def convert_years_to_date(data_array, start_date):
    """
    This function converts the years since a start date to a date. For format dat_array as ds.time, start_date should 'YYYY-MM-DD'
    """
    array_dates = np.empty(len(data_array)).astype("datetime64[ns]")
    for i in range(len(data_array)):
        array_dates[i] = pd.to_datetime(start_date) + pd.DateOffset(
            years=np.floor(data_array[i]), days=((data_array[i] % 1) * 365)
        )
    return array_dates


def open_compressed_xarray(file_path):
    """
    Opens a compressed netcdf file with xarray
    """
    # Path to the compressed file
    compressed_file = file_path

    # Path to the decompressed file
    decompressed_file = compressed_file.replace(".gz", "")

    # Decompress the file
    with gzip.open(compressed_file, "rb") as f_in:
        with open(decompressed_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Open the decompressed file with xarray
    ds = xr.open_dataset(decompressed_file, engine="netcdf4", decode_times=False)

    # remove the decompressed file
    os.remove(decompressed_file)
    return ds


def process_runoff_data(file_path, start_date, convert_func, ds_adj_sect, ds_masks1k):
    ds_run = open_compressed_xarray(file_path)
    ds_run["time_since_start"] = ds_run.time
    ds_run["time"] = convert_func(ds_run["time_since_start"].values, start_date)  # Convert timedelta to datetime
    ds_run_mean = ds_run.mean(dim=["time"])  # Get masks in the same xarray coordinate system

    # ds_run_mean['section_numbers_adjusted'] = add_variable(ds_run_mean, ds_adj_sect['section_numbers_adjusted'], 'runoffcorr')
    ds_run_mean["section_numbers_adjusted"] = ds_run_mean["runoffcorr"].copy(deep=True)
    ds_run_mean["section_numbers_adjusted"].values = ds_adj_sect["section_numbers_adjusted"].values

    for var in ds_masks1k.data_vars:
        ds_run_mean[var] = ds_run_mean["runoffcorr"].copy(deep=True)
        ds_run_mean[var].values = ds_masks1k[var].values

    # Calculate runoff for each basin yearly in Gt
    ds_run_GrIS_basin = (
        ds_run["runoffcorr"].where(ds_run_mean["GrIS"] == 1).groupby(ds_run_mean["section_numbers_adjusted"]).sum()
        / 1e6
    )
    ds_run_GIC_basin = (
        ds_run["runoffcorr"].where(ds_run_mean["GIC"] == 1).groupby(ds_run_mean["section_numbers_adjusted"]).sum() / 1e6
    )

    return ds_run_GrIS_basin, ds_run_GIC_basin


# Process MAR data
ds_run_MAR_GrIS_basin, ds_run_MAR_GIC_basin = process_runoff_data(
    folder_MARRACMO1km + "runoff.1940-2023.MARv3.14-ERA5.1km.YY.nc.gz",
    "1940-01-15",
    convert_months_to_date,
    ds_adj_sect,
    ds_masks1k,
)

# Update dsSectorSumMAR with MAR data
dsSectorSumMAR = dsSectorSum.copy(deep=True).resample(time="YS").mean() * 12
dsSectorSumMAR["Liquid Runoff Ice Caps"] = (
    ds_run_MAR_GIC_basin.rename({"section_numbers_adjusted": "Basins"}).resample(time="YS").mean(skipna=True)
)
dsSectorSumMAR["Liquid Runoff Ice Sheet"] = (
    ds_run_MAR_GrIS_basin.rename({"section_numbers_adjusted": "Basins"}).resample(time="YS").mean(skipna=True)
)

# Map section numbers to names
mapped_values = np.array([dict_sections.get(val, val) for val in ds_run_MAR_GIC_basin.section_numbers_adjusted.values])
ds_run_MAR_GIC_basin["section_numbers_adjusted"] = xr.DataArray(mapped_values, dims="section_numbers_adjusted")
ds_run_MAR_GrIS_basin["section_numbers_adjusted"] = xr.DataArray(mapped_values, dims="section_numbers_adjusted")

# Process RACMO data
ds_run_RACMO_GrIS_basin, ds_run_RACMO_GIC_basin = process_runoff_data(
    folder_MARRACMO1km + "runoff.1958-2023.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.YY.nc.gz",
    "1958-01-15",
    convert_years_to_date,
    ds_adj_sect,
    ds_masks1k,
)
