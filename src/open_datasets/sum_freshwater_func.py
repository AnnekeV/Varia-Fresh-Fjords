# %%
from scripts.open_preprocess_racmo import *
from scripts.paths import *
from scripts.dicts import *

import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import numpy as np


import pandas as pd

# %%
pathIMAU02 = "/Volumes/imau02/rapid/Anneke/"


dict_coords = {}
dict_vars = {}
time_resolution = "Annual"
spatial_resolution = "1k"
years = ["1940", "2024"]

if spatial_resolution == "1k":
    # dsRunoff = dsRunoff.rename({'runoffcorr': 'runoff'})
    dsRunoffTundra_mmyear = read_RACMO(
        "1k_reproject", time_resolution, years=years, variable="runoff"
    )
    dsRunoffTundra_mmyear = dsRunoffTundra_mmyear.rename({"rlon": "x", "rlat": "y"})
    dsRunoffTundra = volume(
        mask_data(dsRunoffTundra_mmyear["runoff"], "LSM Tundra", spatial_resolution),
        spatial_resolution,
    )
    dict_coords["Tundra"] = dict({"x": "x", "y": "y"})

dsRunoffTundraSum = dsRunoffTundra.sum(dim=["x", "y"]).resample(time="YS").sum()

mask_sections = (
    pathIMAU02 + "RACMO2.3p2/FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_may24.nc"
)
dsmask_sections = xr.open_mfdataset(mask_sections)
# 1km runoff
dsRunoff = read_RACMO(spatial_resolution, time_resolution, years=years, variable="runoff")
# dsRunoff['time'] =  pd.to_datetime("1958-01-01") + pd.to_timedelta(dsRunoff.time.values*365, 'D')

dsRunoff_RACMO_1k_YY_GrIS = (
    mask_data(dsRunoff["runoffcorr"], "PROMICE_GrIS", spatial_resolution).resample(time="YS").sum()
)
dsRunoff_RACMO_1k_YY_GIC = (
    mask_data(dsRunoff["runoffcorr"], "PROMICE_Ice_caps", spatial_resolution)
    .resample(time="YS")
    .sum()
)
dsRunoff_RACMO_1k_YY_GrIS_sum = dsRunoff_RACMO_1k_YY_GrIS.sum(dim=["x", "y"]) / 1e6
dsRunoff_RACMO_1k_YY_GIC_sum = dsRunoff_RACMO_1k_YY_GIC.sum(dim=["x", "y"]) / 1e6

dsRunoff_RACMO_1k_YY_GrIS_sections = (
    dsRunoff_RACMO_1k_YY_GrIS.groupby(dsmask_sections["section_numbers_adjusted"]).sum() / 1e6
)
dsRunoff_RACMO_1k_YY_GIC_sections = (
    dsRunoff_RACMO_1k_YY_GIC.groupby(dsmask_sections["section_numbers_adjusted"]).sum() / 1e6
)


# %%
folder500m = pathIMAU02 + "RACMO2.3p2/FGRN055/Downscaling_GR_500m/"
openMAR = True
dsRunoff500mRACMO = xr.open_mfdataset(
    folder500m + "Annual/runoff.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.YY_copy.nc"
)
# dsRunoff500mMAR   = xr.open_dataset(folder500m + "Annual/runoff.1940-2023.MAR3v14.GrIS.0.5km.YY.nc")

# SUM
if os.path.exists(
    folder500m
    + "Annual/sum/runoff_yearly_sum.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.YY.nc"
):
    dsRunoff500mRACMO_sum = xr.open_mfdataset(
        folder500m
        + "Annual/sum/runoff_yearly_sum.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.YY.nc"
    )
else:
    dsRunoff500mRACMO_sum = dsRunoff500mRACMO.runoffcorr.sum(dim=["x", "y"]) / 4
if openMAR:
    if os.path.exists(
        folder500m + "Annual/sum/runoff_yearly_sum.1940-2023.MAR3v14.GrIS.0.5km.YY.nc"
    ):
        dsRunoff500mMAR_sum = xr.open_dataset(
            folder500m + "Annual/sum/runoff_yearly_sum.1940-2023.MAR3v14.GrIS.0.5km.YY.nc"
        )
    else:
        dsRunoff500mMAR_sum = ds500mMAR.runoffcorr.sum(dim=["x", "y"]) / 4
    dsRunoff500mMAR_sum = dsRunoff500mMAR_sum / 1e6
    dsRunoff500mMAR_sum.attrs["units"] = "km3 w.e."
    dsRunoff500mMAR_sum.resample(time="YS").sum()

dsRunoff500mRACMO_sum = dsRunoff500mRACMO_sum / 1e6  # mm/km2 to km3
dsRunoff500mRACMO_sum.attrs["units"] = "km3 w.e."
dsRunoff500mRACMO_sum.resample(time="YS").sum()
# %%
path_sums_masks_500m = pathIMAU02 + "RACMO2.3p2/FGRN055/Downscaling_GR_500m/Annual/Sums and masks/"

file_paths = {
    "dsRunoff500mRACMO_GIC": path_sums_masks_500m
    + "runoff.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.GIC.YY.nc",
    "dsRunoff500mRACMO_GrIS": path_sums_masks_500m
    + "runoff.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.GrIS.YY.nc",
    "dsRunoff500mRACMO_GIC_sum": path_sums_masks_500m
    + "runoff_yearly_sum.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.GIC.YY.nc",
    "dsRunoff500mRACMO_GrIS_sum": path_sums_masks_500m
    + "runoff_yearly_sum.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.GrIS.YY.nc",
}
datasets_RACMO_runoff_500m = {}

for name, file_path in file_paths.items():
    if os.path.exists(file_path):
        datasets_RACMO_runoff_500m[name] = xr.open_dataset(file_path)
        print(f"Opened {file_path}")

dsRunoff500mRACMO_GIC_sum = datasets_RACMO_runoff_500m["dsRunoff500mRACMO_GIC_sum"] / 1e6
dsRunoff500mRACMO_GrIS_sum = datasets_RACMO_runoff_500m["dsRunoff500mRACMO_GrIS_sum"] / 1e6
dsRunoff500mRACMO_GIC = datasets_RACMO_runoff_500m["dsRunoff500mRACMO_GIC"]
dsRunoff500mRACMO_GrIS = datasets_RACMO_runoff_500m["dsRunoff500mRACMO_GrIS"]


fpaths_MAR = {
    "dsRunoff500mMAR_GIC": path_sums_masks_500m + "runoff.1940-2023.MAR3v14.GrIS.0.5km.GIC.YY.nc",
    "dsRunoff500mMAR_GrIS": path_sums_masks_500m + "runoff.1940-2023.MAR3v14.GrIS.0.5km.GrIS.YY.nc",
    "dsRunoff500mMAR_GIC_sum": path_sums_masks_500m
    + "runoff_yearly_sum.1940-2023.MAR3v14.GrIS.0.5km.GIC.YY.nc",
    "dsRunoff500mMAR_GrIS_sum": path_sums_masks_500m
    + "runoff_yearly_sum.1940-2023.MAR3v14.GrIS.0.5km.GrIS.YY.nc",
}

datasets_MAR_runoff_500m = {}

for name, file_path in fpaths_MAR.items():
    if os.path.exists(file_path):
        datasets_MAR_runoff_500m[name] = xr.open_dataset(file_path)
        print(f"Opened {file_path}")
dsRunoff500mMAR_GIC_sum = datasets_MAR_runoff_500m["dsRunoff500mMAR_GIC_sum"] / 1e6
dsRunoff500mMAR_GrIS_sum = datasets_MAR_runoff_500m["dsRunoff500mMAR_GrIS_sum"] / 1e6
dsRunoff500mMAR_GIC = datasets_MAR_runoff_500m["dsRunoff500mMAR_GIC"]
dsRunoff500mMAR_GrIS = datasets_MAR_runoff_500m["dsRunoff500mMAR_GrIS"]

for dataset in [
    dsRunoff500mRACMO_GIC_sum,
    dsRunoff500mRACMO_GrIS_sum,
    dsRunoff500mMAR_GIC_sum,
    dsRunoff500mMAR_GrIS_sum,
]:
    dataset.attrs["units"] = "km3 w.e."
    # for variable in dataset
    for variable in dataset.variables:
        # if 'units' in dataset[variable].attrs:
        # if dtype is float
        if "units" in dataset[variable].attrs:
            dataset[variable].values = dataset[variable].values / 1e6
            dataset = dataset.resample(time="YS").sum()

    # dataset.values = dataset.values/1e6
# %%
# data/temp/sections_500m_compressed.nc
mask_0_5km = xr.open_dataset("../../data/temp/sections_500m_compressed.nc")
mask_0_5_background = xr.open_dataset(
    pathIMAU02 + "RACMO2.3p2/FGRN055/Downscaling_GR_500m/GrIS_topo_icemask_lsm_lon_lat_0.5km.nc"
)
masks1k = xr.open_dataset("/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/masks1k.nc")

# %% SOLID DISCHARGE
dfGISDMankoff = pd.read_csv(path_Mankoff2020Solid + "GIS_D.csv", index_col=0, parse_dates=True)
dfErrorGISDMankoff = pd.read_csv(
    path_Mankoff2020Solid + "GIS_err.csv", index_col=0, parse_dates=True
)
dfGIScoverageMankoff = pd.read_csv(
    path_Mankoff2020Solid + "GIS_coverage.csv", index_col=0, parse_dates=True
)
filterCov = (dfGIScoverageMankoff < 0.5).values


if time_resolution == "Monthly":
    dfGISDMankoff = dfGISDMankoff / 12
    dfErrorGISDMankoff = dfErrorGISDMankoff / 12
elif time_resolution == "Annual":
    dfGISDMankoff = dfGISDMankoff.groupby(dfGISDMankoff.index.year).mean()
    dfErrorGISDMankoff = dfErrorGISDMankoff.groupby(dfErrorGISDMankoff.index.year).mean()


# %%  PRECIPITATION RACMO
file_annual_fjord = "../../data/temp/RACMO2.3p2_1km_precip_fjords_Annual_1958_2023.nc"


dsPrecipFjordsVol = xr.open_dataset(file_annual_fjord)
dsPrecipFjordsVolSum = dsPrecipFjordsVol.sum(dim=["x", "y"]).resample(time="YS").sum()

# %% PRECIPITATION CARRA
ds_precip_carra_1991_2008_sum = xr.open_dataset(
    pathIMAU02
    + "CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2008.1km.YY.fjords_only.sum_per_basin.nc"
)
ds_precip_carra_2009_2023_sum = xr.open_dataset(
    pathIMAU02
    + "CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.2009-2023.1km.YY.fjords_only.sum_per_basin.nc"
)
# concatenate
ds_precip_carra_1991_2023_sum = xr.concat(
    [ds_precip_carra_1991_2008_sum, ds_precip_carra_2009_2023_sum], dim="time"
)

dsPrecipFjordsCARRA_Annual_Sum = (
    (ds_precip_carra_1991_2023_sum.sum(dim="section_numbers_adjusted") / 1e6)
    .resample(time="YS")
    .sum()
)

# %%
for ds in [
    dsRunoff500mRACMO_GrIS_sum,
    dsRunoff500mRACMO_GIC_sum,
    dsRunoff500mMAR_GrIS_sum,
    dsRunoff500mMAR_GIC_sum,
]:
    ds = ds.resample(time="YS").sum()
# %%
dfRunoffTundra = (
    dsRunoffTundraSum.squeeze().to_dataframe(name="Runoff Tundra").drop(columns="height")
)
# dfRunoffIce = dsRunoffIce.sum(dim=['rlon', 'rlat']).squeeze().to_dataframe(name='Runoff Ice').drop(columns='height')
dfRunoffRACMOGrIS = (
    dsRunoff500mRACMO_GrIS_sum["runoffcorr"].squeeze().to_dataframe(name="Runoff GrIS")
)
dfRunoffRACMOGIC = (
    dsRunoff500mRACMO_GIC_sum["runoffcorr"].squeeze().to_dataframe(name="Runoff Ice Caps")
)
dfPrecipFjordsVol = dsPrecipFjordsVolSum["precipcorr"].to_dataframe(name="Precipitation Fjords")
dfSolidMankoff = dfGISDMankoff.copy(deep=True)

if time_resolution == "Annual":
    for df in [dfRunoffTundra, dfRunoffRACMOGrIS, dfRunoffRACMOGIC]:
        df.index = df.index.year
    dfPrecipFjordsVol = dfPrecipFjordsVol.groupby(dfPrecipFjordsVol.index.year).mean()

elif time_resolution == "Monthly":
    for df in [
        dfRunoffTundra,
        dfRunoffIce,
        dfPrecipFjordsVol,
    ]:
        df.index = df.index.strftime("%Y-%m")
    dfSolidMankoff.index = pd.to_datetime(dfSolidMankoff.index)
    dfSolidMankoff = dfSolidMankoff.groupby(dfSolidMankoff.index.strftime("%Y-%m")).mean()
#    dfPrecipFjordsVol = dfPrecipFjordsVol.groupby(dfPrecipFjordsVol.index.strftime('%Y-%m')).mean()

dfRunoff500mRACMO_GIC = (
    dsRunoff500mRACMO_GIC_sum["runoffcorr"].squeeze().to_dataframe(name="Runoff Ice Caps (RACMO 500m)")
)
dfRunoff500mMAR_GIC_sum = dsRunoff500mMAR_GIC_sum["runoffcorr"].squeeze().to_dataframe(name="Runoff Ice Caps (MAR 500m)")
dfRunoff_RACMO_1k_YY_GIC_sum = dsRunoff_RACMO_1k_YY_GIC_sum.squeeze().to_dataframe(name="Runoff Ice Caps (RACMO 1km)")
dfRunoff500mMAR_GrIS_sum = dsRunoff500mMAR_GrIS_sum["runoffcorr"].squeeze().to_dataframe(name="Runoff GrIS (MAR 500m)")
dfRunoff500mRACMO_GrIS_sum = dsRunoff500mRACMO_GrIS_sum["runoffcorr"].squeeze().to_dataframe(name="Runoff GrIS (RACMO 500m)")
dfRunoff_RACMO_1k_YY_GrIS_sum = dsRunoff_RACMO_1k_YY_GrIS_sum.squeeze().to_dataframe(name="Runoff GrIS (RACMO 1km)")
dfPrecipFjordsCARRA = dsPrecipFjordsCARRA_Annual_Sum["precip"].squeeze().to_dataframe(name="Precipitation Fjords (CARRA)")

dfPrecipFjordsCARRA.index = dfPrecipFjordsCARRA.index.year
dfRunoff500mMAR_GIC_sum.index = dfRunoff500mMAR_GIC_sum.index.year
dfRunoff500mRACMO_GIC.index = dfRunoff500mRACMO_GIC.index.year
dfRunoff_RACMO_1k_YY_GIC_sum.index = dfRunoff_RACMO_1k_YY_GIC_sum.index.year
dfRunoff500mMAR_GrIS_sum.index = dfRunoff500mMAR_GrIS_sum.index.year
dfRunoff500mRACMO_GrIS_sum.index = dfRunoff500mRACMO_GrIS_sum.index.year
dfRunoff_RACMO_1k_YY_GrIS_sum.index = dfRunoff_RACMO_1k_YY_GrIS_sum.index.year


df_sum_GIS_55 = pd.concat(
    [
        dfSolidMankoff,
        dfRunoff_RACMO_1k_YY_GIC_sum,
        dfRunoff_RACMO_1k_YY_GrIS_sum,
        dfRunoffTundra,
        dfPrecipFjordsVol,
        dfPrecipFjordsCARRA,
        dfRunoff500mMAR_GIC_sum,
        dfRunoff500mMAR_GrIS_sum,
        dfRunoff500mRACMO_GrIS_sum,
        dfRunoffRACMOGrIS,
        dfRunoffRACMOGIC,


    ],
    axis=1,
).sort_index()

if time_resolution == "Monthly":
    df_sum_GIS_55.index = pd.to_datetime(df_sum_GIS_55.index, format="%Y-%m")
else:
    df_sum_GIS_55.index = pd.to_datetime(df_sum_GIS_55.index, format="%Y")

filterNanMonths = (
    df_sum_GIS_55.isna().sum(axis=1) > 0
)  # find months with NaN values in any of the variables
df_sum_GIS_55Relative = (df_sum_GIS_55.iloc[:, :5].mask(filterNanMonths).T / df_sum_GIS_55.iloc[:, :5].sum(axis=1)).T


for ds in [
    dsRunoff500mRACMO_GrIS_sum,
    dsRunoff500mRACMO_GIC_sum,
    dsRunoff500mMAR_GrIS_sum,
    dsRunoff500mMAR_GIC_sum,
]:
    ds = ds.resample(time="YS").sum()


# %% =============
# PER SECTOR
# =============
def kgperm2_to_Gt(ds):
    ds = ds / 1e6
    ds.attrs["units"] = "Gt"
    return ds


def mask_MougBasins_ice(ds, IceOrTundra):
    """Mask the RACMO data with the Mouginot basins and sum the values for each basin
    IceOrTundra: "Ice" or "Tundra" to select the mask to use"""
    if spatial_resolution == "5_5k":

        mask55 = open_mask_5_5k(spatial_resolution)
        mask55["rlat"] = ds["rlat"]
        mask55["rlon"] = ds["rlon"]

        if IceOrTundra == "Ice":
            ds["Basins"] = mask55["Basins"]
            dsMougTime = ds.groupby("Basins").sum().isel(height=0)
        elif IceOrTundra == "Tundra":
            ds["Basins"] = mask55["Mouginot_Tundra"]
            dsMougTime = ds.groupby("Basins").sum().isel(height=0)
        ds.attrs["Description"] = (
            f"Sum per sector of the Mouginot basins for {spatial_resolution} RACMO2.3p2, for the {IceOrTundra} "
        )
        return dsMougTime

    elif spatial_resolution == "1k":
        if not "masks1k" in globals():
            mask1k = open_mask_1k()
        path_mask_1k_with_tundra = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/GrIS_topo_icemask_lsm_tundra_basins_lon_lat_1km.nc"
        mask_1k_with_tundra = xr.open_dataset(path_mask_1k_with_tundra)
        if IceOrTundra == "Ice":
            ds["Basins"] = mask_1k_with_tundra["Basins_All_Greenland"]
            dsMougTime = ds.groupby("Mouginot_basins").sum()
        elif IceOrTundra == "Tundra":
            ds["Mouginot_basins"] = mask_1k_with_tundra["Tundra_basins"]
            dsMougTime = ds.groupby("Mouginot_basins").sum()
        # give a attribute
        ds.attrs["Description"] = (
            f"Sum per sector of the Mouginot basins for {spatial_resolution} RACMO2.3p2, for the {IceOrTundra} "
        )
        return dsMougTime


# %% loading datasets


if "dsRunoff500mMAR_GIC_per_section" not in locals():
    print("Opening the data")
    dsRunoff500mMAR_GIC_per_section = xr.open_dataset(
        path_sums_masks_500m + "runoff.1940-2023.MAR3v14.GIC.0.5km.GIC.YY.per_section.nc"
    )
dsRunoff500mMAR_GIC_per_section

# add units
dsRunoff500mMAR_GIC_per_section.attrs["units"] = "km3 w.e. per year"
dsRunoff500mMAR_GIC_per_section.attrs["Description"] = (
    "Sum per sector of the Greenland Ice Caps for the 500m MAR3.14 data"
)
#  GrIS
dsRunoffIceSectormm = xr.open_dataset(
    pathIMAU02
    + "RACMO2.3p2/FGRN055/Downscaling_GR/Monthly/runoff_GrIS.1990-2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.GrIS.section_sum.nc"
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

if "dsRunoff500mMAR_GrIS_per_section" not in locals():
    print("Opening the data")
    dsRunoff500mMAR_GrIS_per_section = xr.open_dataset(
        path_sums_masks_500m + "runoff.1940-2023.MAR3v14.GrIS.0.5km.GrIS.YY.per_section.nc"
    )
# (dsRunoff500mMAR_GrIS_per_section/1e6).plot(hue = 'section_numbers_adjusted')
# add attributes km3 w.e. per year
dsRunoff500mMAR_GrIS_per_section.attrs["units"] = "km3 w.e. per year"
dsRunoff500mMAR_GrIS_per_section.attrs["Description"] = (
    "Sum per sector of the Greenland Ice Sheet for the 500m MAR3.14 data"
)

# Ice caps
dfRunoff500mMAR_GIC_per_section = (
    dsRunoff500mMAR_GIC_per_section["runoffcorr"].to_dataframe().reset_index()
)
dfRunoff500mMAR_GIC_per_section["section_numbers_adjusted"] = dfRunoff500mMAR_GIC_per_section[
    "section_numbers_adjusted"
].map(dict_sections)
dfRunoff500mMAR_GIC_per_section = (
    dfRunoff500mMAR_GIC_per_section.set_index("time").pivot(columns="section_numbers_adjusted")
    / 1e6
)
dfRunoff500mMAR_GIC_per_section.columns = dfRunoff500mMAR_GIC_per_section.columns.get_level_values(
    1
)

dsRunoffIceCapSectormm = xr.open_dataset("/Volumes/imau02/rapid/Anneke/RACMO2.3p2/FGRN055/Downscaling_GR/Monthly/runoff_GIC.1990-2023.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.GIC.section_sum.nc")['runoffcorr']
dsRunoffIceCapSector = kgperm2_to_Gt(dsRunoffIceCapSectormm)
dfRunoffIceCapSector = dsRunoffIceCapSector.to_dataframe(name='Liquid Runoff Ice Sheet').reset_index().set_index('time').pivot(columns='section_numbers_adjusted').rename(columns=dict_sections)
dfRunoffIceCapSector.columns = dfRunoffIceCapSector.columns.get_level_values(1)


# Tundra
fname = (
    pathIMAU02
    + "RACMO2.3p2/FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_may24_copy.nc"
)
mask = xr.open_dataset(fname)
fname = pathIMAU02 + "RACMO2.3p2/FGRN055/Downscaling_GR/masks1k.nc"
mask1 = xr.open_dataset(fname)
mask_tundra = mask1["LSM"].where(mask1["Promicemask"] == 0)
mask["Tundra_sections"] = mask["section_numbers_adjusted"].where(mask_tundra.values == 1)
mask_tundra_sections = mask["Tundra_sections"]
mask_tundra_sections = mask_tundra_sections.rename({"x": "rlon", "y": "rlat"})


time_resolution = "Monthly"
years = list(range(2012, 2023))

dsRunoffTundra_mmyear = read_RACMO("1k_reproject", time_resolution, years=years, variable="runoff")

# rename x and y to rlon and rlat
dsRunoffTundra_mmyear["Tundra_basins"] = mask_tundra_sections
path_monthly_racmo = pathAnnekeFolderIMAU02 + "/Downscaling_GR/Monthly/"
fname_runofftundraMM = (
    f"runoff_tundra.{years[0]}-{years[1]}.RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.tundra.section_sum.nc"
)

if os.path.isfile(path_monthly_racmo + fname_runofftundraMM):
    dsRunoffTundra_mmyear_SECTOR = xr.open_dataset(path_monthly_racmo + fname_runofftundraMM)
else:
    dsRunoffTundra_mmyear_SECTOR = (
        dsRunoffTundra_mmyear["runoff"]
        .groupby(dsRunoffTundra_mmyear["Tundra_basins"])
        .sum()
        .drop("height")
    )
    dsRunoffTundra_mmyear_SECTOR.to_netcdf(path_monthly_racmo + fname_runofftundraMM)


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
# dfRunoffTundraSector['2010':'2012'].plot(kind='area',stacked=True, color=colors_hex_Moug, figsize=(15, 4), title = "Liquid Runoff Tundra")
#
#

dfSectionD = pd.read_csv(
    path_Mankoff2020Solid_adjusted + "section_D.csv", index_col=0, parse_dates=True
)
dfErrorSectionDMankoff = pd.read_csv(
    path_Mankoff2020Solid_adjusted + "section_err.csv", index_col=0, parse_dates=True
)
dfSectionCovManoff = pd.read_csv(
    path_Mankoff2020Solid_adjusted + "section_coverage.csv", index_col=0, parse_dates=True
)
# filterCovReg =  (dfSectionCovManoff < 0.5).values

dfSectionD = dfSectionD / 12
dfErrorSectionDMankoff = dfErrorSectionDMankoff / 12
dfSectionD.columns.name = "Basins"


# open ds_precip_carra_1991_2023_sum.to_netcdf(pathIMAU02 +"CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2023.1km.YY.fjords_only.sum_per_basin.nc")
ds_precip_carra_1991_2023_sum = (
    xr.open_mfdataset(
        [
            pathIMAU02
            + "CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2008.1km.YY.fjords_only.sum_per_basin.nc",
            pathIMAU02
            + "CARRA/Yearly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.2009-2023.1km.YY.fjords_only.sum_per_basin.nc",
        ]
    )
    / 1e6
)
ds_precip_carra_1991_2023_sum.resample(time="YS").sum()


dsPrecipFjords = xr.open_mfdataset(
    pathDataTemp + "RACMO2.3p2_1km_precip_fjords_Annual_1958_2023.nc"
)
dsPrecipFjordsAnnualRACMO_58_23 = dsPrecipFjords.copy(deep=True)

dsmask_sections = xr.open_mfdataset(mask_sections)
dsPrecipFjordsSectormm_sum = dsPrecipFjords.groupby(
    dsmask_sections["section_numbers_adjusted"]
).sum()
dsPrecipFjordsSectormm_sum = dsPrecipFjordsSectormm_sum.resample(time="YS").sum()
# %%
monthly_precip_cara1 = (
    pathIMAU02
    + "CARRA/Monthly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.2009-2023.1km.MM.fjords_only.sum_per_basin.nc"
)
ds_precip_carra_2009_2023_month = xr.open_dataset(monthly_precip_cara1) / 1e6

ds_precip_carra_1991_2008_month = (
    xr.open_dataset(
        pathIMAU02
        + "CARRA/Monthly/RACMOgrid/fjords_only/total_precipitation.CARRA.west_domain.1991-2008.1km.MM.fjords_only.sum_per_basin.nc"
    )
    / 1e6
)

ds_precip_carra_1991_2023_month = xr.concat(
    [ds_precip_carra_1991_2008_month, ds_precip_carra_2009_2023_month], dim="time"
)

ds_precip_racmo_1990_2023_sum_monthly = xr.open_mfdataset(
    pathIMAU02
    + "RACMO2.3p2/FGRN055/Downscaling_GR/Monthly/fjords_only/precip.1990-2023.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.fjords_only.sum_per_basin.nc",
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
dfPrecipFjordsSector_RACMO_year.columns = dfPrecipFjordsSector_RACMO_year.columns.get_level_values(
    1
)

dfPrecipFjordsSector_CARRA_year = (
    ds_precip_carra_1991_2023_sum.to_dataframe()
    .reset_index()
    .set_index("time")
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)
dfPrecipFjordsSector_CARRA_year.columns = dfPrecipFjordsSector_CARRA_year.columns.get_level_values(
    1
)


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
dfPrecipFjordsSector_RACMO_month.columns = (
    dfPrecipFjordsSector_RACMO_month.columns.get_level_values(1)
)

dfPrecipFjordsSector_CARRA_month = (
    ds_precip_carra_1991_2023_month.to_dataframe()
    .reset_index()
    .set_index("time")
    .pivot(columns="section_numbers_adjusted")
    .rename(columns=dict_sections)
)
dfPrecipFjordsSector_CARRA_month.columns = (
    dfPrecipFjordsSector_CARRA_month.columns.get_level_values(1)
)
# %%
open_previous = True

if open_previous:
    dsSectorSum = xr.open_dataset(pathDataTemp + "RACMO2.3p2_1k_sector_sum_2024_06_12.nc")

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
            .rename({"section_numbers_adjusted": "Basins"}),
            dfPrecipFjordsSector_CARRA_month.stack()
            .to_xarray()
            .to_dataset(name="Precipitation Fjords CARRA")
            .rename({"section_numbers_adjusted": "Basins"}),
        ]
    )
    # dsSectorSumD['Solid Ice discharge original'] = dfRegionDMankoff.stack().to_xarray()
    # dsSectorSum['Solid Ice discharge (interp)'] = dsSectorSum['Solid Ice discharge original'].interp(Date=dsSectorSum.time, method='linear')
    # # groupby month and mean
    # dsSectorSum['Solid Ice discharge (mean)'] = dsSectorSum['Solid Ice discharge original'].groupby('Date.month').mean()

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

    def compute_monthly_mean(data):
        return data.resample(Date="1ME").mean()

    dsSectorSum["Solid Ice discharge (weighted mean)"] = (
        compute_weighted_monthly_mean(dfSectionD.stack().to_xarray())
        .interp(Date=dsSectorSum.time, method="nearest")
        .drop("Date")
    )

    #
    # export dsSectorSum to csv
    start = dsSectorSum.time[0].dt.strftime("%Y").values
    end = dsSectorSum.time[-1].dt.strftime("%Y").values
    dsSectorSum.to_netcdf(pathDataTemp + f"RACMO2.3p2_1k_sector_sum_2024_06_13.nc")

dsSectorSum = dsSectorSum.where(dsSectorSum != 0)


# dsSectorSum.resample(time='MS').mean().to_dataframe().reset_index().set_index('time').to_csv(pathDataTemp + f"RACMO2.3p2_1k_sector_sum_2024_06_12.csv")

(dsSectorSum.resample(time="YS").mean() * 12).sum(
    dim="Basins"
).to_dataframe().reset_index().set_index("time").to_csv(
    pathDataTemp + f"Greenland_Sum_2024_06_12.csv"
)

# Calculate the mean of each DataArray in the Dataset
means = {var: dsSectorSum[var].mean().item() for var in dsSectorSum.data_vars}

# Convert the means to a DataFrame
variables_sorted = (
    pd.DataFrame(list(means.items()), columns=["Variable", "Mean"])
    .sort_values(by="Mean", ascending=False)["Variable"]
    .values
)

dsSectorSum = dsSectorSum[variables_sorted]

# plot ds solid discharge
col_order_abs = [
    "Solid Ice Discharge",
    "Precipitation Fjords",
    "Liquid Runoff Tundra",
    "Liquid Runoff Ice Caps",
    "Liquid Runoff Ice Sheet",
]
col_order_rel = [
    "Solid Ice Discharge",
    "Liquid Runoff Ice Sheet",
    "Liquid Runoff Ice Caps",
    "Liquid Runoff Tundra",
    "Precipitation Fjords",
]


# make a list of all the names of data variables in dsSectorSum
period1 = {"start": "1990", "end": "2004"}
period2 = {"start": "2005", "end": "2022"}

try:
    dsMonthlyGr = dsSectorSum.copy(deep=True).drop("winter_year").drop("summer_year")
except:
    dsMonthlyGr = dsSectorSum.copy(deep=True)
# make all 0 values nan
dsMonthlyGr = dsMonthlyGr.where(dsMonthlyGr > 0)
# rename to short names
dsMonthlyGr = dsMonthlyGr.rename(dict_shorter_name)
# only select col_order_rel
dsMonthlyGr = dsMonthlyGr[col_order_rel]
dsMonthlyGr = dsMonthlyGr.resample(time="MS").mean()

# Calculate the mean and standard deviation for each month
data_mean_seasonal_period2 = (
    dsMonthlyGr.sel(time=slice(period2["start"], period2["end"]))
    .sum(dim="Basins")
    .groupby("time.month")
    .mean()
)
data_std_seasonal_period2 = (
    dsMonthlyGr.sel(time=slice(period2["start"], period2["end"]))
    .sum(dim="Basins")
    .groupby("time.month")
    .std()
)
# do the same for 2002-2012
data_mean_seasonal_period1 = (
    dsMonthlyGr.sel(time=slice(period1["start"], period1["end"]))
    .sum(dim="Basins")
    .groupby("time.month")
    .mean()
)
data_std_seasonal_period1 = (
    dsMonthlyGr.sel(time=slice(period1["start"], period1["end"]))
    .sum(dim="Basins")
    .groupby("time.month")
    .std()
)

dfGIS_period1 = pd.read_csv(path_Mankoff2020Solid + "GIS_D.csv", index_col=0, parse_dates=True)[
    period1["start"] : period1["end"]
]
dfGIS_period2 = pd.read_csv(path_Mankoff2020Solid + "GIS_D.csv", index_col=0, parse_dates=True)[
    period2["start"] : period2["end"]
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
# %% MAR 1 KM RUNOFF

import gzip
import shutil

def convert_months_to_date(data_array, start_date): 
    '''
    This function converts the months since a start date to a date. For format dat_array as ds.time, start_date should 'YYYY-MM-DD'
    '''
    array_dates = np.empty(len(data_array)).astype('datetime64[ns]')
    for i in range(len(data_array)):
        array_dates[i] = pd.to_datetime(start_date) + pd.DateOffset(months=np.floor(data_array[i]), days=((data_array[i]%1)*30))
    return array_dates


def convert_years_to_date(data_array, start_date): 
    '''
    This function converts the years since a start date to a date. For format dat_array as ds.time, start_date should 'YYYY-MM-DD'
    '''
    array_dates = np.empty(len(data_array)).astype('datetime64[ns]')
    for i in range(len(data_array)):
        array_dates[i] = pd.to_datetime(start_date) + pd.DateOffset(years=np.floor(data_array[i]), days=((data_array[i]%1)*365))
    return array_dates


def open_compressed_xarray(file_path):
    '''
    Opens a compressed netcdf file with xarray
    '''
    # Path to the compressed file
    compressed_file = file_path

    # Path to the decompressed file
    decompressed_file = compressed_file.replace('.gz', '')

    # Decompress the file
    with gzip.open(compressed_file, 'rb') as f_in:
        with open(decompressed_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Open the decompressed file with xarray
    ds = xr.open_dataset(decompressed_file, engine='netcdf4', decode_times=False)
    return ds

fpath_adj_sect  = '/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/FWclean/data/temp/adjusted_section_numbers_slater.nc'
fpath_masks1k = '/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/FWclean/data/temp/masks1k.nc'
folder_MARRACMO1km = "/Users/annek/Library/CloudStorage/OneDrive-SharedLibraries-NIOZ/PhD Anneke Vries - General/FWclean/data/raw/liquid/"

# open all files and align time
ds_masks1k = xr.open_dataset(fpath_masks1k)
ds_adj_sect = xr.open_dataset(fpath_adj_sect)

ds_run_MAR = open_compressed_xarray(folder_MARRACMO1km + "runoff.1940-2023.MARv3.14-ERA5.1km.YY.nc.gz")
ds_run_MAR['years_since_19400115'] = ds_run_MAR.time
ds_run_MAR['time'] = convert_months_to_date(ds_run_MAR['years_since_19400115'], '1940-01-15')

# get masks in same xarray coordinate system
ds_run_MAR_mean = ds_run_MAR.mean(dim=['time'])
ds_run_MAR_mean['section_numbers_adjusted'] = ds_run_MAR_mean['runoffcorr'].copy(deep=True)
ds_run_MAR_mean['section_numbers_adjusted'].values = ds_adj_sect['section_numbers_adjusted'].values
for var in ds_masks1k.data_vars:
    ds_run_MAR_mean[var] = ds_run_MAR_mean['runoffcorr'].copy(deep=True)
    ds_run_MAR_mean[var].values = ds_masks1k[var].values


# calculate runoff for each basin yearly in Gt
ds_run_MAR_GrIS_basin = (ds_run_MAR['runoffcorr'].where(ds_run_MAR_mean['GrIS']==1).groupby(ds_run_MAR_mean['section_numbers_adjusted']).sum()/1e6)
ds_run_MAR_GIC_basin = (ds_run_MAR['runoffcorr'].where(ds_run_MAR_mean['GIC']==1).groupby(ds_run_MAR_mean['section_numbers_adjusted']).sum()/1e6)


ds_run_RACMO = open_compressed_xarray(folder_MARRACMO1km+ "runoff.1958-2023.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.YY.nc.gz")
ds_run_RACMO['years_since_19580115'] = ds_run_RACMO.time
ds_run_RACMO['time'] = convert_years_to_date(ds_run_RACMO['years_since_19580115'], '1958-01-15')
ds_runoff_RACMO_mean = ds_run_RACMO.sum(dim=['time'])
ds_runoff_RACMO_mean['section_numbers_adjusted'] = ds_runoff_RACMO_mean['runoffcorr'].copy(deep=True)
ds_runoff_RACMO_mean['section_numbers_adjusted'].values = ds_adj_sect['section_numbers_adjusted'].values

for var in ds_masks1k.data_vars:
    ds_runoff_RACMO_mean[var] = ds_runoff_RACMO_mean['runoffcorr'].copy(deep=True)
    ds_runoff_RACMO_mean[var].values = ds_masks1k[var].values


ds_run_RACMO_GrIS_basin = (ds_run_RACMO['runoffcorr'].where(ds_runoff_RACMO_mean['GrIS']==1).groupby(ds_runoff_RACMO_mean['section_numbers_adjusted']).sum()/1e6)
ds_run_RACMO_GIC_basin = (ds_run_RACMO['runoffcorr'].where(ds_runoff_RACMO_mean['GIC']==1).groupby(ds_runoff_RACMO_mean['section_numbers_adjusted']).sum()/1e6)
