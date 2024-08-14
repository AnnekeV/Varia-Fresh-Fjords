import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

try:
    from scripts.find_largest_islands import *
    from scripts.paths import *
except:
    from find_largest_islands import *
    from paths import *

print(f"Pathdata: {pathdata}")


pathlocal = "/Users/annek/Documents/RACMO2.3p2/FGRN055/"

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
FW_type = ["Solid", "IceRun", "Tundra", "Precip"]
dictFWcolor = dict(zip(FW_type, colors))


def pick_file_from_list(path, string1,string2=None):
    """Returns a list of files in the specified path that contain the specified string"""
    # Define the strings to search for

    # Find files that contain both strings
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if string2 != None:
                if string1 in filename and string2 in filename:
                    matches.append(os.path.join(root, filename))
            else:
                if string1 in filename:
                    matches.append(os.path.join(root, filename))
    files_containing_specific_string = matches

    if len(files_containing_specific_string) == 0:
        print(f"No file found for {string1}, {string2} in {path}. \nPlease specify the file you want to open.")
        return
    print(f"Opening file: \n {files_containing_specific_string[0]}\n... \n{files_containing_specific_string[-1]}")
    
    return files_containing_specific_string


def read_RACMO(downscaling_type, time_resolution, years, variable, since_1940 =True):
    """Reads RACMO data from the specified path and returns a xarray dataset
    downscaling type = "1k" or "5_5k" or "1k_reproject" 
    1k_reproject is the 5_5k data reprojected to 1k and should only be used for runoff on tundra
    time_resolution = "Annual", "Monthly", "Daily" or "Yearly"
    If time_resolution is "Yearly", you don't need to specify the years
    years = list of years, e.g. ["2010", "2020"]
    variable = e.g. "runoff", "precipitation", "temperature"
    returns a xarray dataset 
    """
    if variable == "snowfrac":
        path = "/Users/annek/Documents/RACMO2.3p2/FGRN055"
    else:
        path = pathdata

    if downscaling_type == "1k":
        if time_resolution == "Annual" or time_resolution == "Yearly":
            path1k= os.path.join(path,"Downscaling_GR", "Annual", )
            fpath = pick_file_from_list(path1k, variable)
            if len(fpath) > 1:
                print("More than one file found. Please specify the file you want to open.")
                return
            elif len(fpath) == 0:
                print(f"No file found for {variable}. Please specify the file you want to open.")
                return
            if variable == "precipitation":
                fpath=  "/Volumes/imau01/rapid/RACMO2.3p2/FGRN055/Downscaling_GR/Annual/precip.1958-2023.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.YY.nc.gz"
                precipitation_annual = "/Volumes/imau01/rapid/RACMO2.3p2/FGRN055/Downscaling_GR/Annual/precip.1958-2023.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.YY.nc.gz"
                ds_precip = xr.open_mfdataset(precipitation_annual, decode_times=False)
                ds_precip['time'] = pd.to_datetime(decadal_year, format='%Y') +pd.to_timedelta(decadal_year%1*365, unit='d')
                return ds_precip
            try: 
                ds = xr.open_mfdataset(fpath[0])
                ds = ds.sel(time=slice(str(years[0]), str(years[-1])))

            except:
                ds = xr.open_mfdataset(fpath[0], decode_times=False)
                ds['time'] = pd.to_datetime(ds.time.values + 1958, format='%Y')


        elif time_resolution == "Monthly":
            print(" it might be hard to process the data for the whole period, please be patient, or choose a lower time resolution or spatial resolution")
            list_ds_paths = []
            ds_list = []
            print(path)
            for year in tqdm(np.arange(int(years[0]), int(years[-1])+1)):
                path1k= os.path.join(path,   "Downscaling_GR", time_resolution, str(year))
                fpath = pick_file_from_list(path1k, variable)
                list_ds_paths.append(fpath[0])

                ds_single = xr.open_mfdataset(fpath[0])
                if years[0] == years[-1]:
                    ds = ds_single
                    return ds
                else:
                    ds_list.append(ds_single)
            # print("Open mfdataset")
            # ds = xr.open_mfdataset(list_ds_paths, combine='by_coords')
            if len(ds_list)>10: print("More than 10 files, this might take a while to concatenate")
            ds = xr.concat(ds_list, dim="time")
        elif time_resolution == "Daily":
            print(" it might be hard to process the data for the whole period, please be patient, or choose a lower time resolution or spatial resolution. Ideally do not open more than 1 year at a time")
            ds_list = []
            path1k= os.path.join(path,   "Downscaling_GR", time_resolution, variable)
            for year in tqdm(np.arange(int(years[0]), int(years[-1])+1)):
                files = os.listdir(path1k)
                # open all files that contain year in the name
                for file in files:
                    if year in file:
                        ds = xr.open_mfdataset(os.path.join(path_1k, file))
                        ds_list.append(ds)
                if years[0] == years[-1]:
                    ds = ds_single
                    return ds
                else:
                    ds_list.append(ds_single)
            ds = xr.concat(ds_list, dim="time")

        else:
            print("Time resolution not supported. Please choose 'Annual' or 'Monthly, or 'Daily'")

    elif downscaling_type == "1k_reproject":
        ''' use locally reprojected data for tundra runoff'''
        print("Reading locally reprojected 1k data")
        if time_resolution == 'Monthly':
            path1k= os.path.join(pathlocal,"Downscaling_GR/Monthly/runoff_monthlyS_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303_1km.nc")
            ds = xr.open_mfdataset(path1k)

        elif time_resolution == 'Annual'or time_resolution == 'Yearly':
            path1k= os.path.join(pathlocal,"Downscaling_GR/Annual/runoff_yearlyS_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_1940_2022_1km.nc")
            ds = xr.open_mfdataset(path1k)
        print("Finished reading locally reprojected 1k data")

        
    elif downscaling_type=="5_5k":
        if time_resolution == "Annual" or time_resolution == "Yearly":
            if since_1940: time_resolution = "Yearly"
            if time_resolution =="Yearly":
                path5_5k = os.path.join(path+"_1940", time_resolution)
                fpath = pick_file_from_list(path5_5k, variable)
                ds = xr.open_mfdataset(fpath[0])
            else:    
                path5_5k = os.path.join(path, "BN_RACMO2.3p2_ERA5_3h_FGRN055", time_resolution)
                if years[0] == years[-1]:
                    fpath = pick_file_from_list(path5_5k, variable, years[0])
                    ds = xr.open_mfdataset(fpath[0])
                else:
                    ds_list = []
                    for year in tqdm(np.arange(int(years[0]), int(years[-1])+1)):
                        fpath = pick_file_from_list(path5_5k, variable, str(year))
                        ds_single = xr.open_mfdataset(fpath[0])
                        ds_list.append(ds_single)
                    ds = xr.concat(ds_list, dim="time")
        

        elif time_resolution == "Monthly":
            if since_1940:
                path5_5k = os.path.join(path+"_1940", time_resolution)
                fpath = pick_file_from_list(path5_5k, variable)
                print(fpath)
                ds = xr.open_mfdataset(fpath[0])
            elif years[0] == years[-1]:  # not since 1940
                path5_5k = os.path.join(path, "BN_RACMO2.3p2_ERA5_3h_FGRN055", time_resolution, years[0])
                fpath = pick_file_from_list(path5_5k, variable, years[0])
                ds = xr.open_mfdataset(fpath[0])
            elif years[0] != years[-1]:
                ds_list = []
                for year in tqdm(np.arange(int(years[0]), int(years[-1])+1)):
                    path5_5k = os.path.join(path, "BN_RACMO2.3p2_ERA5_3h_FGRN055", time_resolution, str(year))
                    fpath = pick_file_from_list(path5_5k, variable, str(year))
                    ds_single = xr.open_mfdataset(fpath[0])
                    ds_list.append(ds_single)
                ds = xr.concat(ds_list, dim="time")
        print(f"Finished reading {fpath[0]}")
    return ds


# ============ 5_5k ============

def open_mask_5_5k():
    pathmasks5_5k = pathdata +"/BN_RACMO2.3p2_ERA5_3h_FGRN055/FGRN055_Masks.nc"
    masks5_5k = xr.open_mfdataset(pathmasks5_5k, engine='netcdf4')
    # Tundra Mouginot region mask
    path_mask_incl_Tundra = pathdata + "/BN_RACMO2.3p2_ERA5_3h_FGRN055/FGRN055_drainage_basins_all_land.nc"
    mask_incl_Tundra = xr.open_mfdataset(path_mask_incl_Tundra)
    masks5_5k['Mouginot_AllGreenland'] = masks5_5k['Mouginot_basins'].copy(deep=True)
    masks5_5k['Mouginot_AllGreenland'].values = mask_incl_Tundra['Basins'].values
    masks5_5k['Mouginot_Tundra'] = masks5_5k['Mouginot_AllGreenland'].copy(deep=True)
    masks5_5k['Mouginot_Tundra'].values[np.where(masks5_5k['Mouginot_basins'].values != 0)] = np.nan

    return masks5_5k

def mask_data_5_5k(ds, masktype):
    masks5_5k = open_mask_5_5k()
    masks5_5k['rlat'] = ds['rlat']
    masks5_5k['rlon'] = ds['rlon']
    
    if masktype == "PROMICE_GrIS":
        # 4 GrIS
        ds_new = ds.where(masks5_5k['Promicemask'] == 4)
    elif masktype == "PROMICE_Ice_Total":
        # 2,3 ice caps and 4 GrIS
        ds_new = ds.where(masks5_5k['Promicemask'] >=2 ) 
    elif masktype == "LSM All Greenland":
        ds_new = ds.where(masks5_5k['LSM_GR'] == 1)
    elif masktype == "LSM Tundra":
        # LSM mask without ice sheet from PROMICE
        ds_new = ds.where((masks5_5k['LSM_GR']== 1) & (masks5_5k['Promicemask'] == 0))
    
    else: 
        print("Mask not found")
        print("Choose between PROMICE_GrIS, PROMICE_Ice_Total, LSM All Greenland, LSM No Ice Sheet")
        return
    return ds_new

# ============ 1k ============

def open_fjord_mask():
    dsFjords = xr.open_mfdataset(pathFjordMask, engine='netcdf4')
    return dsFjords["fjords"]
def make_LSM_mask_1k(masks1k, elevation_threshold=1):
  ''' 
  Function to create a land sea mask based on topography, standard elevation threshold is 1 m, also exclude weird values in top right corner
  '''
  topographyBinary = xr.where(((masks1k['Topography']>=elevation_threshold)
                              & ~((masks1k.x > 6000)&(masks1k.y>12000))),
                                1, 0)
  masks1k['LSM'] =topographyBinary
  return masks1k

def openLSMmask1k(ds):
    '''Open 1k land sea mask as send by Brice to Anneke on March 7 2024'''
    dsLSMmask1k = xr.open_mfdataset(pathlocal +"Downscaling_GR/GrIS_topo_icemask_lsm_tundra_basins_ocean_lon_lat_1km.nc")
    ds['LSM'] = ds['Topography'].copy(deep=True)
    ds['LSM'].values = dsLSMmask1k['LSM'].values
    return ds

def openOtherMasks1k(ds):
    '''Open 1k masks made by Anneke on March 7 2024, such as basins for Ocean and  Tundra on 1k Tundra is made by reprojecting the 5_5k mask to 1k with nearest neighbour interpolation.'''  
    dsOtherMasks1k = xr.open_mfdataset(pathlocal + "Downscaling_GR/GrIS_topo_icemask_lsm_tundra_basins_ocean_lon_lat_1km.nc")
    # for all variables that are not yet in the dataset add them
    for var in dsOtherMasks1k:
        if var not in ds:
            ds[var] = dsOtherMasks1k[var]
    return ds

def make_mask_largest_islands(masks1k):
    '''Imports a mask of the biggest and second biggest island in the LSM mask, and adds it to the dataset as LSM_no_islands and LSM_Qeqertarsuaq.
    '''
    # if not file path_mask_islands
    if os.path.isfile(path_mask_islands):
        # open path_mask_islands
        mask_islands = xr.open_mfdataset(path_mask_islands)

        masks1k['LSM_no_islands'] = masks1k['LSM'].copy(deep=True)
        masks1k['LSM_no_islands'].values = mask_islands['LSM_no_islands'].values

        masks1k['LSM_Qeqertarsuaq'] = masks1k['LSM'].copy(deep=True)
        masks1k['LSM_Qeqertarsuaq'].values = mask_islands['LSM_Qeqertarsuaq'].values
        return masks1k
    else:
        print("Mask for largest islands not found, please run find_largest_islands.py first")

def open_mask_1k():
    pathmasks1k =  os.path.join(pathAnnekeFolderIMAU02, "Downscaling_GR",  "Icemask_Topo_Iceclasses_lon_lat_average_1km.nc")
    masks1k = xr.open_mfdataset(pathmasks1k, engine='netcdf4')
    masks1k = openLSMmask1k(masks1k)
    masks1k = make_mask_largest_islands(masks1k)
    masks1k['fjords'] = open_fjord_mask()
    return masks1k

def mask_data_1k(ds, masktype):
    masks1k = open_mask_1k()
    
    masks1k['x'] = ds['x']
    masks1k['y'] = ds['y']
    
    if masktype == "PROMICE_GrIS":
        # 3 GrIS
        ds_new = ds.where(masks1k['Promicemask'] == 3)
    elif masktype == "PROMICE_Ice_Total":
        # 2,1 ice caps and 3GrIS
        ds_new = ds.where(masks1k['Promicemask'] >=1 ) 
    elif masktype == "PROMICE_Ice_caps":
        # 1,2 ice caps where is 1 or 2 
        ds_new = ds.where((masks1k['Promicemask'] == 1)| (masks1k['Promicemask'] == 2))
    elif masktype == "LSM All Greenland":
        ds_new = ds.where(masks1k['LSM'] == 1)
    elif masktype == "LSM Tundra":
        # LSM mask without ice sheet from PROMICE
        ds_new = ds.where((masks1k['LSM']== 1) & (masks1k['Promicemask'] == 0))
    elif masktype == "Fjords":
        ds_new = ds.where(masks1k['fjords'] == 1)
    else: 
        print("Mask not found")
        print("Choose between PROMICE_GrIS, PROMICE_Ice_Total, LSM All Greenland, LSM No Ice Sheet")
        return
    return ds_new

def mask_data(ds, masktype, spatial_resolution):
    if spatial_resolution == "1k":
        return mask_data_1k(ds, masktype)
    elif spatial_resolution == "5_5k":
        return mask_data_5_5k(ds, masktype)
    

# ========================================================================
#  Area and volume calculations
# =========================================================================
def area(ds, scale):
    if scale == "5_5k":
        masks5_5k = open_mask_5_5k()
        masks5_5k['rlat'] = ds['rlat']
        masks5_5k['rlon'] = ds['rlon']
        surface_area_per_gridcell = masks5_5k['Area']*1e6
    elif scale == "1k":
        surface_area_per_gridcell = 1e6 # in m2
    elif scale == "500m":
        surface_area_per_gridcell = 500*500.
    else:
        print("Scale not found, choose between 5_5k and 1k")
        return
    return surface_area_per_gridcell

def volume(ds, scale):
        '''Calculate volume of runoff in km3, given the dataset and the scale of the mask, first computes area per gridcell in m2, then multiplies with the runoff in mm w.e. and converts to km3.'''
        # add attribute
        ds.attrs['units'] = "km3 w.e."
        ds.attrs['long_name'] = "Runoff volume calculated from mm w.e. runoff and area per gridcell"
        ds.values = ds.values* area(ds, scale)/1000./1.e9
        return ds

if __name__ == '__main__':
    """Executed from the command line"""
    read_RACMO("1k", "Monthly",years=["2010", "2020"], variable="runoff")


else:
    """Executed on import"""


    # path = "/Users/annek/Documents/RACMO2.3p2/FGRN055"
    pass
