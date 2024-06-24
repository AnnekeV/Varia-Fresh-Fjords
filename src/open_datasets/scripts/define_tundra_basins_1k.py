import xarray as xr


#======= TUNDRA ===========
path_mask_1k_with_LSM = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/GrIS_topo_icemask_lsm_lon_lat_1km.nc"
mask_1k_with_LSM = xr.open_dataset(path_mask_1k_with_LSM)
mask_1k_with_LSM['Tundra'] = mask_1k_with_LSM['LSM']- mask_1k_with_LSM['Icemask']

ds_Basins_1k = xr.open_dataset("/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/FGRN055_drainage_basins_all_land_1km.nc")
mask_1k_with_LSM['Tundra_basins'] = ds_Basins_1k['Basins'].values * mask_1k_with_LSM['Tundra']
mask_1k_with_LSM['Basins_All_Greenland'] = ds_Basins_1k['Basins'].values * mask_1k_with_LSM['LSM']
mask_1k_with_LSM['Basins_All_Greenland'] = mask_1k_with_LSM['Basins_All_Greenland'].where(mask_1k_with_LSM['LSM'] == 1)



#======= OCEAN BASINS ===========
# also in RACMO_local_example

y = mask_1k_with_LSM.y
x = mask_1k_with_LSM.x

border_CE_SE = (y>-2.55e6) & (x>0)
border_CE_NE = (y<-1.85e6)
border_CE = border_CE_SE & border_CE_NE
# SE
border_SE_SW = (x>5e4) 
border_SE_CE = (y<-2.55e6)
border_SE = border_SE_SW & border_SE_CE
# SW 
border_SW_CW = (y<-2.3e6) 
border_SW_SE =  (x<5e4)
border_SW = border_SW_CW & border_SW_SE
# CW
border_CW_NW = (y<-2e6) & (x<-2e5)
border_CW_SW = (y>-2.3e6) & (x<-2e5)
border_CW = border_CW_NW & border_CW_SW
# NW
border_NW_NO = (y<-1.15e6) & (x<-2e5)
border_NW_CW = (y>-2e6) & (x<-2e5)
border_NW = border_NW_NO & border_NW_CW
# North: 
border_NO_NW = (y>-1.15e6)  & (x <0)
border_NO_NE = (y>-0.85e6) & (x<3.5e5)
border_NO = border_NO_NW | border_NO_NE
# NE
border_NE_NO = (x>3.5e5) 
border_NE_CE = (y>-1.85e6) 
border_NE = border_NE_NO & border_NE_CE

mask_1k_with_LSM['Basins_ocean'] = 1-  mask_1k_with_LSM['LSM']
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_NE), 2)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_CE), 3)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_SE), 4)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_SW), 5)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_CW), 6)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_NW), 7)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(~(border_NO), 1)
mask_1k_with_LSM['Basins_ocean'] = mask_1k_with_LSM['Basins_ocean'].where(mask_1k_with_LSM['LSM']==0, 0)


# now export to 
path_mask_1k_with_LSM = "/Users/annek/Documents/RACMO2.3p2/FGRN055/Downscaling_GR/GrIS_topo_icemask_lsm_tundra_basins_ocean_lon_lat_1km.nc"
mask_1k_with_LSM.to_netcdf(path_mask_1k_with_LSM)