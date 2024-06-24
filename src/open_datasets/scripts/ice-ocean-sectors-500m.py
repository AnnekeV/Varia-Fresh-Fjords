import pandas as pd
import numpy as np
import xarray as xr
import scipy.io
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import geopandas as gpd


path_run_gemin = '/science/projects/imau02/rapid/Anneke/' 
folder500m = path_run_gemin + "RACMO2.3p2/FGRN055/Downscaling_GR_500m/"

fileRACMO = folder500m + "Annual/runoff.1940-2023.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.GrIS.0.5km.YY.nc"
mask_500m = folder500m + "GrIS_topo_icemask_lsm_lon_lat_0.5km.nc"
PathSectorsIceOceanPath = path_run_gemin+ "Other/fjords_GRL_2022/ice_ocean_sectors.mat"
shapefile_CE = path_run_gemin+"Other/new_CE_new_proj.shp"



dsBrice = xr.open_dataset(mask_500m)

# Load the data from the .mat file
data = scipy.io.loadmat(PathSectorsIceOceanPath)

# Extract the 'ice' data
ice_data = data['regions']['ice'][0][0]

# Extract the 'x' and 'y' coordinates
x_ice = ice_data['x'][0][0]
y_ice = ice_data['y'][0][0]


# Repeat the process for the 'ocean' data
ocean_data = data['regions']['ocean'][0][0]
x_ocean = ocean_data['x'][0][0]
y_ocean = ocean_data['y'][0][0]




section_numbers = xr.DataArray(np.full((len(dsBrice['y'].values), len(dsBrice['x'].values)), np.nan), 
                               coords={'y': dsBrice['y'].values, 'x': dsBrice['x'].values}, 
                               dims=['y', 'x'])

# Create a list of polygons for the ocean sections
ocean_polygons = [Polygon(list(zip(np.squeeze(section['x'][0][0]), np.squeeze(section['y'][0][0])))) for section in ocean_data]

# Create a list of points from the 'x' and 'y' coordinates in dsBrice
points = [(x, y) for x in dsBrice['x'].values for y in dsBrice['y'].values]



# Loop over all points
for point_tuple in tqdm(points):
    point_obj = Point(point_tuple)

    # Check which ocean section the point is in
    for i, polygon in enumerate(ocean_polygons):
        if polygon.contains(point_obj):
            # Assign the section number to the corresponding location in the data array
            section_numbers.loc[dict(x=point_tuple[0], y=point_tuple[1])] = i + 1


section_numbers.to_dataset(name='section_numbers_slater').to_netcdf(os.path.join(folder500m, 'section_numbers_slater_0.5km.nc'))


# Load the shapefile
gdf = gpd.read_file(shapefile_CE)
gdf_CE_extend = gdf
# Print the geometry
print(gdf.geometry)




# Create a copy of section_numbers
adjusted_sections = section_numbers.copy()
# Update the values in adjusted_sections with the values in extended_CE
adjusted_sections.values[extended_CE.values > 0] = 3
adjusted_sections.values[daMougCE_Ice.values==3] = 3

adjusted_sections.to_dataset(name='section_numbers_adjusted').to_netcdf(os.path.join(folder500m, 'adjusted_section_numbers_slater_0.5km.nc'))