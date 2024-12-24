import numpy as np
from scipy.spatial import ConvexHull, Delaunay
import matplotlib.pyplot as plt
import xarray as xr
import shapely.geometry as geometry
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from tqdm import tqdm

from paths import *
from concave_hull import concave_hull, concave_hull_indexes
from find_largest_islands import *
from open_preprocess_racmo import *

mode = "including_islands"

ofolder = f'{pathIMAU02}RACMO2.3p2/FGRN055/Downscaling_GR/'

ds_no_islands = xr.open_dataset(os.path.join(pathDataTemp, "Mask_1km_largest_islands.nc"))
masks1k = xr.open_dataset(f'{pathIMAU02}RACMO2.3p2/FGRN055/Downscaling_GR/masks1k_copy.nc')
masks1k['LSM_no_islands'] = masks1k['LSM'].copy(data=ds_no_islands['LSM_no_islands'].values)
masks1k['LSM_Qeqertarsuaq'] = masks1k['LSM_no_islands'].copy(data=ds_no_islands['LSM_Qeqertarsuaq'].values)


if mode == "Normal":

    da = masks1k['LSM_no_islands']
    coords = np.argwhere(da.values == 1)
    points = coords

    print("Finding edge points")
    idxes = concave_hull_indexes(
        points[:, :2],
        length_threshold=10,
        concavity=2
    )

    print("Checking which points are inside the polygon")
    points_in_polygon = []
    polygon = geometry.Polygon(points[idxes])
    mask_greenland_boundary= np.zeros_like(da.values)
    x,y = np.meshgrid(np.arange(len(da.x)), np.arange(len(da.y)))
    y,x  = x.flatten(), y.flatten()
    for i in tqdm(np.arange(len(x))):
        point = geometry.Point(x[i], y[i])
        if polygon.contains(point):
            points_in_polygon.append([x[i], y[i]])
            mask_greenland_boundary[int(x[i]), int(y[i])] = 1

    masks1k['fjords'] = da.copy(deep=True) 
    masks1k['fjords'].values = mask_greenland_boundary - masks1k['LSM_no_islands'].values
    masks1k['fjords'].values = masks1k['fjords'].values.clip(min=0)
    masks1k['ocean'] = da.copy(deep=True)
    masks1k['ocean'].values = 1 - masks1k['LSM']
    masks1k['fjords'] = masks1k['fjords']*masks1k['ocean']



# #%% =============================================================================
#  now do the same but for Qeqertarsuaq
# =============================================================================

if mode == "Qeqertarsuaq":
    da = masks1k['LSM_Qeqertarsuaq']
    coords = np.argwhere(da.values == 1)
    points = coords

    print("Finding edge points")
    idxes = concave_hull_indexes(
        points[:, :2],
        length_threshold=10,
        concavity=2
    )

    print("Checking which points are inside the polygon")
    points_in_polygon = []
    polygon = geometry.Polygon(points[idxes])
    mask_greenland_boundary= np.zeros_like(da.values)
    x,y = np.meshgrid(np.arange(len(da.x)), np.arange(len(da.y)))
    y,x  = x.flatten(), y.flatten()
    for i in tqdm(np.arange(len(x))):
        point = geometry.Point(x[i], y[i])
        if polygon.contains(point):
            points_in_polygon.append([x[i], y[i]])
            mask_greenland_boundary[int(x[i]), int(y[i])] = 1

    fjords_grid_Qeqertarsuaq = mask_greenland_boundary - masks1k['LSM_Qeqertarsuaq'].values
    masks1k['fjords_Qeqertarsuaq'] = da.copy(data=fjords_grid_Qeqertarsuaq)
    masks1k['fjords_Qeqertarsuaq'].values = masks1k['fjords_Qeqertarsuaq'].values.clip(min=0)
    masks1k['ocean'] = da.copy(deep=True)
    masks1k['ocean'].values = 1 - masks1k['LSM']
    masks1k['fjords_Qeqertarsuaq'] = masks1k['fjords_Qeqertarsuaq']*masks1k['ocean']
    masks1k['fjords_Qeqertarsuaq'] = masks1k['fjords_Qeqertarsuaq'].fillna(0)
    masks1k['fjords_Qeqertarsuaq'].to_netcdf(ofolder + 'fjords_Qeqertarsuaq.nc')
    print("Checked for every point if it is inside the polygon and saved the mask to fjords.nc")


if mode == 'including_islands':
    print("Including islands")
    ofile = "fjords_incl_islands.nc"
    da = masks1k['LSM']
    coords = np.argwhere(da.values == 1)
    points = coords

    print("Finding edge points")
    idxes = concave_hull_indexes(
        points[:, :2],
        length_threshold=10,
        concavity=2
    )

    print("Checking which points are inside the polygon")
    points_in_polygon = []
    polygon = geometry.Polygon(points[idxes])
    mask_greenland_boundary= np.zeros_like(da.values)
    x,y = np.meshgrid(np.arange(len(da.x)), np.arange(len(da.y)))
    y,x  = x.flatten(), y.flatten()
    for i in tqdm(np.arange(len(x))):
        point = geometry.Point(x[i], y[i])
        if polygon.contains(point):
            points_in_polygon.append([x[i], y[i]])
            mask_greenland_boundary[int(x[i]), int(y[i])] = 1

    fjords_grid= mask_greenland_boundary - masks1k['LSM'].values
    masks1k['fjords_incl_islands'] = da.copy(data=fjords_grid)
    masks1k['fjords_incl_islands'].values = masks1k['fjords_incl_islands'].values.clip(min=0)   
    masks1k['ocean'] = da.copy(deep=True)
    masks1k['ocean'].values = 1 - masks1k['LSM']
    masks1k['fjords_incl_islands'] = masks1k['fjords_incl_islands']*masks1k['ocean']
    masks1k['fjords_incl_islands'] = masks1k['fjords_incl_islands'].fillna(0)

    masks1k['fjords_incl_islands'].to_netcdf(ofile)
    print("Checked for every point if it is inside the polygon and saved the mask to "+ ofile)
