
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import xarray as xr
import shapely.geometry as geometry
import pylab as pl
import numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
from tqdm import tqdm
from open_preprocess_racmo import *
from paths import *
time_resolution = "Monthly"

print("Reading RACMO 1km")

dsPrecipMonth = read_RACMO("1k",time_resolution, ['1990', '2022'], variable="precip")
masks1k = open_mask_1k()
masks1k['x'] = dsPrecipMonth['x']
masks1k['y'] = dsPrecipMonth['y']

print("Selecting fjords")
dsPrecipFjords = dsPrecipMonth["precipcorr"].where(masks1k["fjords"], np.nan)
print("Saving to: ")
print(pathDataTemp + "RACMO2.3p2_1km_precip_fjords.nc")
dsPrecipFjords.to_netcdf(pathDataTemp + "RACMO2.3p2_1km_precip_fjords.nc")