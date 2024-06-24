import numpy as np
import xarray as xr
from tqdm import tqdm
import os

try:
    from scripts.paths import *
except:
    from paths import *

path_mask_islands = os.path.join(pathDataTemp, "Mask_1km_largest_islands.nc")
def largest_islands(grid):
    ''' Calculate the size of the largest and second largest island in a grid of 1's and 0's.
    Makes a mask of the biggest and second biggest island in the LSM mask, and adds it to the dataset as LSM_no_islands and LSM_Qeqertarsuaq.
    input is a dataset with a LSM mask, output is the same dataset with two new masks.
    '''
    def dfs(row, col):
        stack = [(row, col)]
        island_coords = []  # List to store island coordinates
        size = 0
        while stack:
            r, c = stack.pop()
            if 0 <= r < rows and 0 <= c < cols and grid[r][c] == 1 and not visited[r][c]:
                visited[r][c] = True
                island_coords.append([r, c])  # Add current coordinate to island_coords
                size += 1
                for dr, dc in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    stack.append((r + dr, c + dc))
        return size, island_coords

    rows, cols = len(grid), len(grid[0])
    visited = [[False] * cols for _ in range(rows)]
    max_size = 0
    max_island_coords = []
    second_max_size = 0
    second_max_island_coords = []

    for i in tqdm(range(rows)):
        for j in range(cols):
            if grid[i][j] == 1 and not visited[i][j]:
                island_size, island_coords = dfs(i, j)
                if island_size > max_size:
                    second_max_size = max_size
                    second_max_island_coords = max_island_coords
                    max_size = island_size
                    max_island_coords = island_coords
                elif island_size > second_max_size:
                    second_max_size = island_size
                    second_max_island_coords = island_coords

    return max_size, np.array(max_island_coords), second_max_size, np.array(second_max_island_coords)


if __name__ == "__main__":
    print("opening 1km masks")
    pathmasks1k =  os.path.join(pathdata, "Downscaling_GR",  "Icemask_Topo_Iceclasses_lon_lat_average_1km.nc")

    masks1k = xr.open_dataset(pathmasks1k, engine='netcdf4')


    print("opening 1km LSM mask")
    dsLSMmask1k = xr.open_dataset(pathlocal +"Downscaling_GR/GrIS_topo_icemask_lsm_lon_lat_1km.nc")


    '''Makes a mask of the biggest and second biggest island in the LSM mask, and adds it to the dataset as LSM_no_islands and LSM_Qeqertarsuaq.
    input is a dataset with a LSM mask, output is the same dataset with two new masks.
    '''

    print("calculating largest islands")
    grid = dsLSMmask1k['LSM'].values
    max_size, max_island_coords, second_max_size, second_max_island_coords = largest_islands(grid)
    max_island_coords = np.array(max_island_coords)

    print("assign main island")
    masks1k['LSM_no_islands'] = masks1k['Topography'].copy(deep=True)
    array_no_islands = np.zeros_like(masks1k['Topography'].values)
    array_no_islands[max_island_coords[:, 0], max_island_coords[:, 1]] = 1
    masks1k['LSM_no_islands'].values = array_no_islands

    print("assign second island Qeqertarsuaq")
    array_disko = np.zeros_like(masks1k['Topography'].values)
    array_disko[second_max_island_coords[:, 0], second_max_island_coords[:, 1]] = 1
    masks1k['LSM_Qeqertarsuaq'] = masks1k['Topography'].copy(deep=True)
    masks1k['LSM_Qeqertarsuaq'].values = array_disko

    # now only keep variables LSM_no_islands and LSM_Qeqertarsuaq
    masks1k = masks1k[['LSM_no_islands', 'LSM_Qeqertarsuaq']]

    print("saving to ",  path_mask_islands)

    masks1k.to_netcdf(path_mask_islands)
    print("done")



