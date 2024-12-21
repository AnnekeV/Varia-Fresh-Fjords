from scripts.paths import *
import os


# Old
# mask_sections = (
#     pathIMAU02 + "RACMO2.3p2/FGRN055/Downscaling_GR/Mask_adjusted_section_numbers_slater_may24.nc")

# New
mask_sections = pathDataTemp + "adjusted_section_numbers_slater.nc"

path_mask_islands = os.path.join(pathDataTemp, "Mask_1km_largest_islands.nc")
pathmasks1k =  os.path.join(pathDataTemp, "Icemask_Topo_Iceclasses_lon_lat_average_1km.nc")
