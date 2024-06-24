dict_Moug = {0:"Not GRIS", 1: "NO", 2: "NE", 3: "CE", 4: "SE", 5: "SW", 6: "CW", 7: "NW"}

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmcrameri.cm as cmc
import numpy as np


# make a list of all the names of data variables in dsSectorSum
data_varsFW = ['Liquid Runoff Ice Sheet',
 'Liquid Runoff Ice Caps',
 'Liquid Runoff Tundra',
 'Precipitation Fjords',
 'Solid Ice discharge (weighted mean)']
colors_FW = ['tab:orange', 'plum', 'tab:green', 'tab:red', 'tab:blue']

dict_colors_FW = dict(zip(data_varsFW, colors_FW))
dict_colors_FW['Solid Ice Discharge'] = 'tab:blue'
dict_colors_FW['Precipitation Fjords CARRA'] = 'tab:red'

# Convert the colors to RGB format and store them in the new dictionary
dict_colors_FW_rgb = {}
for source, color in dict_colors_FW.items():
    rgb = mcolors.to_rgb(color)
    dict_colors_FW_rgb[source] = rgb




colors = cmc.batlowS.colors
colors = [colors[i] for i in range(0, len(colors), len(colors)//8)]
colors_hex_region = [mcolors.rgb2hex(color) for color in colors]
dict_Moug_colors_region = dict(zip(np.arange(8), colors))
dict_region_colors = dict_Moug_colors_region


dict_Moug = {0:"Not GRIS", 1: "NO", 2: "NE", 3: "CE", 4: "SE", 5: "SW", 6: "CW", 7: "NW"}
inverted_dict_Moug = {v: k for k, v in dict_Moug.items()}

dict_regionName_colors = {k: dict_region_colors[v] for k, v in inverted_dict_Moug.items()}

dict_sections = dict(zip(np.arange(1,8),['SE', 'SW', 'CE', 'CW', 'NE', 'NW', 'NO']))
# {1: 'SE', 2: 'SW', 3: 'CE', 4: 'CW', 5: 'NE', 6: 'NW', 7: 'NO'}