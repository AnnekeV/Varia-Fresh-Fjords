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
# cmc.roma.colors
# get 8 colors from the roma colormap
colors = cmc.roma.colors
colors = [colors[i] for i in range(0, len(colors), len(colors)//8)]
colors_hex_Moug = [mcolors.rgb2hex(color) for color in colors]
dict_Moug_colors = dict(zip(np.arange(8), colors))



dict_Moug = {0:"Not GRIS", 1: "NO", 2: "NE", 3: "CE", 4: "SE", 5: "SW", 6: "CW", 7: "NW"}
inverted_dict_Moug = {v: k for k, v in dict_Moug.items()}

dict_regionName_colors = {k: dict_region_colors[v] for k, v in inverted_dict_Moug.items()}

dict_sections = dict(zip(np.arange(1,8),['SE', 'SW', 'CE', 'CW', 'NE', 'NW', 'NO']))
# {1: 'SE', 2: 'SW', 3: 'CE', 4: 'CW', 5: 'NE', 6: 'NW', 7: 'NO'}

dict_linestyle_FW = {'Liquid Runoff Ice Sheet': '-',
 'Liquid Runoff Ice Caps': '-',
 'Liquid Runoff Tundra': '-',
 'Precipitation Fjords': '-',
 'Solid Ice discharge (weighted mean)': '-',
 'Solid Ice Discharge': '-',
 'Precipitation Fjords CARRA': ':'}

 # change labels of legend from Solid Ice discharge (weighted mean) to Solid Ice discharge 
dict_shorter_name = {'Solid Ice discharge (weighted mean)': 'Solid Ice Discharge'}

# for legend
# plot ds solid discharge
col_order_abs = ['Solid Ice Discharge', 'Precipitation Fjords',  'Liquid Runoff Tundra', 'Liquid Runoff Ice Caps','Liquid Runoff Ice Sheet', ]
col_order_rel = ['Solid Ice Discharge',   'Liquid Runoff Ice Sheet', 'Liquid Runoff Ice Caps','Liquid Runoff Tundra','Precipitation Fjords', ]

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
FW_type_long = ['Discharge [Gt yr-1]', 'Runoff GrIS', 'Runoff Tundra',
       'Precipitation Fjords', 'Runoff Ice Caps']
dictFWcolor_long = dict(zip(FW_type_long, colors))
dictFWcolor_long['Runoff Ice Caps'] = 'plum'

FW_type = ["Solid", "IceRun", "Tundra", "Precip", "IceCap"]
dictFWcolor = dict(zip(FW_type, colors))
dictFWcolor['IceCap'] = 'plum'
