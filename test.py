from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import pandas as pd
import matplotlib

from biopandas.mol2 import PandasMol2
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from sklearn.cluster import KMeans
from math import sqrt, asin, atan2, log, pi, tan

from alignment import align
# Read molecules in mol2 format
pd.options.mode.chained_assignment = None
mol2 = PandasMol2().read_mol2('./mol/3nbfC02-1.mol2')
atoms = mol2.df[['subst_id', 'subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z']]
atoms.columns = ['res_id', 'residue_type', 'atom_type', 'atom_name', 'x', 'y', 'z']
atoms['residue_type'] = atoms['residue_type'].apply(lambda x: x[0:3])
print(atoms)

dataset = {'ALA':1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS':2.5,'GLN':-3.5,'GLU':-3.5,
       'GLY':-0.4,'HIS':-3.2,'ILE':4.5,'LEU':3.8,'LYS':-3.9,'MET':1.9,'PHE':2.8,
       'PRO':-1.6,'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL':4.2}

#red_cyan example colormap
red_cyan = ((255,0,0),(0,255,255))
colorarray = np.asarray(red_cyan)
colorarray = colorarray/255
red_cyan = np.array(colorarray).tolist()
a=np.outer(np.arange(0,1,0.001),np.ones(500))
red_cyan_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('red_cyan_cmap', red_cyan, N=256)

def RGB_to_hex(RGB):
    '''(RGB) -> #FFFFFF'''
    RGB = [int(x) for x in RGB]
    return "#" + "".join(["0{0:x}".format(v).upper() if v < 16 else
                          "{0:x}".format(v).upper() for v in RGB])

def colorgen(data,cmap):
    '''dictionary of data -> maps normalized values to colormap'''

    #normalize data
    valnorm_lst = list()
    rgbcolor_lst = list()
    for val in data.values():
        val = float(val)
        valnorm = ((val-min(data.values()))/(max(data.values())-min(data.values())))
        valnorm_lst.append(valnorm)

    #apply colormap
    for val in valnorm_lst:
        color = cmap(val)
        color = color[0:3]
        rgbcolor_lst.append(color)

    rgb_array = np.asarray(rgbcolor_lst)*255
    rgbcolor_lst = np.array(rgb_array).tolist()

    #convert rgb to hex
    hexcolor_lst = list()
    for color in rgbcolor_lst:
        hex = RGB_to_hex(color)
        hexcolor_lst.append(hex)

    #create color dictionary
    color_map = dict(zip(data.keys(), hexcolor_lst))
    names = ['code', 'color']
    dtype = dict(names=names)
    hexcolor_array = np.asarray(list(color_map.items()))
    colormap = {code: {"color": color} for code, color in hexcolor_array}
    colors = [color_map[_type]["color"] for _type in atoms['residue_type']]

    return colors

print(colorgen(dataset,red_cyan_cmap))








#color_map = {code: {"color": color} for code, color in hexcolorarray}
#print(color_map)


#colors = [color_map["color"] for _type in atoms['residue_type']]
#print(atoms['residue_type'])

#colorgen(dataset,red_cyan_cmap)


'''plt.imshow(a,cmap=colormap)
    plt.axis('off')
    plt.show()
    plt.savefig('{}.png'.format(color), bbox_inches="tight",transparent="True", pad_inches=0)
    plt.close()'''

#other usable colors
'''#opposite colors
red_cyan = ((255,0,0),(0,255,255))
orange_bluecyan = ((255,127,0),(0,127,255))
yellow_blue = ((255,255,0),(0,0,255))
greenyellow_bluemagenta = ((127,255,0),(127,0,255))
green_magenta = ((0,255,0),(255,0,255))
greencyan_redmagenta = ((0,255,127),(255,0,127))

#neighboring colors
red_orange = ((255,0,0),(255,127,0))
yellow_yellowgreen = ((255,255,0),(127,255,0))
green_greencyan = ((0,255,0),(0,255,127))
cyan_bluecyan = ((0,255,255),(0,127,255))
blue_bluemagenta = ((0,0,255),(127,0,255))
magenta_redmagenta = ((255,0,255),(255,0,127))

colorlist = (red_cyan,orange_bluecyan,yellow_blue,
             greenyellow_bluemagenta,green_magenta,greencyan_redmagenta,
             red_orange,yellow_yellowgreen,green_greencyan,cyan_bluecyan,
             blue_bluemagenta,magenta_redmagenta)

'''
codes = ['ALA','ARG','ASN','ASP','CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
definitions = ['ALA','ARG','ASN','ASP','CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
