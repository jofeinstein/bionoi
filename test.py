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
mol2 = PandasMol2().read_mol2("./mol/3nbfC02-1.mol2")
atoms = mol2.df[['charge']]
atoms.columns = ['charge']
charge_list = atoms['charge'].tolist()
atoms['charge'] = atoms['charge'].astype(str)
charge_data = dict(zip(charge_list, charge_list))


orange_bluecyan = ((255,127,0),(0,127,255))
colorarray = np.asarray(orange_bluecyan)
colorarray = colorarray/255
orange_bluecyan = np.array(colorarray).tolist()
orange_bluecyan_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('orange_bluecyan_cmap', orange_bluecyan, N=256)

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
    colors = [colormap[_type]["color"] for _type in atoms['charge']]

    return colors



'''plt.imshow(a,cmap=colormap)
    plt.axis('off')
    plt.show()
    plt.savefig('{}.png'.format(color), bbox_inches="tight",transparent="True", pad_inches=0)
    plt.close()'''



def normalizer(data,colorby):
    #normalize data
    valnorm_lst = []
    if colorby == "residue_type":
        for val in data.values():
            val = float(val)
            valnorm = ((val-min(data.values()))/(max(data.values())-min(data.values())))
            valnorm_lst.append(valnorm)
    elif colorby == "charge":
        for val in data.values():
            val = float(val)
            valnorm = ((val + 0.4807) / 1.0)
            print(valnorm)
            valnorm_lst.append(valnorm)
    return
print(normalizer(charge_data,"charge"))

def colorgen(valnorm_lst,cmap,data):
    '''dictionary of data -> maps normalized values to colormap'''
    rgbcolor_lst = []
    #apply colormap
    for val in valnorm_lst:
        color = cmap(val)
        color = color[0:3]
        rgbcolor_lst.append(color)

    rgb_array = np.asarray(rgbcolor_lst)*255
    rgbcolor_lst = np.array(rgb_array).tolist()

    #convert rgb to hex
    hexcolor_lst = []
    for color in rgbcolor_lst:
        hex = RGB_to_hex(color)
        hexcolor_lst.append(hex)

    #create color dictionary
    color_map = dict(zip(data.keys(), hexcolor_lst))
    names = ['code', 'color']
    dtype = dict(names=names)
    hexcolor_array = np.asarray(list(color_map.items()))
    color_map = {code: {"color": color} for code, color in hexcolor_array}
    return color_map