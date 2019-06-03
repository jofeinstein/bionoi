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

def custom_colormap(colorscale):
    '''takes two hex colors and creates a linear colormap'''
    if colorscale == "red_cyan":
        colorlist = ("#ff0000","#00ffff")
    elif colorscale == "orange_bluecyan":
        colorlist = ("#ff7f00","#007fff")
    elif colorscale == "yellow_blue":
        colorlist = ("#ffff00","#0000ff")
    elif colorscale == "greenyellow_bluemagenta":
        colorlist = ("#7fff00","#7f00ff")
    elif colorscale == "green_magenta":
        colorlist = ("#00ff00","#ff00ff")
    elif colorscale == "greencyan_redmagenta":
        colorlist = ("#00ff7f","#ff007f")
    elif colorscale == "red_orange":
        colorlist = ("#ff0000","#ff7f00")
    elif colorscale == "yellow_yellowgreen":
        colorlist = ("#ffff00","#7fff00")
    elif colorscale == "green_greencyan":
        colorlist = ("#00ff00","#00ff7f")
    elif colorscale == "cyan_bluecyan":
        colorlist = ("#00ffff","#007fff")
    elif colorscale == "blue_bluemagenta":
        colorlist = ("#0000ff","#7f00ff")
    elif colorscale == "magenta_redmagenta":
        colorlist = ("#ff00ff","#ff007f")

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap1', colorlist, N=256)

    return cmap

red_cyan = custom_colormap("blue_bluemagenta")
x = (1,2,3,4,5)
for val in x:
    color = red_cyan(val)
    color = matplotlib.colors.rgb2hex(color)
    print(color)

    , "yellow_blue", "greenyellow_bluemagenta",
    "green_magenta", "greencyan_redmagenta", "red_orange", "yellow_yellowgreen",
    "green_greencyan", "cyan_bluecyan", "blue_bluemagenta", "magenta_redmagenta"