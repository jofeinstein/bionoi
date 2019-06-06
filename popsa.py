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


def extract_sasa_data(mol,pop):

    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['subst_name']]
    atoms.columns = ['residue_type']
    siteresidue_list = atoms['residue_type'].tolist()

    residue_list = []
    qsasa_list = []
    with open(pop) as popsa:
        for line in popsa:
            list = line.split()
            if len(list) == 12:
                residue_type = list[2]+list[4]
                qsasa = list[7]
                residue_list.append(residue_type)
                qsasa_list.append(qsasa)
        residue_list = residue_list[1:]
        qsasa_list = qsasa_list[1:]

    fullprotein_data = dict(zip(residue_list,qsasa_list))
    qsasa_data = {k: float(fullprotein_data[k]) for k in siteresidue_list if k in fullprotein_data}

    return qsasa_data

def custom_colormap(color_scale):
    '''takes two hex colors and creates a linear colormap'''
    if color_scale == "red_cyan":
        colorlist = ("#ff0000","#00ffff")
    elif color_scale == "orange_bluecyan":
        colorlist = ("#ff7f00","#007fff")
    elif color_scale == "yellow_blue":
        colorlist = ("#ffff00","#0000ff")
    elif color_scale == "greenyellow_bluemagenta":
        colorlist = ("#7fff00","#7f00ff")
    elif color_scale == "green_magenta":
        colorlist = ("#00ff00","#ff00ff")
    elif color_scale == "greencyan_redmagenta":
        colorlist = ("#00ff7f","#ff007f")

    try:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap1', colorlist, N=256)
    except:
        cmap = None

    return cmap


def normalizer(dataset,colorby):
    '''normalizes data depending on the data given'''

    valnorm_lst = []

    #relative normalization
    if colorby in ["hydrophobicity","binding_prob","center_dist"]:
        for val in dataset.values():
            val = float(val)
            valnorm = ((val-min(dataset.values()))/(max(dataset.values())-min(dataset.values())))
            valnorm_lst.append(valnorm)

    #normalization based on given min/max values
    elif colorby == "charge":
        for val in dataset.values():
            val = float(val)
            valnorm = ((val + 0.4807) / 1.02)
            valnorm_lst.append(valnorm)

    elif colorby == "sasa":
        valnorm_lst = dataset.values()

    return valnorm_lst


def colorgen(colorby,valnorm_lst,cmap,dataset):
    '''dictionary of data -> maps normalized values to colormap'''

    if colorby in ["atom_type", "residue_type"]:
        color_map = "./cmaps/atom_cmap.csv" if colorby == "atom_type" else "./cmaps/res_hydro_cmap.csv"

        # Check for color mapping file, make dict
        with open(color_map, "rt") as color_mapF:
            # Parse color map file
            color_map = np.array(
                [line.replace("\n", "").split(";") for line in color_mapF.readlines() if not line.startswith("#")])
            # To dict
            color_map = {code: {"color": color, "definition": definition} for code, definition, color in color_map}
            return color_map
    else:
        color_lst = []

        #apply colormap
        for val in valnorm_lst:
            color = cmap(val)
            color = matplotlib.colors.rgb2hex(color)
            color_lst.append(color)

        #create color dictionary
        color_map = dict(zip(dataset.keys(), color_lst))
        names = ['code', 'color']
        dtype = dict(names=names)
        hexcolor_array = np.asarray(list(color_map.items()))
        color_map = {code: {"color": color} for code, color in hexcolor_array}
        return color_map

colorby = "sasa"
dataset = extract_sasa_data("./mol/5upxA00.mol2", "./popsa/5upxA.pdb")
cmap = custom_colormap("red_cyan")
valnorm_lst = normalizer(dataset,colorby)
color_map = colorgen("sasa",valnorm_lst,cmap,dataset)
print(color_map)

pd.options.mode.chained_assignment = None

# Read molecules in mol2 format
mol2 = PandasMol2().read_mol2("./mol/5upxA00.mol2")
atoms = mol2.df[['atom_id', 'subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']]
atoms.columns = ['atom_id', colorby, 'atom_type', 'atom_name', 'x', 'y', 'z', 'relative_charge']
atoms['atom_id'] = atoms['atom_id'].astype(str)
if colorby != "sasa":
    atoms[colorby] = atoms[colorby].apply(lambda x: x[0:3])


if colorby in ["charge", "center_dist"]:
    colors = [color_map[_type]["color"] for _type in atoms['atom_id']]
else:
    colors = [color_map[_type]["color"] for _type in atoms[colorby]]
atoms["color"] = colors
print(atoms.head())