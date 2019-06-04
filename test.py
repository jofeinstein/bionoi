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

def colorby_con(colorby):
    if colorby in ["hydrophobicity", "binding_prob"]:
        colorby = colorby
    else:
        colorby = "residue_type"
    return colorby




def custom_colormap(color_scale):
    #takes two hex colors and creates a linear colormap
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
    elif color_scale == "red_orange":
        colorlist = ("#ff0000","#ff7f00")
    elif color_scale == "yellow_yellowgreen":
        colorlist = ("#ffff00","#7fff00")
    elif color_scale == "green_greencyan":
        colorlist = ("#00ff00","#00ff7f")
    elif color_scale == "cyan_bluecyan":
        colorlist = ("#00ffff","#007fff")
    elif color_scale == "blue_bluemagenta":
        colorlist = ("#0000ff","#7f00ff")
    elif color_scale == "magenta_redmagenta":
        colorlist = ("#ff00ff","#ff007f")

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap1', colorlist, N=256)

    return cmap


def normalizer(dataset,colorby):
    #normalize data
    valnorm_lst = []
    if colorby in ["hydropathicity","binding_probability"]:
        for val in dataset.values():
            val = float(val)
            valnorm = ((val-min(dataset.values()))/(max(dataset.values())-min(dataset.values())))
            valnorm_lst.append(valnorm)
    elif colorby == "charge":
        for val in dataset.values():
            val = float(val)
            valnorm = ((val + 0.4807) / 1.02)
            valnorm_lst.append(valnorm)
    return valnorm_lst


def colorgen(colorby,valnorm_lst,cmap,dataset):
    #dictionary of data -> maps normalized values to colormap
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
def extract_charge_data(mol):
    # Read molecules in mol2 format
    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['charge']]
    atoms.columns = ['charge']
    charge_list = atoms['charge'].tolist()
    atoms['charge'] = atoms['charge'].astype(str)
    charge_data = dict(zip(charge_list, charge_list))

    return charge_data

'''def colorby_converter(colorby):
    if colorby in ["hydropathicity", "binding_probability"]:
        if colorby == "hydropathicity":
            colorby = "residue_type"
        elif colorby == "binding_probability":
            colorby = "residue_type"
        elif colorby == "charge":
            colorby = "charge"
        return colorby'''


binding_probability_data = {'ALA':0.701,'ARG':0.916,'ASN':0.811,'ASP':1.015,
                             'CYS':1.650,'GLN':0.669,'GLU':0.956,'GLY':0.788,
                             'HIS':2.286,'ILE':1.006,'LEU':1.045,'LYS':0.468,
                             'MET':1.894,'PHE':1.952,'PRO':0.212,'SER':0.883,
                             'THR':0.730,'TRP':3.084,'TYR':1.672,'VAL':0.884}


colorby = "binding_prob"

pd.options.mode.chained_assignment = None
mol2 = PandasMol2().read_mol2("./mol/3nbfC02-1.mol2")
atoms = mol2.df[['subst_id', 'subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']]
atoms.columns = ['res_id', 'residue_type', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']
atoms['residue_type'] = atoms['residue_type'].apply(lambda x: x[0:3])
if colorby == "hydrophobicity":
    atoms.columns = ['res_id', 'hydrophobicity', 'atom_type', 'atom_name', 'x', 'y', 'z','charge']
elif colorby == "binding_prob":
    atoms.columns = ['res_id', 'binding_prob', 'atom_type', 'atom_name', 'x', 'y', 'z','charge']
atoms['charge'] = atoms['charge'].astype(str)

valnorm_lst = normalizer(binding_probability_data,colorby)
cmap = custom_colormap("red_cyan")
color_map = colorgen("binding_prob",valnorm_lst,cmap,binding_probability_data)
colors = [color_map[_type]["color"] for _type in atoms[colorby]]
atoms["color"] = colors
