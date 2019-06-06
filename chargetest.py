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
def alignment(pocket, proj_direction):
    """Principal Axes Alignment
    Returns transformation coordinates(matrix: X*3)"""
    pocket_coords = np.array([pocket.x, pocket.y, pocket.z]).T
    pocket_center = np.mean(pocket_coords, axis=0)  # calculate mean of each column
    pocket_coords = pocket_coords - pocket_center  # Centralization
    inertia = np.cov(pocket_coords.T)  # get covariance matrix (of centralized data)
    e_values, e_vectors = np.linalg.eig(inertia)  # linear algebra eigenvalue eigenvector
    sorted_index = np.argsort(e_values)[::-1]  # sort eigenvalues (increase)and reverse (decrease)
    sorted_vectors = e_vectors[:, sorted_index]

    transformation_matrix = align(sorted_vectors, proj_direction)
    transformed_coords = (np.matmul(transformation_matrix, pocket_coords.T)).T

    return transformed_coords
def extract_charge_data(mol):
    '''extracts and formats charge data from mol2 file'''

    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['x', 'y', 'z', 'charge']]
    atoms.columns = ['x', 'y', 'z', 'charge']
    atoms['position'] = atoms['x'].astype(str) + atoms['y'].astype(str) + atoms['z'].astype(str)
    charge_list = atoms['charge'].tolist()
    position_list = atoms['position'].tolist()
    charge_data = dict(zip(position_list, charge_list))

    return charge_data

def extract_centerdistance_data(mol,proj_direction):
    '''extracts and formats center distance from mol2 file after alignment to principal axes'''

    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['x', 'y', 'z']]
    atoms.columns = ['x', 'y', 'z']
    atoms['position'] = atoms['x'].astype(str) + atoms['y'].astype(str) + atoms['z'].astype(str)
    trans_coords = alignment(atoms, proj_direction)  # get the transformation coordinate
    atoms['x'] = trans_coords[:, 0]
    atoms['y'] = trans_coords[:, 1]
    atoms['z'] = trans_coords[:, 2]

    position_list = atoms['position'].tolist()
    coordinate_list = atoms.values.tolist()

    center_dist_list = []
    for xyz in coordinate_list:
        center_dist = ((xyz[0]) ** 2 + (xyz[1]) ** 2 + (xyz[2]) ** 2) ** .5
        center_dist_list.append(center_dist)
    center_dist_data = dict(zip(position_list, center_dist_list))

    return center_dist_data

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

def normalizer(dataset,colorby):
    '''normalizes data depending on the data given'''

    valnorm_lst = []

    #relative normalization
    if colorby in ["hydrophobicity","binding_prob","center_distance"]:
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
    return valnorm_lst

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

colorby = "charge"
dataset = extract_charge_data("./mol/4v94E.mol2")
valnorm_lst = normalizer(dataset,colorby)
cmap = custom_colormap("red_cyan")
color_map = colorgen(colorby,valnorm_lst,cmap,dataset)
print(color_map)

pd.options.mode.chained_assignment = None

# Read molecules in mol2 format
mol2 = PandasMol2().read_mol2("./mol/4v94E.mol2")
atoms = mol2.df[['atom_id', 'subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']]
atoms.columns = ['atom_id', colorby, 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']
atoms[colorby] = atoms[colorby].apply(lambda x: x[0:3])
atoms['position'] = atoms['x'].astype(str) + atoms['y'].astype(str) + atoms['z'].astype(str)

'''atoms = mol2.df[['subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']]
atoms.columns = [colorby, 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']
atoms[colorby] = atoms[colorby].apply(lambda x: x[0:3])
atoms['position'] = atoms['x'].astype(str) + atoms['y'].astype(str) + atoms['z'].astype(str)'''

if colorby in ["charge", "center_distance"]:
    colors = [color_map[_type]["color"] for _type in atoms['position']]
else:
    colors = [color_map[_type]["color"] for _type in atoms[colorby]]
atoms["color"] = colors

print(atoms.head())