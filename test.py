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

hydropathicty_data = {'ALA':1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,
                      'CYS':2.5,'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,
                      'HIS':-3.2,'ILE':4.5,'LEU':3.8,'LYS':-3.9,
                      'MET':1.9,'PHE':2.8,'PRO':-1.6,'SER':-0.8,
                      'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL':4.2}

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


def extract_centerdistance_data(mol,proj_direction):
    '''extracts and formats center distance from mol2 file after alignment to principal axes'''

    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['atom_id', 'x', 'y', 'z']]
    atoms.columns = ['atom_id', 'x', 'y', 'z']
    trans_coords = alignment(atoms, proj_direction)  # get the transformation coordinate
    atoms['x'] = trans_coords[:, 0]
    atoms['y'] = trans_coords[:, 1]
    atoms['z'] = trans_coords[:, 2]

    atomid_list = atoms['atom_id'].tolist()
    coordinate_list = atoms.values.tolist()

    center_dist_list = []
    for xyz in coordinate_list:
        center_dist = ((xyz[1]) ** 2 + (xyz[2]) ** 2 + (xyz[3]) ** 2) ** .5
        center_dist_list.append(center_dist)
    center_dist_data = dict(zip(atomid_list, center_dist_list))

    return center_dist_data

def extract_charge_data(mol):
    '''extracts and formats charge data from mol2 file'''

    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['atom_id','charge']]
    atoms.columns = ['atom_id','charge']
    charge_list = atoms['charge'].tolist()
    atomid_list = atoms['atom_id'].tolist()
    atoms['charge'] = atoms['charge'].astype(str)
    charge_data = dict(zip(atomid_list, charge_list))

    return charge_data

def normalizer(dataset,colorby):
    #normalize data
    valnorm_lst = []
    if colorby in ["hydrophobicity","binding_probability","center_distance"]:
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


colorby = "center_distance"
dataset = extract_centerdistance_data("./mol/3nbfC02-1.mol2",proj_direction=1)
valnorm_lst = normalizer(dataset,colorby)
cmap = custom_colormap("red_cyan")
color_map = colorgen(colorby,valnorm_lst,cmap,dataset)
pd.options.mode.chained_assignment = None
print(color_map)
    # Read molecules in mol2 format
mol2 = PandasMol2().read_mol2("./mol/3nbfC02-1.mol2")
atoms = mol2.df[['atom_id', 'subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z','charge']]
if colorby == "hydrophobicity":
    atoms.columns = ['atom_id', 'hydrophobicity', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']
    atoms['hydrophobicity'] = atoms['hydrophobicity'].apply(lambda x: x[0:3])
elif colorby == "binding_prob":
    atoms.columns = ['atom_id', 'binding_prob', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']
    atoms['binding_prob'] = atoms['binding_prob'].apply(lambda x: x[0:3])
else:
    atoms.columns = ['atom_id', 'residue_type', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']
    atoms['residue_type'] = atoms['residue_type'].apply(lambda x: x[0:3])
atoms['atom_id']=atoms['atom_id'].astype(str)



if colorby in ["hydrophobicity","binding_prob"]:
    colors = [color_map[_type]["color"] for _type in atoms[colorby]]
else:
    colors = [color_map[_type]["color"] for _type in atoms['atom_id']]
atoms["color"] = colors
print(atoms.head())