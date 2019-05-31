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