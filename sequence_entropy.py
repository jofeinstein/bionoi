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

pd.options.mode.chained_assignment = None
mol2 = PandasMol2().read_mol2(mol)
atoms = mol2.df[['subst_name']]
atoms.columns = ['residue_type']
siteresidue_list = atoms['residue_type'].tolist()

residue_list = []
entropy_list = []
with open(pop) as profile:
    for line in profile:
        line_list = line.split()

        residue_type = line_list[0]
        prob_data = line_list[1:-1]
        print(prob_data)
        print(residue_type)


'''fullprotein_data = dict(zip(residue_list,qsasa_list))
qsasa_data = {k: float(fullprotein_data[k]) for k in siteresidue_list if k in fullprotein_data}'''