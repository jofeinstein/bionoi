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


def k_different_colors(k: int):
    colors = dict(**mcolors.CSS4_COLORS)

    rgb = lambda color: mcolors.to_rgba(color)[:3]
    hsv = lambda color: mcolors.rgb_to_hsv(color)

    col_dict = [(k, rgb(k)) for c, k in colors.items()]
    X = np.array([j for i, j in col_dict])

    # Perform kmeans on rqb vectors
    kmeans = KMeans(n_clusters=k)
    kmeans = kmeans.fit(X)
    # Getting the cluster labels
    labels = kmeans.predict(X)
    # Centroid values
    C = kmeans.cluster_centers_

    # Find one color near each of the k cluster centers
    closest_colors = np.array([np.sum((X - C[i]) ** 2, axis=1) for i in range(C.shape[0])])
    keys = sorted(closest_colors.argmin(axis=1))

    return [col_dict[i][0] for i in keys]


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    Source
    -------
    Copied from https://gist.github.com/pv/8036995
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max() * 2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def fig_to_numpy(fig, alpha=1) -> np.ndarray:
    '''
    Converts matplotlib figure to a numpy array.

    Source
    ------
    Adapted from https://stackoverflow.com/questions/7821518/matplotlib-save-plot-to-numpy-array
    '''

    # Setup figure
    fig.patch.set_alpha(alpha)
    fig.canvas.draw()

    # Now we can save it to a numpy array.
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return data


def miller(x, y, z):
    radius = sqrt(x ** 2 + y ** 2 + z ** 2)
    latitude = asin(z / radius)
    longitude = atan2(y, x)
    lat = 5 / 4 * log(tan(pi / 4 + 2 / 5 * latitude))
    return lat, longitude


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


def voronoi_atoms(bs, color_map, colorby, bs_out=None, size=None, dpi=None, alpha=1, save_fig=True,
                  projection=miller, proj_direction=None):
    # Suppresses warning
    pd.options.mode.chained_assignment = None

    # Read molecules in mol2 format
    mol2 = PandasMol2().read_mol2(bs)
    atoms = mol2.df[['atom_id','subst_name', 'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']]
    atoms.columns = ['atom_id',colorby_conv(colorby), 'atom_type', 'atom_name', 'x', 'y', 'z', 'relative_charge']
    atoms['atom_id'] = atoms['atom_id'].astype(str)
    if colorby in ["hydrophobicity","binding_prob","residue_type"]:
        atoms[colorby] = atoms[colorby].apply(lambda x: x[0:3])

    # Align to principal Axis
    trans_coords = alignment(atoms, proj_direction)  # get the transformation coordinate
    atoms['x'] = trans_coords[:, 0]
    atoms['y'] = trans_coords[:, 1]
    atoms['z'] = trans_coords[:, 2]

    # convert 3D  to 2D
    atoms["P(x)"] = atoms[['x', 'y', 'z']].apply(lambda coord: projection(coord.x, coord.y, coord.z)[0], axis=1)
    atoms["P(y)"] = atoms[['x', 'y', 'z']].apply(lambda coord: projection(coord.x, coord.y, coord.z)[1], axis=1)

    # setting output image size, labels off, set 120 dpi w x h
    size = 128 if size is None else size
    dpi = 120 if dpi is None else dpi

    figure = plt.figure(figsize=(int(size) / int(dpi), int(size) / int(dpi)), dpi=int(dpi))

    # figsize is in inches, dpi is the resolution of the figure
    ax = plt.subplot(111)  # default is (111)

    ax.axis('off')
    ax.tick_params(axis='both', bottom=False, left=False, right=False,
                   labelleft=False, labeltop=False,
                   labelright=False, labelbottom=False)

    # Compute Voronoi tesselation
    vor = Voronoi(atoms[['P(x)', 'P(y)']])
    regions, vertices = voronoi_finite_polygons_2d(vor)
    polygons = []
    for reg in regions:
        polygon = vertices[reg]
        polygons.append(polygon)
    atoms.loc[:, 'polygons'] = polygons

    # Check alpha
    alpha = float(alpha)

    # Colors color_map
    if colorby in ["charge","center_dist","sasa"]:
        colors = [color_map[_type]["color"] for _type in atoms['atom_id']]
    else:
        colors = [color_map[_type]["color"] for _type in atoms[colorby]]
    atoms["color"] = colors

    for i, row in atoms.iterrows():
        colored_cell = matplotlib.patches.Polygon(row["polygons"],
                                                  facecolor=row['color'],
                                                  edgecolor=row['color'],
                                                  alpha=alpha,
                                                  linewidth=0.2)
        ax.add_patch(colored_cell)

    # Set limits
    ax.set_xlim(vor.min_bound[0], vor.max_bound[0])
    ax.set_ylim(vor.min_bound[1], vor.max_bound[1])

    # Output image saving in any format; default jpg
    bs_out = 'out.jpg' if bs_out is None else bs_out

    # Get image as numpy array
    figure.tight_layout(pad=0)
    img = fig_to_numpy(figure, alpha=alpha)

    if save_fig:
        plt.subplots_adjust(bottom=0, top=1, left=0, right=1)
        plt.savefig(bs_out, frameon=False, pad_inches=False)

    plt.close(fig=figure)

    return atoms, vor, img


def colorby_conv(colorby):
    '''change the colorby parameter to match the dictionary data keys and the column of the atoms dataframe'''

    if colorby in ["atom_type", "charge", "center_dist"]:
        color_by = "residue_type"
    else:
        color_by = colorby
    return color_by


def custom_colormap(color_scale):
    '''takes two hex colors and creates a linear colormap'''

    color_dict = {"red_cyan" : ("#ff0000","#00ffff"), "orange_bluecyan" : ("#ff7f00","#007fff"),
                  "yellow_blue" : ("#ffff00","#0000ff"), "greenyellow_bluemagenta" : ("#7fff00","#7f00ff"),
                  "green_magenta" : ("#00ff00","#ff00ff"), "greencyan_redmagenta" : ("#00ff7f","#ff007f")}

    try:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap1', color_dict[color_scale], N=256)

    except:
        cmap = None     # used for residue_type and atom_type because of predetermined color maps

    return cmap


def normalizer(dataset,colorby):
    '''normalizes dataset using min and max values'''

    valnorm_lst = []
    if colorby not in ["atom_type", "residue_type"]:
        for val in dataset.values():
            val = float(val)
            if max(dataset.values()) == min(dataset.values()):      # used if all values in the set are the same
                valnorm = 0.0
            else:
                valnorm = ((val-min(dataset.values()))/(max(dataset.values())-min(dataset.values())))
            valnorm_lst.append(valnorm)

    return valnorm_lst


def colorgen(colorby,valnorm_lst,cmap,dataset):
    '''creates a new dictionary that contains the color of the given key'''

    # atom type and residue type colors are predetermined
    if colorby in ["atom_type", "residue_type"]:
        color_map = "./cmaps/atom_cmap.csv" if colorby == "atom_type" else "./cmaps/res_cmap.csv"

        # Check for color mapping file, make dict
        with open(color_map, "rt") as color_mapF:
            # Parse color map file
            color_map = np.array(
                [line.replace("\n", "").split(";") for line in color_mapF.readlines() if not line.startswith("#")])
            # Create color dictionary
            color_map = {code: {"color": color, "definition": definition} for code, definition, color in color_map}
            return color_map
    else:
        color_lst = []

        # Apply colormap to the normalized data
        for val in valnorm_lst:
            color = cmap(val)
            color = matplotlib.colors.rgb2hex(color)
            color_lst.append(color)

        # Create color dictionary
        color_map = dict(zip(dataset.keys(), color_lst))
        names = ['code', 'color']
        dtype = dict(names=names)
        hexcolor_array = np.asarray(list(color_map.items()))
        color_map = {code: {"color": color} for code, color in hexcolor_array}
        return color_map


def extract_charge_data(mol):
    '''extracts and formats charge data from mol2 file'''

    # Extracting data from mol2
    pd.options.mode.chained_assignment = None   # Suppress warning
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['atom_id','charge']]       # Only need atom_id and charge data
    atoms.columns = ['atom_id', 'charge']

    # Create dictionary
    charge_list = atoms['charge'].tolist()
    atomid_list = atoms['atom_id'].tolist()
    charge_data = dict(zip(atomid_list, charge_list))

    return charge_data


def extract_centerdistance_data(mol,proj_direction):
    '''extracts and formats center distance from mol2 file after alignment to principal axes'''

    # Extracting data from mol2
    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['atom_id', 'x', 'y', 'z']]
    atoms.columns = ['atom_id', 'x', 'y', 'z']

    # Aligning to principal axes so that origin is the center of pocket
    trans_coords = alignment(atoms, proj_direction)  # get the transformation coordinate
    atoms['x'] = trans_coords[:, 0]
    atoms['y'] = trans_coords[:, 1]
    atoms['z'] = trans_coords[:, 2]

    atomid_list = atoms['atom_id'].tolist()
    coordinate_list = atoms.values.tolist()

    # Calculating the distance to the center of the pocket and creating dictionary
    center_dist_list = []
    for xyz in coordinate_list:
        center_dist = ((xyz[0]) ** 2 + (xyz[1]) ** 2 + (xyz[2]) ** 2) ** .5
        center_dist_list.append(center_dist)
    center_dist_data = dict(zip(atomid_list, center_dist_list))

    return center_dist_data


def extract_sasa_data(mol,pop):
    '''extracts accessible surface area data from .out file generated by POPSlegacy.

        then matches the data in the .out file to the binding site in the mol2 file.

        Used POPSlegacy https://github.com/Fraternalilab/POPSlegacy '''

    # Extracting data from mol2 file
    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['subst_name']]     # only need subst_name for matching. Other data comes from .out file
    atoms.columns = ['residue_type']
    siteresidue_list = atoms['residue_type'].tolist()

    # Extracting sasa data from .out file
    residue_list = []
    qsasa_list = []
    with open(pop) as popsa:  # opening .out file
        for line in popsa:
            line_list = line.split()
            if len(line_list) == 12:  # extracting relevant information
                residue_type = line_list[2] + line_list[4]
                if residue_type in siteresidue_list:
                    qsasa = line_list[7]
                    residue_list.append(residue_type)
                    qsasa_list.append(qsasa)

    # Matching amino acids from .mol2 and .out files and creating dictionary
    qsasa_data = {}
    fullprotein_data = list(zip(residue_list, qsasa_list))
    for i in range(len(fullprotein_data)):
        if fullprotein_data[i][0] in siteresidue_list:
            qsasa_data[i + 1] = float(fullprotein_data[i][1])

    return qsasa_data


def amino_single_to_triple(single):
    '''converts the single letter amino acid abbreviation to the triple letter abbreviation'''

    single_to_triple_dict = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                             'G': 'GLY', 'Q': 'GLN', 'E': 'GLU', 'H': 'HIS', 'I': 'ILE',
                             'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                             'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    for i in single_to_triple_dict.keys():
        if i == single:
            triple = single_to_triple_dict[i]

    return triple


def extract_seq_entropy_data(profile, mol):
    '''extracts sequence entropy data from .profile'''

    # Extracting data from mol2
    pd.options.mode.chained_assignment = None
    mol2 = PandasMol2().read_mol2(mol)
    atoms = mol2.df[['subst_name']]
    atoms.columns = ['residue_type']
    siteresidue_list = atoms['residue_type'].tolist()

    # Opening and formatting lists of the probabilities and residues
    with open(profile) as profile:  # opening .profile file
        ressingle_list = []
        probdata_list = []
        for line in profile:    # extracting relevant information
            line_list = line.split()
            residue_type = line_list[0]
            prob_data = line_list[1:]
            prob_data = list(map(float, prob_data))
            ressingle_list.append(residue_type)
            probdata_list.append(prob_data)

    ressingle_list = ressingle_list[1:]
    probdata_list = probdata_list[1:]

    # Changing single letter amino acid to triple letter with its corresponding number
    count = 0
    restriple_list = []
    for res in ressingle_list:
        newres = res.replace(res, amino_single_to_triple(res))
        count += 1
        restriple_list.append(newres + str(count))

    # Calculating information entropy
    with np.errstate(divide='ignore'):      # suppress warning
        prob_array = np.asarray(probdata_list)
        log_array = np.log2(prob_array)
        log_array[~np.isfinite(log_array)] = 0  # change all infinite values to 0
        entropy_array = log_array * prob_array
        entropydata_array = np.sum(a=entropy_array, axis=1) * -1
        entropydata_list = entropydata_array.tolist()

    # Matching amino acids from .mol2 and .profile files and creating dictionary
    fullprotein_data = dict(zip(restriple_list, entropydata_list))
    seq_entropy_data = {k: float(fullprotein_data[k]) for k in siteresidue_list if k in fullprotein_data}

    return seq_entropy_data

# Hard coded datasets
hydrophobicity_data = {'ALA':1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,
                      'CYS':2.5,'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,
                      'HIS':-3.2,'ILE':4.5,'LEU':3.8,'LYS':-3.9,
                      'MET':1.9,'PHE':2.8,'PRO':-1.6,'SER':-0.8,
                      'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL':4.2}

binding_prob_data = {'ALA':0.701,'ARG':0.916,'ASN':0.811,'ASP':1.015,
                             'CYS':1.650,'GLN':0.669,'GLU':0.956,'GLY':0.788,
                             'HIS':2.286,'ILE':1.006,'LEU':1.045,'LYS':0.468,
                             'MET':1.894,'PHE':1.952,'PRO':0.212,'SER':0.883,
                             'THR':0.730,'TRP':3.084,'TYR':1.672,'VAL':0.884}


def Bionoi(mol, pop, profile, bs_out, size, colorby, dpi, alpha, proj_direction):

    # Dataset and colorscale determined by colorby
    if colorby in ["atom_type","residue_type"]:
        dataset = None
        colorscale = None
    elif colorby == "hydrophobicity":
        dataset = hydrophobicity_data
        colorscale = "red_cyan"
    elif colorby == "charge":
        dataset = extract_charge_data(mol)
        colorscale = "orange_bluecyan"
    elif colorby == "binding_prob":
        dataset = binding_prob_data
        colorscale = "greencyan_redmagenta"
    elif colorby == "center_dist":
        dataset = extract_centerdistance_data(mol,proj_direction)
        colorscale = "yellow_blue"
    elif colorby == "sasa":
        dataset = extract_sasa_data(mol,pop)
        colorscale = "greenyellow_bluemagenta"
    elif colorby == "seq_entropy":
        dataset = extract_seq_entropy_data(profile,mol)
        colorscale = "green_magenta"


    # Run
    cmap = custom_colormap(colorscale)      # create colormap

    valnorm_lst = normalizer(dataset,colorby)   # normalize dataset

    color_map = colorgen(colorby,valnorm_lst,cmap,dataset)      # apply colormap to normalized data

    atoms, vor, img = voronoi_atoms(mol, color_map, colorby, bs_out=bs_out, size=size, dpi=dpi, alpha=alpha,
                                    save_fig=False, proj_direction=proj_direction)      # create voronoi diagram

    return atoms, vor, img