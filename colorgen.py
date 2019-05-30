from matplotlib import pyplot as plt
import matplotlib
import numpy as np


red_cyan = ((255,0,0),(0,255,255))
colorarray = np.asarray(red_cyan)
colorarray = colorarray/255
red_cyan = np.array(colorarray).tolist()
a=np.outer(np.arange(0,1,0.001),np.ones(500))
red_cyan_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('red_cyan_cmap', red_cyan, N=256)


dataset = {'ALA':1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS':2.5,'GLN':-3.5,'GLU':-3.5,
       'GLY':-0.4,'HIS':-3.2,'ILE':4.5,'LEU':3.8,'LYS':-3.9,'MET':1.9,'PHE':2.8,
       'PRO':-1.6,'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL':4.2}

def RGB_to_hex(RGB):
    '''(RGB) -> #FFFFFF'''
    RGB = [int(x) for x in RGB]
    return "#" + "".join(["0{0:x}".format(v).upper() if v < 16 else
                          "{0:x}".format(v).upper() for v in RGB])

def colorgen(data,cmap): #dictionary of data -> maps normalized values to colormap
    #normalize data
    valnorm_lst = list()
    rgbcolor_lst = list()
    for val in data.values():
        val = float(val)
        valnorm = ((val-min(data.values()))/(max(data.values())-min(data.values())))
        valnorm_lst.append(valnorm)
        valnorm_lst.sort()

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
    print(hexcolor_lst)



colorgen(dataset, red_cyan_cmap)