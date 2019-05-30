from matplotlib import pyplot as plt
import matplotlib
import numpy as np


red_cyan = ((255,0,0),(0,255,255))
colorarray = np.asarray(red_cyan)
colorarray = colorarray/255
red_cyan = np.array(colorarray).tolist()
a=np.outer(np.arange(0,1,0.001),np.ones(500))
red_cyan_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('red_cyan_cmap', red_cyan, N=256)

lst = (1,2,3.1,4,1.2,-5)


def RGB_to_hex(RGB):
    '''(RGB) -> #FFFFFF'''
    RGB = [int(x) for x in RGB]
    return "#" + "".join(["0{0:x}".format(v).upper() if v < 16 else
                          "{0:x}".format(v).upper() for v in RGB])

def colorgen(data,cmap): #data list -> maps normalized values to colormap

    #normalize data
    valnorm_lst = list()
    rgbcolor_lst = list()
    for val in data:
        val = float(val)
        valnorm = ((val-min(data))/(max(data)-min(data)))
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

colorgen(lst, red_cyan_cmap)

