from matplotlib import pyplot as plt
import matplotlib
import numpy as np

#opposite colors
red_cyan = ((255,0,0),(0,255,255))
orange_bluecyan = ((255,127,0),(0,127,255))
yellow_blue = ((255,255,0),(0,0,255))
greenyellow_bluemagenta = ((127,255,0),(127,0,255))
green_magenta = ((0,255,0),(255,0,255))
greencyan_redmagenta = ((0,255,127),(255,0,127))


colorlist = (red_cyan,orange_bluecyan,yellow_blue,
             greenyellow_bluemagenta,green_magenta,greencyan_redmagenta)


for color in colorlist:
    colorarray = np.asarray(color)
    colorarray = colorarray/255
    color = np.array(colorarray).tolist()
    a=np.outer(np.arange(0,1,0.001),np.ones(300))
    colormap = matplotlib.colors.LinearSegmentedColormap.from_list('red_cyan_cmap', color, N=256)
    plt.imshow(a,cmap=colormap)
    plt.axis(False)
    plt.show()
    plt.savefig('{}.png'.format(color), bbox_inches="tight",transparent="True", pad_inches=0)
    plt.close()