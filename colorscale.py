import numpy as np
import matplotlib
import matplotlib.pyplot as plt

color_dict = {"red_cyan" : ("#ff0000","#00ffff"), "orange_bluecyan" : ("#ff7f00","#007fff"),
                  "yellow_blue" : ("#ffff00","#0000ff"), "greenyellow_bluemagenta" : ("#7fff00","#7f00ff"),
                  "green_magenta" : ("#00ff00","#ff00ff"), "greencyan_redmagenta" : ("#00ff7f","#ff007f")}
color_scale = "red_cyan"
a=np.outer(np.arange(0,1,0.001),np.ones(500))
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap1', color_dict[color_scale], N=256)
plt.imshow(a, cmap=cmap)
plt.yticks(np.arange(0,1000,250), (1.0,0.75,0.5,0.25,0.0,0.0,0))
plt.xticks(np.arange(0,500,500), (0.0,00))
plt.savefig('{}.png'.format(color_scale), dpi=4000, bbox_inches="tight",transparent="True", pad_inches=0)
# plt.show()

