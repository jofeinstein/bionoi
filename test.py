#opposite colors
red_cyan = ((255,0,0),(0,255,255))
orange_bluecyan = ((255,127,0),(0,127,255))
yellow_blue = ((255,255,0),(0,0,255))
greenyellow_bluemagenta = ((127,255,0),(127,0,255))
green_magenta = ((0,255,0),(255,0,255))
greencyan_redmagenta = ((0,255,127),(255,0,127))

#neighboring colors
red_orange = ((255,0,0),(255,127,0))
yellow_yellowgreen = ((255,255,0),(127,255,0))
green_greencyan = ((0,255,0),(0,255,127))
cyan_bluecyan = ((0,255,255),(0,127,255))
blue_bluemagenta = ((0,0,255),(127,0,255))
magenta_redmagenta = ((255,0,255),(255,0,127))

colorlist = (red_cyan,orange_bluecyan,yellow_blue,
             greenyellow_bluemagenta,green_magenta,greencyan_redmagenta,
             red_orange,yellow_yellowgreen,green_greencyan,cyan_bluecyan,
             blue_bluemagenta,magenta_redmagenta)



colorarray = np.asarray(red_cyan)
colorarray = colorarray/255
red_cyan = np.array(colorarray).tolist()
a=np.outer(np.arange(0,1,0.001),np.ones(500))
red_cyan_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('red_cyan_cmap', red_cyan, N=256)



















'''plt.imshow(a,cmap=colormap)
    plt.axis('off')
    plt.show()
    plt.savefig('{}.png'.format(color), bbox_inches="tight",transparent="True", pad_inches=0)
    plt.close()'''