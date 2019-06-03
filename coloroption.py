def custom_colormap(colorscale):
    '''takes two hex colors and creates a linear colormap'''
    if colorscale == "red_cyan":
        colorlist = ("#ff0000","#00ffff")
    elif colorscale == "orange_bluecyan":
        colorlist = ("#ff7f00","#007fff")
    elif colorscale == "yellow_blue":
        colorlist = ("#ffff00","#0000ff")
    elif colorscale == "greenyellow_bluemagenta":
        colorlist = ("#7fff00","#7f00ff")
    elif colorscale == "green_magenta":
        colorlist = ("#00ff00","#ff00ff")
    elif colorscale == "greencyan_redmagenta":
        colorlist = ("#00ff7f","#ff007f")
    elif colorscale == "red_orange":
        colorlist = ("#ff0000","#ff7f00")
    elif colorscale == "yellow_yellowgreen":
        colorlist = ("#ffff00","#7fff00")
    elif colorscale == "green_greencyan":
        colorlist = ("#00ff00","#00ff7f")
    elif colorscale == "cyan_bluecyan":
        colorlist = ("#00ffff","#007fff")
    elif colorscale == "blue_bluemagenta":
        colorlist = ("#0000ff","#7f00ff")
    elif colorscale == "magenta_redmagenta":
        colorlist = ("#ff00ff","#ff007f")

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cmap1', colorlist, N=256)

    return cmap