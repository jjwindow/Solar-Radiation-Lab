"""
Jack J. Window
Physics Undergraduate
Imperial College London
03.2020

Standard figure setup for use in academic reports.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

def plotSetup(xlabel, ylabel, xscale = 'log', yscale = 'log'):
    """
    Sets up the standard plot format used for the report. Saves repeating setup
    for every type of plot.

    PARAMS
    ----------------------------------------------------------------------------
    xlabel      :   str, the parameter passed to plt.xlabel.
    ylabel      :   str, the parameter passed to plt.ylabel.
    xscale      :   str, parameter passed to plt.xscale.
    yscale      :   str, parameter passed to plt.yscale.

    RETURNS
    ----------------------------------------------------------------------------
    font_prop   :   FrontProperties object, contains the properties of the 
                    Computer Modern Roman Serif font used in the plot. Returned
                    so that it may be used on the other plot elements outside 
                    this function.
    """
    # set tick width
    mpl.rcParams['xtick.major.size'] = 7
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 1
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 1
    # Get font properties - font chosen is Computer Modern Serif Roman (LaTeX font).
    font_path = 'cmunrm.ttf'
    font_prop = font_manager.FontProperties(fname=font_path, size=18)
    # Set up plot
    plt.grid()
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel(xlabel, fontproperties = font_prop)
    plt.ylabel(ylabel, fontproperties = font_prop)
    # Font has no '-' glyph, so matplotlib serif family used for tick labels.
    plt.xticks(fontfamily = 'serif', fontsize = 14)
    plt.yticks(fontfamily = 'serif', fontsize = 14)
    return font_prop