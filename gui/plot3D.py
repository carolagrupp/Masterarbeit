# ------------------------------------------------------------------------------
# Description:  Plotting 3-D Surface on one figure
# Author:       Carola Grupp
# Created:      2022-02-28
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

import matplotlib.pyplot as plt
from matplotlib import cm
import pickle as pkl
import numpy as np
import os
import sys

# ----------------------------------------------------------------------
# Imported functions
# ----------------------------------------------------------------------

import pyLEK.plotters.plotStyle.colorCycler as colorCycler
import pyLEK.plotters.plotStyle.mplStyle as mplStyle
import pyLEK.plotters.plotHelpers as plotHelpers

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------


def plot3D(x, y, z, xlabel=None, ylabel=None, zlabel = None, title=None, legend=None,
           dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
           xlim=[], ylim=[], zlim = [], xscale='linear', yscale='linear', zscale = 'linear',
           style_dict={}, mpl='_2D', colorScheme='Monochrome', variation='color', 
           cmap = cm.coolwarm, linewidth = 0, antialised = False,
           customCycler=None, savePlt=False, savePkl=False, showPlt=False, saveTex=False,
           fig=None, ax=None,annotate=[]):
    """Plotting 2-D Lines (x,y-plot) on one figure in a uniform style
    :param x: list w/ data to plot, with shape [n_row, datapoints]
    :param y: list w/ data to plot, with shape [n_row, datapoints]
    :param xlabel: string w/ labels for x axis
    :param xlabel: string w/ labels for y axis
    :param title: string w/ plot title
    :param legend: list w/ legends [n]
    :param dir_fileName: string w/ Directory / Filename to save to,  
                         must be specified when savePlt is specified
    :param vLines: list w/ floats on where to add vertical line
    :param vTexts: list w/ strings for the vertical lines 
    :param hLines: list w/ floats on where to add horizontal line
    :param hTexts: list w/ strings for the horizontal lines 
    :param xlim: list w/ limits  for x axis [xmin, xmax]
    :param ylim: list w/ limits  for y axis [ymin, ymax]
    :param xscale: string w/ scales acc. to matplotlib
    :param yscale: string w/ scales acc. to matplotlib
    :param style_dict: dict w/ settings to overwrite mplstyle-template
    :param mpl: string w/ name of the mplstyle-sheet
    :param colorScheme: string ('Monochrome', 'UniS')
    :param variation: string ('color', 'linestyle')
    :param customCycler: cycler, instead of colorScheme & variation cycler can be passed 
    :param savePĺt: bool true to save plot 
    :param savePkl: bool true to save as .pickle
    :param saveTex: bool true to save as .pdf_tex
    :param showPlt: bool true show plot in interactive mode
    :param fig: fig object to be overwritten 
    :param ax: ax object to be overwritten 
    :rtype fig: modified fig object
    :rtype ax: modified ax object
    """

    # Find plot styles
    #mplPath = mplStyle.findPlotStyle(mpl)


    # Modify plot styles
    #mplStyle.modifyPlotStyle(style_dict, mplPath)

    # Get the plot styles
    #mplStyle.retrievePlotStyle(style_dict, mplPath)

    # Check font
    #plotHelpers.fontChecker()

    # Prepare Plots
    x = np.transpose(x)
    y = np.transpose(y)

    # An empty figure with one axe
    if fig is None:
        fig, ax = plt.subplots()

    # Setting the title of the axe-object
    if not (title is None):
        ax.set_title(title)

    # Setting the x-axis / y-axis label of the axe-object
    if not (xlabel is None):
        ax.set_xlabel(xlabel)
    if not (ylabel is None):
        ax.set_ylabel(ylabel)
    

    # Setting the x-axis / y-axis limits of the axe-object
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    # Setting the x-axis / y-axis scale of the axe-object
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    #Anmerkungen für Punkte im Plot
    if annotate:
        for a in annotate:
            ax.annotate(a[0],a[1])
        
    # Create color / linestyles
    if customCycler is None:
        customCycler = colorCycler.createCycler(colorScheme, variation)

    # Setting the cycler
    ax.set_prop_cycle(customCycler)

    # 2D - Plot of the axe-object
    #ax.plot(x, y, label='label')

    # Correctly ordering legend entries by replacing labels with entries
    # from the legend list. If this is not done there is not order in
    # the legend see: https://matplotlib.org/users/legend_guide.html
    if not (legend is None):
        handles, labels = ax.get_legend_handles_labels()

        for j, element in enumerate(legend, start=0):
            labels[j] = element

        # Setting the legend
        leg = ax.legend(labels)

        # Set line thickness / transparency of the legend to standart
        for line in leg.get_lines():
            line.set_linewidth(1.5)
            # Do not set for only marker plots
            if "lines.linewidth" in style_dict:
                if style_dict["lines.linewidth"] == 0:
                    line.set_linewidth(0)
            line.set_alpha(1.0)

    # Add vertical line
    if vLines:
        for vLine in vLines:
            # Add vertical line to ax
            ax.axvline(vLine, linestyle="--", linewidth=1.0, marker="None", color='black',
                       zorder=max(len(x), len(y))+1)

    if vTexts:
        for vLine, vText in zip(vLines, vTexts):
            # Add Text to lines
            ax.text(vLine, 0.98 * (ax.get_ylim()[1] - ax.get_ylim()[0]) + ax.get_ylim()[0],
                    vText, rotation=90, rotation_mode='anchor',
                    horizontalalignment='right', verticalalignment='bottom',
                    fontsize='x-small')

    # Add horizontal line
    if hLines:
        for hLine in hLines:
            # Add vertical line to ax
            ax.axhline(hLine, linestyle="--", linewidth=1.0, marker="None", color='black',
                       zorder=max(len(x), len(y))+1)

    if hTexts:
        for hLine, hText in zip(hLines, hTexts):
            # Add Text to lines
            ax.text(0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0]) + ax.get_xlim()[0], hLine,
                    hText, rotation=0, rotation_mode='anchor',
                    horizontalalignment='left', verticalalignment='bottom',
                    fontsize='x-small')
    # Save plot
    if savePlt == True:
        try:
            plt.savefig(dir_fileName)
        except ValueError:
            print("Error saving plot: To save plot specify a file name")

    # Save plot with pickle
    if savePkl == True:
        try:
            pkl.dump(ax, open(dir_fileName + ".pickle", "wb"))
        except TypeError:
            print("Error dumping pickle: To dump pickle specify a file name")

    # Save plot as pdf_tex
    if saveTex == True:
        try:
            plotHelpers.savePdf_tex(fig, dir_fileName)
        except TypeError:
            print("Error saving .pdf_tex: To save pdf_tex specify a file name")

    # Show plot in interactive mode
    if showPlt == True:
        plt.show()

    # Clean up mplstyles
    #mplStyle.cleanPlotStyle(mplPath)

    # Clean up everything
    if fig is None:
        plt.clf()
        plt.close()

    return fig, ax
