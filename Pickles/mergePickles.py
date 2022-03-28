#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:   Merge Plots which were saved as .pickle
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2019-03-06
# Execution:    Import functions / collections (from pyLek.helpers import util)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
import pickle as pkl
import matplotlib.pyplot as plt
from pyLEK.plotters import plot2D
from pyLEK.plotters import plotHelpers
import os
import sys
import numpy as np
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

# 1) Create plots with the plotters using the option savePkl = True

# ------------------------------------------------------------------------------
# Note: You wont need this part of the script in your setup, it just serves
# to create the sample .pickle plots


# 2) Copy all plots you want to combine into a folder
# ------------------------------------------------------------------------------
# Done already

# 3) Copy mergePicklePlots.py into the folder and customize the following to your needs
# ------------------------------------------------------------------------------

def getDataFromPickle(ax):
    line = ax.lines[0]
    x = line.get_xdata()
    y = line.get_ydata()
    legend = ax.get_title()

    return x, y, legend


def parsePlots():
    # Specify the file name
    # fname = str(input("Specify File Name: "))
    fname1 = "plot_as_pickle_01"
    fname2 = "plot_as_pickle_02"
    #fname3 = "plot_as_pickle_03"
    #fname4 = "plot_as_pickle_04"
    #fname5 = "plot_as_pickle_05"

    # Change to current file location
    os.chdir(os.path.dirname(sys.argv[0]))

    # Load the ax-objects
    ax1 = pkl.load(open(fname1 + '.pickle', 'rb'))
    ax2 = pkl.load(open(fname2 + '.pickle', 'rb'))
    #ax3 = pkl.load(open(fname3 + '.pickle', 'rb'))
    #ax4 = pkl.load(open(fname4 + '.pickle', 'rb'))
    #ax5 = pkl.load(open(fname5 + '.pickle', 'rb'))

    # Get data from the ax-objects
    x1, y1, legend1 = getDataFromPickle(ax1)
    x2, y2, legend2 = getDataFromPickle(ax2)
    #x3, y3, legend3 = getDataFromPickle(ax3)
    #x4, y4, legend4 = getDataFromPickle(ax4)
    #x5, y5, legend5 = getDataFromPickle(ax5)

    # Find y lim
    y_min = min(min(y1), min(y2))#, min(y3), min(y4), min(y5))
    y_max = max(max(y1), max(y2))#, max(y3), max(y4), max(y5))
    ylim = [y_min,y_max]

    # Clean up
    plt.close('all')

    # Merged data for plotting
    x = [x1, x2]#, x3, x4, x5]
    y = [y1, y2]#, y3, y4, y5]
    legend = [legend1, legend2]#, legend3, legend4, legend5]

    return x, y, legend, ylim

# 4) Use again the plotters to plot the figures loaded with mergePicklePlots.py
# ------------------------------------------------------------------------------

def replotMerged(x, y, legend, ylim):

    xlabel = 'Schlankheit [h/b]'
    ylabel = 'Maximales Biegemoment [MNm]'
    
    plot2D.plot2D(x, y, title="Maximales Biegemoment",
                  legend=legend, xlabel=xlabel, ylabel=ylabel, ylim=ylim, showPlt=True)

# Execution routine:
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    # 1) Create plots with the plotters using the option savePkl = True

    # 2) Copy all plots you want to combine into a folder
    pass

    # 3) Copy mergePicklePlots.py into the folder and customize the following to your needs
    x, y, legend, ylim = parsePlots()

    # 4) Use again the plotters to plot the figures loaded with mergePicklePlots.py
    replotMerged(x, y, legend, ylim)
