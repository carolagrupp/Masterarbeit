# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Description:  Framework for custom application
# Author:       ** Add here author's e-mail adress **
# Created:      ** Add here the date of creation **
# Execution:    Import functions / collections (from helpers import util)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
# Literature / Website ressources
# https://stackoverflow.com/questions/43947318/plotting-matplotlib-figure-inside-qwidget-using-qt-designer-form-and-pyqt5/44029435#44029435
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
# Contains all imported modules / functions

from ast import increment_lineno
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import cm
from matplotlib import interactive
from matplotlib.ticker import LinearLocator
import matplotlib

# pyLEK/helpers
import pyLEK.plotters.plot2D as plot2D
import pyLEK.plotters.plotBarChart as plotBarChart
import pyLEK.plotters.plotPieChart as plotPieChart
import pyLEK.plotters.plotSurface as plotSurface

# ------------------------------------------------------------------------------
# Functions / Classes
# ------------------------------------------------------------------------------

# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')

# Matplotlib canvas class to create figure


class MplCanvas(Canvas):
    def __init__(self):
        self.fig, self.ax = plt.subplots(tight_layout=True, projection = '3d')
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(
            self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        Canvas.updateGeometry(self)
        #self.ax = self.fig.add_subplot(projection = '3d')

# Matplotlib widget left

class DreiDMplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        # Inherit from QWidget
        QtWidgets.QWidget.__init__(self, parent)
        self.fig = Figure()
        self.canvas = Canvas(self.fig)                           # Create canvas object
        self.canvas.ax = self.fig.add_subplot(projection = '3d')
        # Create corresponding mpl toolbar
    
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.vbl = QtWidgets.QVBoxLayout()                  # Set box for plotting
        self.vbl.addWidget(self.toolbar)
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)

    def plot3DScatter(self, x, y, z, xlabel=None, ylabel=None, zlabel = None, title=None, legend=None,
               dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
               xlim=[], ylim=[], zlim = [], xscale='linear', yscale='linear', zscale = 'linear',
               style_dict={}, mpl='plotStyle_plot2D', colorScheme='Monochrome', variation='color', 
               customCycler=None, linewidth = 0, antialised = False,
               savePlt=False, savePkl=False, saveTex=False):
        
        # Clear
        self.canvas.ax.clear()
        self.canvas.ax.cla()

        #self.canvas.ax = self.canvas.fig.add_subplot(projection = '3d')
        

        #self.canvas.ax = self.canvas.fig.subplots(projection = '3d')
        #self.canvas.fig, self.canvas.ax = plt.subplots(projection = '3d')

        #self.canvas.ax(projection = '3d')
        #self.canvas.ax(subplot_kw={'projection': '3d'})

        # Plotting: Pass the canvas.fig / ..ax object to the sample plot defined
        # in plot2D.py. The sample plot then calls the plot2D-function, which modifies and returns the
        # MplWidget.canvas.fig / ..ax object

        # self.canvas.fig, self.canvas.ax = plot2D.sample_2(
        #     showPlt=False, fig=self.canvas.fig, ax=self.canvas.ax)

        # In order to set up your own plot, use the following scheme
        # self.canvas.fig, self.canvas.ax = plot2D.plot2D(
        #     x, y, *keyargs, fig=self.canvas.fig, ax=self.canvas.ax)
        #colours = [[0,0,0],[62,68,76],[159,153,152],[185,186,187]]
        #colours = [[0,0.75,1.0],[0,0.5,1.0],[0,0.25,1.0],[0,0,1.0]]
        #colours = ['#8000FF', '#4000FF', '#0000FF', '#0040FF', '#0080FF','#00BFFF']

        surf = self.canvas.ax.scatter(x,y,z, norm=Normalize)   #c=colours,
        #self.fig.colorbar(surf, shrink = 0.5, aspect = 5)
        #cmap = cmap,linewidth = linewidth, antialised = antialised)

        self.canvas.ax.set_zlim(zlim)
        #self.canvas.ax.zaxis.set_major_locator()
        #self.canvas.ax.zaxis.set_major_formatter()

        self.canvas.ax.set_xlabel(xlabel)
        self.canvas.ax.set_ylabel(ylabel)
        self.canvas.ax.set_zlabel(zlabel)

        # Update the canvas
        self.canvas.draw()

        # In order to save the plot, recall the plot function but without passing fig, ax
        # Like this, the plot style will be according to the desired .mplstyle

        # For demonstration plot2D.sample_2() will be used. Here savePlt is set to True
        # The plot will be saved in the current working directory

        # For your own plot:
        # if savePlt == True:
        #     plot2D.plot2D(x, y, *keyargs, dir_fileName='Your_file_name',
        #                   savePlt=savePlt, savePkl=savePkl, saveTex=saveTex, fig=None, ax=None)

        #if (savePlt == True) or (savePkl == True) or (saveTex == True):
            #plot2D.plot2D(
                #x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
                #dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts,  hLines=hLines, hTexts=hTexts,
                #xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale,
                #style_dict=style_dict, mpl=mpl, colorScheme=colorScheme, variation=variation, customCycler=customCycler,
                #savePlt=savePlt, savePkl=savePkl, saveTex=saveTex,
                #fig=None, ax=None)

    def plot3DSurface(self, x, y, z, xlabel=None, ylabel=None, zlabel = None, title=None, legend=None,
               dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
               xlim=[], ylim=[], zlim = [], xscale='linear', yscale='linear', zscale = 'linear',
               style_dict={}, mpl='plotStyle_plot2D', colorScheme='Monochrome', variation='color', 
               cmap = cm.coolwarm, customCycler=None, linewidth = 0, antialised = False,
               savePlt=False, savePkl=False, saveTex=False):
        # Clear
        self.canvas.ax.clear()
        self.canvas.ax.cla()


        #self.canvas.ax = self.canvas.fig.add_subplot(projection = '3d')

        #self.canvas.ax(projection = '3d')
        #self.canvas.ax(subplot_kw={'projection': '3d'})

        x,y = np.meshgrid(x,y)
      
        surf = self.canvas.ax.plot_surface(x,y,z, cmap = cm.coolwarm, vmin = zlim[0], vmax = zlim[1], 
            rstride =1, cstride=1, linewidth=linewidth)
        #cmap = cmap,linewidth = linewidth, antialised = antialised, cmap = plt.cm.YlGnBu_r)

        self.canvas.ax.set_zlim(zlim)
        self.canvas.ax.zaxis.set_major_locator(LinearLocator(10))
        self.canvas.ax.zaxis.set_major_formatter('{x:.00f}')

        self.canvas.ax.set_xlabel(xlabel)
        self.canvas.ax.set_ylabel(ylabel)
        self.canvas.ax.set_zlabel(zlabel)



        cbar = self.fig.colorbar(surf, shrink = 0.25, aspect = 5)
        cbar.remove()
        #self.canvas.fig.colorbar.remove()
        cbar = self.fig.colorbar(surf, shrink = 0.25, aspect = 5)

        #self.canvas.ax.view_init(30,130)

        #for angle in range(0,360):
            #self.canvas.ax.subplot.view_init(30, angle)
            #self.canvas.draw()
            #self.canvas.pause(0.001)

        # Update the canvas
        self.canvas.draw()

        # For your own plot:
        # if savePlt == True:
        #     plot2D.plot2D(x, y, *keyargs, dir_fileName='Your_file_name',
        #                   savePlt=savePlt, savePkl=savePkl, saveTex=saveTex, fig=None, ax=None)
