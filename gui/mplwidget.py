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
        self.fig, self.ax = plt.subplots(tight_layout=True)
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(
            self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        Canvas.updateGeometry(self)

# Matplotlib widget left


class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        # Inherit from QWidget
        QtWidgets.QWidget.__init__(self, parent)
        #if QtWidgets.QWidget.objectName == "MplWidget_470":
            #self.fig = Figure()
            #self.canvas = Canvas(self.fig)                           # Create canvas object
            #self.axes = self.fig.add_subplot(projection = '3d')
        #else:
        self.canvas = MplCanvas()                           # Create canvas object
        # Create corresponding mpl toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.vbl = QtWidgets.QVBoxLayout()                  # Set box for plotting
        self.vbl.addWidget(self.toolbar)
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)

    # Using the pyLEK-plotter functions

    def plot2D(self, x, y, *, xlabel=None, ylabel=None, title=None, legend=None,
               dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
               xlim=[], ylim=[], xscale='linear', yscale='linear',
               style_dict={}, mpl='plotStyle_plot2D', colorScheme='Monochrome', variation='color', customCycler=None,
               savePlt=False, savePkl=False, saveTex=False):
        # Clear
        self.canvas.ax.clear()
        self.canvas.ax.cla()


        # Plotting: Pass the canvas.fig / ..ax object to the sample plot defined
        # in plot2D.py. The sample plot then calls the plot2D-function, which modifies and returns the
        # MplWidget.canvas.fig / ..ax object

        # self.canvas.fig, self.canvas.ax = plot2D.sample_2(
        #     showPlt=False, fig=self.canvas.fig, ax=self.canvas.ax)

        # In order to set up your own plot, use the following scheme
        # self.canvas.fig, self.canvas.ax = plot2D.plot2D(
        #     x, y, *keyargs, fig=self.canvas.fig, ax=self.canvas.ax)
        self.canvas.fig, self.canvas.ax = plot2D.plot2D(
            x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
            vLines=vLines, vTexts=vTexts,  hLines=hLines, hTexts=hTexts,
            xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale,
            style_dict=style_dict, mpl=mpl, colorScheme=colorScheme, variation=variation, customCycler=customCycler,
            fig=self.canvas.fig, ax=self.canvas.ax)

        

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

        if (savePlt == True) or (savePkl == True) or (saveTex == True):
            plot2D.plot2D(
                x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
                dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts,  hLines=hLines, hTexts=hTexts,
                xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale,
                style_dict=style_dict, mpl=mpl, colorScheme=colorScheme, variation=variation, customCycler=customCycler,
                savePlt=savePlt, savePkl=savePkl, saveTex=saveTex,
                fig=None, ax=None)

    


    def plotBarChart(self, y, *, xlabel=None, ylabel=None, title=None, legend=None,
                     xticks=None, xticklabels=None, xticksrotation=None,
                     yticks=None, yticklabels=None, yticksrotation=None,
                     barChart='stacked', annotations=None, annotations_position='above',
                     orientation='vertical',
                     dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
                     xlim=[], ylim=[], xscale='linear', yscale='linear',
                     style_dict={}, mpl='plotStyle_barChart', colorScheme='Monochrome', variation='color', customCycler=None,
                     savePlt=False, savePkl=False, showPlt=False, saveTex=False):
        # Clear
        self.canvas.ax.clear()
        self.canvas.ax.cla()


        # Plotting: Pass the canvas.fig / ..ax object to the sample plot defined
        # in plotBarChart.py. The sample plot then calls the plotBarChart-function, which modifies and returns the
        # MplWidget.canvas.fig / ..ax object

        #self.canvas.fig, self.canvas.ax = plotBarChart.sample_grouped(
        #    showPlt=False, fig=self.canvas.fig, ax=self.canvas.ax)

        # In order to set up your own plot, use the following scheme
        # self.canvas.fig, self.canvas.ax = plotBarChart.plotBarChart(
        #     y, *keyargs, fig=self.canvas.fig, ax=self.gui.canvas.ax)

        self.canvas.fig, self.canvas.ax = plotBarChart.plotBarChart(
            y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend, xticks=xticks,
            xticklabels=xticklabels, xticksrotation=xticksrotation,yticks=yticks, yticklabels=yticklabels, yticksrotation=yticksrotation,
            barChart=barChart, annotations=annotations, annotations_position=annotations_position,
            orientation=orientation, dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts, hLines=hLines, hTexts=hTexts,
            xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale, style_dict=style_dict, mpl=mpl, colorScheme=colorScheme, variation=variation,
            customCycler=customCycler, savePlt=savePlt, savePkl=savePkl, showPlt=showPlt, saveTex=saveTex,
            fig=self.canvas.fig, ax=self.canvas.ax)

        # Update the canvas
        self.canvas.draw()

        # In order to save the plot, recall the plot function but without passing fig, ax
        # Like this, the plot style will be according to the desired .mplstyle

        # For your own plot:
        if (savePlt == True) or (savePkl == True) or (saveTex == True):
            plotBarChart.plotBarChart(
                y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend, xticks=xticks,
                xticklabels=xticklabels, xticksrotation=xticksrotation,yticks=yticks, yticklabels=yticklabels, yticksrotation=yticksrotation,
                barChart=barChart, annotations=annotations, annotations_position=annotations_position,
                orientation=orientation, dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts, hLines=hLines, hTexts=hTexts,
                xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale, style_dict=style_dict, mpl=mpl, colorScheme=colorScheme, variation=variation,
                customCycler=customCycler, savePlt=savePlt, savePkl=savePkl, showPlt=showPlt, saveTex=saveTex,
                fig=None, ax=None)
        # if savePlt == True:
        #     plotBarChart.plotBarChart(y, *keyargs, dir_fileName='Your_file_name', savePlt=True, fig=None, ax=None)

    def plotPieChart(self, y, *, xlabel=None, ylabel=None, title=None, legend=None,
                     xticks=None, xticklabels=None, xticksrotation=None,
                     yticks=None, yticklabels=None, yticksrotation=None,
                     barChart='stacked', annotations=None, annotations_position='above',
                     orientation='vertical',
                     dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
                     xlim=[], ylim=[], xscale='linear', yscale='linear',
                     style_dict={}, mpl='default', colorScheme='Monochrome', variation='color', customCycler=None,
                     savePlt=False, savePkl=False, showPlt=False, saveTex=False):
        # Clear
        self.canvas.ax.clear()
        self.canvas.ax.cla()


        # Plotting: Pass the canvas.fig / ..ax object to the sample plot defined
        # in plotBarChart.py. The sample plot then calls the plotBarChart-function, which modifies and returns the
        # MplWidget.canvas.fig / ..ax object

        self.canvas.fig, self.canvas.ax = plotPieChart.sample1(
            showPlt=False, fig=self.canvas.fig, ax=self.canvas.ax)

        # In order to set up your own plot, use the following scheme
        # self.canvas.fig, self.canvas.ax = plotPieChart.plotPieChart(
        #     y, *keyargs, fig=self.canvas.fig, ax=self.gui.canvas.ax)

        # Update the canvas
        self.canvas.draw()

        # In order to save the plot, recall the plot function but without passing fig, ax
        # Like this, the plot style will be according to the desired .mplstyle

        # For your own plot:
        # if savePlt == True:
        #     plotPieChart.plotPieChart(y, *keyargs, dir_fileName='Your_file_name', savePlt=True, fig=None, ax=None)

    def plotSurface(self, x, y, z, *, xlabel=None, ylabel=None, zlabel=None, title=None, legend=None,
               dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
               xlim=[], ylim=[], zlim=[], xscale='linear', yscale='linear', zscale='linear',
               style_dict={}, mpl='_3D', cmap=cm.coolwarm, rstride=1, cstride=1, linewidth=0, colorScheme='Monochrome', variation='color', customCycler=None,
               savePlt=False, savePkl=False, saveTex=False):
        # Clear
        self.canvas.ax.clear()
        self.canvas.ax.cla()

        if len(zlim) !=0:
            vmin= zlim[0]
            vmax = zlim[1]
        else:
            vmin = None
            vmax = None


        # Plotting: Pass the canvas.fig / ..ax object to the sample plot defined
        # in plot2D.py. The sample plot then calls the plot2D-function, which modifies and returns the
        # MplWidget.canvas.fig / ..ax object

        # self.canvas.fig, self.canvas.ax = plot2D.sample_2(
        #     showPlt=False, fig=self.canvas.fig, ax=self.canvas.ax)

        # In order to set up your own plot, use the following scheme
        # self.canvas.fig, self.canvas.ax = plot2D.plot2D(
        #     x, y, *keyargs, fig=self.canvas.fig, ax=self.canvas.ax)
        self.canvas.fig, self.canvas.ax = plotSurface.plotSurface(x, y, z,
            xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, title=title, legend=legend,
           dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts,  hLines=hLines, hTexts=hTexts,
           xlim=xlim, ylim=ylim, zlim=zlim,  xscale=xscale, yscale=yscale, zscale=zscale,
           style_dict={}, mpl=mpl, cmap=cmap, vmin=vmin, vmax=vmax, rstride=rstride, cstride=cstride, linewidth=linewidth,
           colorScheme='Monochrome', variation='color', customCycler=None,
           savePlt=False, savePkl=False, showPlt=False, saveTex=False,
           fig=self.canvas.fig, ax=self.canvas.ax)

        

        

        self.canvas.ax.view_init(30,40)

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

        if (savePlt == True) or (savePkl == True) or (saveTex == True):
            plotSurface.plotSurface(
                x, y, z, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, title=title, legend=legend,
                dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts,  hLines=hLines, hTexts=hTexts,
                xlim=xlim, ylim=ylim, zlim=zlim,  xscale=xscale, yscale=yscale, zscale=zscale,
                style_dict={}, mpl=mpl, cmap=cmap, vmin=vmin, vmax=vmax, rstride=rstride, cstride=cstride, linewidth=linewidth,
                colorScheme='Monochrome', variation='color', customCycler=None,
                savePlt=savePlt, savePkl=savePkl, saveTex=saveTex,
                fig=None, ax=None)


     #def plot3DScatter(self):#, x, y, z, xlabel=None, ylabel=None, zlabel = None, title=None, legend=None,
    #            dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
    #            xlim=[], ylim=[], zlim = [], xscale='linear', yscale='linear', zscale = 'linear',
    #            style_dict={}, mpl='plotStyle_plot2D', colorScheme='Monochrome', variation='color', 
    #            customCycler=None, linewidth = 0, antialised = False,
    #            savePlt=False, savePkl=False, saveTex=False):
        
    #     # Clear
    #     self.canvas.ax.clear()
    #     self.canvas.ax.cla()

    #     self.canvas.ax.sub = self.canvas.fig.add_subplot(projection = '3d')
        

    #     #self.canvas.ax = self.canvas.fig.subplots(projection = '3d')
    #     #self.canvas.fig, self.canvas.ax = plt.subplots(projection = '3d')

    #     #self.canvas.ax(projection = '3d')
    #     #self.canvas.ax(subplot_kw={'projection': '3d'})

    #     # Plotting: Pass the canvas.fig / ..ax object to the sample plot defined
    #     # in plot2D.py. The sample plot then calls the plot2D-function, which modifies and returns the
    #     # MplWidget.canvas.fig / ..ax object

    #     # self.canvas.fig, self.canvas.ax = plot2D.sample_2(
    #     #     showPlt=False, fig=self.canvas.fig, ax=self.canvas.ax)

    #     # In order to set up your own plot, use the following scheme
    #     # self.canvas.fig, self.canvas.ax = plot2D.plot2D(
    #     #     x, y, *keyargs, fig=self.canvas.fig, ax=self.canvas.ax)
    #     #colours = [[0,0,0],[62,68,76],[159,153,152],[185,186,187]]
    #     #colours = [[0,0.75,1.0],[0,0.5,1.0],[0,0.25,1.0],[0,0,1.0]]
    #     #colours = ['#8000FF', '#4000FF', '#0000FF', '#0040FF', '#0080FF','#00BFFF']

    #     surf = self.canvas.ax.sub.scatter(x,y,z, norm=Normalize)   #c=colours,
    #     #self.fig.colorbar(surf, shrink = 0.5, aspect = 5)
    #     #cmap = cmap,linewidth = linewidth, antialised = antialised)

    #     self.canvas.ax.sub.set_zlim(zlim)
    #     #self.canvas.ax.zaxis.set_major_locator()
    #     #self.canvas.ax.zaxis.set_major_formatter()

    #     self.canvas.ax.sub.set_xlabel(xlabel)
    #     self.canvas.ax.sub.set_ylabel(ylabel)
    #     self.canvas.ax.sub.set_zlabel(zlabel)

    #     # Update the canvas
    #     self.canvas.draw()

    #     # In order to save the plot, recall the plot function but without passing fig, ax
    #     # Like this, the plot style will be according to the desired .mplstyle

    #     # For demonstration plot2D.sample_2() will be used. Here savePlt is set to True
    #     # The plot will be saved in the current working directory

    #     # For your own plot:
    #     # if savePlt == True:
    #     #     plot2D.plot2D(x, y, *keyargs, dir_fileName='Your_file_name',
    #     #                   savePlt=savePlt, savePkl=savePkl, saveTex=saveTex, fig=None, ax=None)

    #     #if (savePlt == True) or (savePkl == True) or (saveTex == True):
    #         #plot2D.plot2D(
    #             #x, y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend,
    #             #dir_fileName=dir_fileName, vLines=vLines, vTexts=vTexts,  hLines=hLines, hTexts=hTexts,
    #             #xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale,
    #             #style_dict=style_dict, mpl=mpl, colorScheme=colorScheme, variation=variation, customCycler=customCycler,
    #             #savePlt=savePlt, savePkl=savePkl, saveTex=saveTex,
    #             #fig=None, ax=None)

    # def plot3DSurface(self, x, y, z, xlabel=None, ylabel=None, zlabel = None, title=None, legend=None,
    #            dir_fileName=None, vLines=None, vTexts=None,  hLines=None, hTexts=None,
    #            xlim=[], ylim=[], zlim = [], xscale='linear', yscale='linear', zscale = 'linear',
    #            style_dict={}, mpl='plotStyle_plot2D', colorScheme='Monochrome', variation='color', 
    #            cmap = cm.coolwarm, customCycler=None, linewidth = 0, antialised = False,
    #            savePlt=False, savePkl=False, saveTex=False):
    #     # Clear
    #     self.canvas.ax.clear()
    #     self.canvas.ax.cla()

    #     self.canvas.ax.sub.clear()
    #     self.canvas.ax.sub.cla()


    #     self.canvas.ax.sub = self.canvas.fig.add_subplot(projection = '3d')

    #     #self.canvas.ax(projection = '3d')
    #     #self.canvas.ax(subplot_kw={'projection': '3d'})

    #     x,y = np.meshgrid(x,y)
      
    #     surf = self.canvas.ax.sub.plot_surface(x,y,z, cmap = cm.coolwarm, vmin = zlim[0], vmax = zlim[1], 
    #         rstride =1, cstride=1, linewidth=linewidth)
    #     #cmap = cmap,linewidth = linewidth, antialised = antialised, cmap = plt.cm.YlGnBu_r)

    #     self.canvas.ax.sub.set_zlim(zlim)
    #     self.canvas.ax.sub.zaxis.set_major_locator(LinearLocator(10))
    #     self.canvas.ax.sub.zaxis.set_major_formatter('{x:.00f}')

    #     self.canvas.fig.colorbar(surf, shrink = 0.25, aspect = 5)

    #     self.canvas.ax.sub.view_init(30,40)

    #     #for angle in range(0,360):
    #         #self.canvas.ax.subplot.view_init(30, angle)
    #         #self.canvas.draw()
    #         #self.canvas.pause(0.001)

    #     # Update the canvas
    #     self.canvas.draw()

    #     # For your own plot:
    #     # if savePlt == True:
    #     #     plot2D.plot2D(x, y, *keyargs, dir_fileName='Your_file_name',
    #     #                   savePlt=savePlt, savePkl=savePkl, saveTex=saveTex, fig=None, ax=None)

    


class MplWidget3D(QtWidgets.QWidget):
    def __init__(self, parent=None):
        # Inherit from QWidget
        QtWidgets.QWidget.__init__(self, parent)
        self.fig = Figure()
        self.canvas = Canvas(self.fig)                           # Create canvas object
        self.axes = self.fig.add_subplot(projection = '3d')
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
        self.axes.clear()
        self.axes.cla()

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
        colours = [[0,0.75,1.0],[0,0.5,1.0],[0,0.25,1.0],[0,0,1.0]]
        #colours = ['#8000FF', '#4000FF', '#0000FF', '#0040FF', '#0080FF','#00BFFF']

        surf = self.axes.scatter(x,y,z, c=colours, norm=Normalize)
        #self.fig.colorbar(surf, shrink = 0.5, aspect = 5)
        #cmap = cmap,linewidth = linewidth, antialised = antialised)

        self.axes.set_zlim(zlim)
        #self.canvas.ax.zaxis.set_major_locator()
        #self.canvas.ax.zaxis.set_major_formatter()

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.axes.set_zlabel(zlabel)

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
        self.axes.clear()
        self.axes.cla()

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


        surf = self.axes.plot_surface(x,y,z, cmap = cmap, linewidth=linewidth, antialised = antialised)
        #cmap = cmap,linewidth = linewidth, antialised = antialised)

        self.axes.set_zlim(zlim)
        #self.canvas.ax.zaxis.set_major_locator()
        #self.canvas.ax.zaxis.set_major_formatter()

        self.fig.colorbar(surf, shrink = 0.5, aspect = 5)

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
