from PyQt5 import QtGui, QtCore  # (the example applies equally well to PySide)
import pyqtgraph as pg
import sys
import os
from astropy.io import fits
from astropy.table import Table, Column
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from astropy.table import Table
import numpy as np
import glob
import argparse

parser = argparse.ArgumentParser(description='FG quasar checking')
parser.add_argument('-i', metavar='i', type=int, help='index to start with', default=0)
parser.add_argument('-xsize', metavar='xsize', type=int, help='xsize in pixels', default=1000)
parser.add_argument('-ysize', metavar='ysize', type=int, help='ysize in pixels', default=1000)

args = parser.parse_args()



def getspec(filename):
   
   spec = Table.read(filename)
   spec.add_column(Column(np.zeros(len(spec)), 'error'))
   spec['error'] = 1/np.sqrt(spec['ivar'])
   spec = spec[spec['ivar'] > 0.0]
   
   return spec

class checkbackground:
   """Definiting the checkbackground GUI for running through background galaxy
      spectra and verifying the classification and redshift."""
   
   
   def keypress_spec(self, event):
      
      #print('{} key pressed'.format(event.text()))
      
      
      if (event.text() == '=') | (event.text() == '+'):
         self.smoothing = self.smoothing + 2
         
         if self.smoothing == 3:
            self.smoothing = 5
               
         self.smooth_spec()
         
      if (event.text() == '-') | (event.text() == '_'):
         self.smoothing = self.smoothing - 2
         
         if self.smoothing < 5:
            self.smoothing = 1
      
         self.smooth_spec()
         
         
      if (event.text() == 'n') | (event.text() == 'N'):
         self.advance(1)
       
      if (event.text() == 'b') | (event.text() == 'B') | (event.text() == 'p') | (event.text() == 'P'):
         self.advance(-1)
         
      
      if (event.text() == '/') | (event.text() == '?')| (event.text() == 'u') | (event.text() == 'U'):
         self.setquality('?')
      
      if (event.text() == 'a') | (event.text() == 'A'):
         self.setquality('accept')
         
      if (event.text() == 'r') | (event.text() == 'R'):
         self.setquality('reject')
   
   def set_spec(self):
      
      self.filename = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_BG'], self.pair['MJD_BG'], self.pair['FIBERID_BG'])
      self.spec = getspec(self.filename)
      self.flux_toplot = self.spec['flux']
      self.error_toplot = self.spec['error']
      self.model_toplot = self.spec['model']
      #print('{} {}'.format(self.index, self.filename))
      self.draw()
      
      
   def advance(self, dIndex):
      
      #print(self.index)
      self.index = self.index + dIndex
      print(self.index)
      if self.index < 0:
         self.index = len(self.pairs)-1
      if self.index > len(self.pairs)-1:
         self.index = 0
         
      self.pair = self.pairs[self.index]
      self.z = self.pair['REDSHIFT_BG']
      self.set_spec()
      self.smooth_spec()
      self.plot_spec.autoRange()
      self.pairs.write('../catalogs/pairs.fits', overwrite=True)
      #self.pairs.write('../catalogs/pairs.txt', format='ascii.fixed_width', overwrite=True)
   
   def draw(self):
      
      self.plot_spec.plot(10.0**self.spec['loglam'], self.flux_toplot,
                          pen=pg.mkPen('w', width=1), clear=True)
      self.plot_spec.plot(10.0**self.spec['loglam'], self.model_toplot,
                          pen=pg.mkPen('r', width=2))
      self.plot_spec.plot(10.0**self.spec['loglam'], self.error_toplot,
                          pen=pg.mkPen('b', width=1))
      
      if self.pair['QUALITY_BG'] == -1:
         qualitystring = '?'  
      if self.pair['QUALITY_BG'] == 0:
         qualitystring = 'rejected'  
      if self.pair['QUALITY_BG'] == 1:
         qualitystring = 'accepted'
      if self.pair['QUALITY_BG'] == -999:
         qualitystring = 'uncertain'                  
      
      self.plot_spec.setTitle('background index={} of {}: {}-{}-{}   z={:0.4f} quality={}'.format(self.index, self.N,
                                                                                        self.pair['PLATE_BG'],
                                                                                        self.pair['MJD_BG'],
                                                                                        self.pair['FIBERID_BG'],
                                                                                        self.z,
                                                                                        qualitystring))
      
      xRange = self.plot_spec.getViewBox().state['viewRange'][0]
      for feature in self.features:

         if (feature['wave']*(1.0 + self.z) >= xRange[0]) & (feature['wave']*(1.0 + self.z) <= xRange[1]):
            self.plot_spec.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z),
                                pen=pg.mkPen('y', width=1, style=QtCore.Qt.DotLine),
                                label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                          
   def setquality(self, quality):
   
      if quality == '?':
         self.pair['QUALITY_BG'] = -999
      
      if quality == 'accept':
         self.pair['QUALITY_BG'] = 1
      
      if quality == 'reject':
         self.pair['QUALITY_BG'] = 0
         
      self.draw()
   
   def smooth_spec(self):
      """Smooth the spectrum using Savitzky-Golay filter."""
      
      if self.smoothing > 1:
         self.flux_toplot = savgol_filter(self.spec['flux'], self.smoothing, 2)
         self.error_toplot = savgol_filter(self.spec['error'], self.smoothing, 2)/np.sqrt(self.smoothing)
         self.model_toplot = savgol_filter(self.spec['model'], self.smoothing, 2)
      if self.smoothing == 1:
         self.flux_toplot = self.spec['flux']
         self.error_toplot = self.spec['error']
         self.model_toplot = self.spec['model']
      
      self.draw()
   
   def __init__(self, index=0, xsize=1000, ysize=1000):
      """Constructor"""
            
      
      # Set initial parameters
      self.z = 0.0
      self.xsize = xsize
      self.ysize = ysize
      self.smoothing = 1
      
      
      # Get the pairs
      self.pairs = Table.read('../catalogs/pairs.fits')
      self.N = len(self.pairs)
      self.index = index
      
      self.pair = self.pairs[self.index]
      self.filename = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_BG'], self.pair['MJD_BG'], self.pair['FIBERID_BG'])
      self.z = self.pair['REDSHIFT_BG']
      # Initialize the gui
      
      self.app = QtGui.QApplication([])       # Always start by initializing Qt (only once per application)
   
      self.widget = QtGui.QWidget()       # Define a top-level widget to hold everything

      # Set the widget size
      self.widget.resize(self.xsize, self.ysize)
      
      # Set the plotting widget
      self.plot_spec = pg.PlotWidget()
      self.plot_spec.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_spec.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_spec.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_spec.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_spec.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_spec.getAxis('top').setStyle(tickLength=-15)
      self.plot_spec.getAxis('left').setStyle(tickLength=-15)
      self.plot_spec.getAxis('right').setStyle(tickLength=-15)
      self.plot_spec.showAxis('right')
      self.plot_spec.showAxis('top')
      self.plot_spec.setLabel('bottom', 'Wavelength [&#8491;]')
      self.plot_spec.setLabel('left', 'Flux')
      
      # Create layout.
      self.layout = QtGui.QGridLayout()
      self.widget.setLayout(self.layout)
      
      # Add keypress handling stuff
      self.plot_spec.keyPressEvent = self.keypress_spec

      # Add plot_spec to the layout
      self.layout.addWidget(self.plot_spec, 0, 0)
      
      self.features = Table.read('quasar.csv', format='ascii')
      
      # Go!
      # Get the spectrum
      self.set_spec()
      
      self.draw()
      self.widget.show()
      self.app.exec_()
      
      

checker = checkbackground(args.i, args.xsize, args.ysize)