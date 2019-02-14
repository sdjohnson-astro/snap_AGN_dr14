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
import lmfit
from scipy.integrate import simps


parser = argparse.ArgumentParser(description='Measure  MgII')
parser.add_argument('-i', metavar='i', type=int, help='index to start with', default=0)
parser.add_argument('-xsize', metavar='xsize', type=int, help='xsize in pixels', default=1000)
parser.add_argument('-ysize', metavar='ysize', type=int, help='ysize in pixels', default=1000)

args = parser.parse_args()

def measureWr_full(wave_obs, flux, error, z):
   
   # Go into the rest-frame
   wave_rest = wave_obs/(1 + z)
   Wr = simps(1-flux, wave_rest)
   dW = np.median(np.roll(wave_rest, -1) - wave_rest)
   WrErr = np.sqrt(simps(error**2, wave_rest))*dW
   
   return Wr, WrErr


def quadratic(wave, a1, a2, a3, a4):
   
   return a1 + a2*wave + a3*wave**2 + a4*wave**3

# Gaussian function
def gaussian(wave, amplitude, centroid, sigma):
   """1-d gaussian: gaussian(x, amp, cen, wid)"""
   return 1 + amplitude*np.exp(-0.5*((wave-centroid)/sigma)**2)


def getspec(filename):
   
   spec = Table.read(filename)
   spec.add_column(Column(np.zeros(len(spec)), 'error'))
   spec.add_column(Column(np.zeros(len(spec)), 'wave'))
   spec['error'] = 1/np.sqrt(spec['ivar'])
   spec['wave'] = 10.0**spec['loglam']
   spec = spec[spec['ivar'] > 0.0]
   
   return spec

class measuremgii:
   """Definiting the checkbackground GUI for running through background galaxy
      spectra and verifying the classification and redshift."""
   
   def mouseMoved_spec(self, pos):
       """Keep track of where the mouse is and update the title with current
          mouse position."""
      
       
       self.mouse_x_spec = self.plot_spec.mapToView(pos).x()
       self.mouse_y_spec = self.plot_spec.mapToView(pos).y()   
       #self.plot_spec.setTitle('{:0.2f}, {:.2E}'.format(self.mouse_x_spec, self.mouse_y_spec))
   
   
   def keypress_spec(self, event):
      
      #print('{} key pressed'.format(event.text()))
      
      # Mark 2796 left 0
      if (event.text() == '1') | (event.text() == '!'):
         
         self.boundary0_2796 = self.mouse_x_spec
         self.pair['BOUNDARY0_2796'] = self.boundary0_2796
         self.draw()
         
      # Mark 2796 left 1
      if (event.text() == '2') | (event.text() == '@'):
         
         self.boundary1_2796 = self.mouse_x_spec
         self.pair['BOUNDARY1_2796'] = self.boundary1_2796
         self.draw()
         
      # Mark 2803 left 0
      if (event.text() == '3') | (event.text() == '#'):
         
         self.boundary0_2803 = self.mouse_x_spec
         self.pair['BOUNDARY0_2803'] = self.boundary0_2803
         self.draw()
         
      # Mark 2803 left 1
      if (event.text() == '4') | (event.text() == '$'):
         
         self.boundary1_2803 = self.mouse_x_spec
         self.pair['BOUNDARY1_2803'] = self.boundary1_2803
         self.draw()
      
      # Mark continuum left 0
      if (event.text() == 'z') | (event.text() == 'Z'):
         
         self.continuum_left0 = self.mouse_x_spec
         self.pair['CONTINUUM_LEFT0'] = self.continuum_left0
         self.draw()
         
      if (event.text() == 'x') | (event.text() == 'X'):
         
         self.continuum_left1 = self.mouse_x_spec
         self.pair['CONTINUUM_LEFT1'] = self.continuum_left1
         self.draw()
         
      if (event.text() == 'c') | (event.text() == 'C'):
         
         self.continuum_right0 = self.mouse_x_spec
         self.pair['CONTINUUM_RIGHT0'] = self.continuum_right0
         self.draw()
      
      if (event.text() == 'v') | (event.text() == 'V'):
         
         self.continuum_right1 = self.mouse_x_spec
         self.pair['CONTINUUM_RIGHT1'] = self.continuum_right1
         self.draw()
      
      # Mark the MgII 2796 redshift
      if (event.text() == 'm') | (event.text() == 'M'):
         
         self.setMgIIredshift()
      
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
         
      if (event.text() == 'l') | (event.text() == 'L'):
         self.setquality('multiple')
         
      if (event.text() == 'r') | (event.text() == 'R'):
         self.setquality('reject')
   
   def setMgIIredshift(self):
      
      wave_MgII_obs = self.mouse_x_spec
      self.z = wave_MgII_obs/2796.35 - 1.0
      self.pair['REDSHIFT_MGII'] = self.z
      self.draw()
   
   def set_spec(self):
      
      self.filename = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_BG'], self.pair['MJD_BG'], self.pair['FIBERID_BG'])
      self.spec = getspec(self.filename)
      
      if self.pair['REDSHIFT_MGII'] == 0.0:
         self.z = self.pair['REDSHIFT_FG']
      else:
         self.z = self.pair['REDSHIFT_MGII']
      
      self.flux_toplot = self.spec['flux']
      self.error_toplot = self.spec['error']
      self.model_toplot = self.spec['model']
      
      wave0 = (2800 - 50)*(1 + self.z)
      wave1 = (2800 + 50)*(1 + self.z)
      self.plot_spec.setXRange(wave0, wave1)
      
      index = np.where((10.0**self.spec['loglam'] > wave0) & (10.0**self.spec['loglam'] < wave1))
      y1 = np.max(self.spec['flux'][index])*1.2
      y0 = -0.1*y1
      self.plot_spec.setYRange(y0, y1)
      
      
      #print('{} {}'.format(self.index, self.filename))
      self.draw()
      
      
   def advance(self, dIndex):
      
      #print(self.index)
      self.index = self.index + dIndex
      if self.index < 0:
         self.index = len(self.pairs)-1
      if self.index > len(self.pairs)-1:
         self.index = 0
         
      self.pair = self.pairs[self.index]
      
      self.set_spec()
      self.smooth_spec()
      #self.plot_spec.autoRange()
      
      
      
      self.pairs.write('../catalogs/pairs_MgII.fits', overwrite=True)
      print('Now looking at MgII index={}'.format(self.pair['INDEX_MGII']))
      
      #self.pairs.write('../catalogs/pairs.txt', format='ascii.fixed_width', overwrite=True)
   
   def draw(self):
      
      self.plot_spec.plot(10.0**self.spec['loglam'], self.flux_toplot,
                          pen=pg.mkPen('w', width=1), clear=True)
      self.plot_spec.plot(10.0**self.spec['loglam'], self.model_toplot,
                          pen=pg.mkPen('r', width=2))
      self.plot_spec.plot(10.0**self.spec['loglam'], self.error_toplot,
                          pen=pg.mkPen('b', width=1))
                          
      
      if self.pair['QUALITY_MGII'] == -1:
         qualitystring = '?'  
      if self.pair['QUALITY_MGII'] == 0:
         qualitystring = 'rejected'  
      if self.pair['QUALITY_MGII'] == 1:
         qualitystring = 'accepted'
      if self.pair['QUALITY_MGII'] == 2:
         qualitystring = 'multiple'
      if self.pair['QUALITY_MGII'] == -999:
         qualitystring = 'uncertain'                  
      
      self.plot_spec.setTitle('MgII index={} of {}: {}-{}-{}   z={:0.4f} quality={}'.format(self.pair['INDEX_MGII'], self.N,
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
      
      
      
      # Draw continuum lines
      if self.pair['CONTINUUM_LEFT0'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['CONTINUUM_LEFT0'],
                             pen=pg.mkPen('g', width=1, style=QtCore.Qt.DotLine)))
                             
      if self.pair['CONTINUUM_LEFT1'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['CONTINUUM_LEFT1'],
                             pen=pg.mkPen('g', width=1, style=QtCore.Qt.DotLine)))
                             
      if self.pair['CONTINUUM_RIGHT0'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['CONTINUUM_RIGHT0'],
                             pen=pg.mkPen('g', width=1, style=QtCore.Qt.DotLine)))
                             
      if self.pair['CONTINUUM_RIGHT1'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['CONTINUUM_RIGHT1'],
                             pen=pg.mkPen('g', width=1, style=QtCore.Qt.DotLine)))
                             
      if self.pair['BOUNDARY0_2796'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['BOUNDARY0_2796'],
                             pen=pg.mkPen('w', width=1, style=QtCore.Qt.DotLine)))
      
      if self.pair['BOUNDARY1_2796'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['BOUNDARY1_2796'],
                             pen=pg.mkPen('w', width=1, style=QtCore.Qt.DotLine)))  
                             
      
      if self.pair['BOUNDARY0_2803'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['BOUNDARY0_2803'],
                             pen=pg.mkPen('w', width=1, style=QtCore.Qt.DotLine)))
      
      if self.pair['BOUNDARY1_2803'] > 0:
         self.plot_spec.addItem(pg.InfiniteLine(self.pair['BOUNDARY1_2803'],
                             pen=pg.mkPen('w', width=1, style=QtCore.Qt.DotLine)))              
      
                             
                             
      # If all the continuum regions are set then fit the continuum and draw
      if (self.pair['CONTINUUM_LEFT0'] > 0) & (self.pair['CONTINUUM_LEFT1'] > 0) & (self.pair['CONTINUUM_RIGHT0'] > 0) & (self.pair['CONTINUUM_RIGHT1'] > 0):
         
         self.fitcontinuum()
      
      if self.pair['CONTINUUM_FIT'] == 1:
         
         continuum = quadratic(self.spec['wave'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'])
         self.plot_spec.plot(self.spec['wave'], continuum,
                             pen=pg.mkPen('g', width=2))
                             
      
      if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['BOUNDARY0_2796'] > 0.0) & (self.pair['BOUNDARY1_2796'] > 0.0):
         
         self.fit2796()
         
      if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['FIT_2796'] == 1):
         
         continuum = quadratic(self.spec['wave'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'])
         model_2796 = continuum*gaussian(self.spec['wave'], self.pair['GAUSS_AMPLITUDE_2796'], self.pair['GAUSS_CENTROID_2796'], self.pair['GAUSS_SIGMA_2796'])
         self.plot_spec.plot(self.spec['wave'], model_2796,
                             pen=pg.mkPen('g', width=2))
      
      if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['BOUNDARY0_2803'] > 0.0) & (self.pair['BOUNDARY1_2803'] > 0.0):
         
         self.fit2803()
                             
      if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['FIT_2803'] == 1):
         
         continuum = quadratic(self.spec['wave'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'])
         model_2803 = continuum*gaussian(self.spec['wave'], self.pair['GAUSS_AMPLITUDE_2803'], self.pair['GAUSS_CENTROID_2803'], self.pair['GAUSS_SIGMA_2803'])
         self.plot_spec.plot(self.spec['wave'], model_2803,
                             pen=pg.mkPen('g', width=2))
         
   def fit2796(self):
      
      index = np.where((self.spec['wave'] >= self.pair['BOUNDARY0_2796']) & (self.spec['wave'] <= self.pair['BOUNDARY1_2796']))
      wave_2796 = self.spec[index]['wave']
      flux_2796 = self.spec[index]['flux']
      error_2796 = self.spec[index]['error']
      continuum_2796 = quadratic(wave_2796, self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'])
      
      flux_normed_2796 = flux_2796/continuum_2796
      error_normed_2796 = error_2796/continuum_2796
      
      model_2796 = lmfit.Model(gaussian)
      parameters = lmfit.Parameters()
      centroid_guess = (self.pair['BOUNDARY0_2796'] + self.pair['BOUNDARY1_2796'])/2.0
      sigma_guess = np.abs(self.pair['BOUNDARY1_2796'] - self.pair['BOUNDARY0_2796'])/5.0
      parameters.add_many(('amplitude',    -0.5,             True,  None, None, None),
                          ('centroid',      centroid_guess,  True,  None, None,  None),
                          ('sigma',         sigma_guess,     True,  None, None,  None))
      
      result2796 = model_2796.fit(flux_normed_2796, wave=wave_2796, params=parameters)
      #print(self.pair['BOUNDARY0_2796'])
      #print(self.pair['BOUNDARY1_2796'])
      #print(result2796.fit_report())
      #print(result2796.success)
      self.pair['FIT_2796'] = 1
      self.pair['GAUSS_CENTROID_2796'] = result2796.best_values['centroid']
      self.pair['GAUSS_CENTROIDERR_2796'] = result2796.params['centroid'].stderr
      self.pair['GAUSS_AMPLITUDE_2796'] = result2796.best_values['amplitude']
      self.pair['GAUSS_AMPLITUDEERR_2796'] = result2796.params['amplitude'].stderr
      self.pair['GAUSS_SIGMA_2796'] = result2796.best_values['sigma']
      self.pair['GAUSS_SIGMAERR_2796'] = result2796.params['sigma'].stderr
      
      #print(np.min(wave_2796))
      #print(np.max(wave_2796))
      
      Wr2796, WrErr2796 = measureWr_full(wave_2796, flux_normed_2796, error_normed_2796, self.pair['REDSHIFT_MGII'])
      #print('Wr = {} +/- {}'.format(Wr2796, WrErr2796))
      #print('z={}'.format(self.pair['REDSHIFT_MGII']))
      self.pair['WR_2796'] = Wr2796
      self.pair['WRERR_2796'] = WrErr2796
      
      
   def fit2803(self):
      
      index = np.where((self.spec['wave'] >= self.pair['BOUNDARY0_2803']) & (self.spec['wave'] <= self.pair['BOUNDARY1_2803']))
      wave_2803 = self.spec[index]['wave']
      flux_2803 = self.spec[index]['flux']
      error_2803 = self.spec[index]['error']
      continuum_2803 = quadratic(wave_2803, self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'])
      
      flux_normed_2803 = flux_2803/continuum_2803
      error_normed_2803 = error_2803/continuum_2803
      
      model_2803 = lmfit.Model(gaussian)
      parameters = lmfit.Parameters()
      centroid_guess = (self.pair['BOUNDARY0_2803'] + self.pair['BOUNDARY1_2803'])/2.0
      sigma_guess = np.abs(self.pair['BOUNDARY1_2803'] - self.pair['BOUNDARY0_2803'])/5.0
      parameters.add_many(('amplitude',    -0.5,             True,  None, None, None),
                          ('centroid',      centroid_guess,  True,  None, None,  None),
                          ('sigma',         sigma_guess,     True,  None, None,  None))
      
      result2803 = model_2803.fit(flux_normed_2803, wave=wave_2803, params=parameters)
      #print(self.pair['BOUNDARY0_2803'])
      #print(self.pair['BOUNDARY1_2803'])
      #print(result2803.fit_report())
      #print(result2803.success)
      self.pair['FIT_2803'] = 1
      self.pair['GAUSS_CENTROID_2803'] = result2803.best_values['centroid']
      self.pair['GAUSS_CENTROIDERR_2803'] = result2803.params['centroid'].stderr
      self.pair['GAUSS_AMPLITUDE_2803'] = result2803.best_values['amplitude']
      self.pair['GAUSS_AMPLITUDEERR_2803'] = result2803.params['amplitude'].stderr
      self.pair['GAUSS_SIGMA_2803'] = result2803.best_values['sigma']
      self.pair['GAUSS_SIGMAERR_2803'] = result2803.params['sigma'].stderr
      
      #print(np.min(wave_2803))
      #print(np.max(wave_2803))
      
      Wr2803, WrErr2803 = measureWr_full(wave_2803, flux_normed_2803, error_normed_2803, self.pair['REDSHIFT_MGII'])
      #print('Wr = {} +/- {}'.format(Wr2803, WrErr2803))
      #print('z={}'.format(self.pair['REDSHIFT_MGII']))
      self.pair['WR_2803'] = Wr2803
      self.pair['WRERR_2803'] = WrErr2803
      
      
   def fitcontinuum(self):
      
      #print('Fitting continuum')
      
      
      index = np.where(((self.spec['wave'] >= self.pair['CONTINUUM_LEFT0']) & (self.spec['wave'] <= self.pair['CONTINUUM_LEFT1'])) | ((self.spec['wave'] >= self.pair['CONTINUUM_RIGHT0']) & (self.spec['wave'] <= self.pair['CONTINUUM_RIGHT1'])))
      
      wave_continuum = self.spec[index]['wave']
      flux_continuum = self.spec[index]['flux']
      error_continuum = self.spec[index]['error']
      
      continuum_model = lmfit.Model(quadratic)
      parameters = lmfit.Parameters()
      parameters.add_many(('a1',     np.median(flux_continuum), True,  None, None, None),
                          ('a2',     0.0,                       True,  None, None,  None),
                          ('a3',     0.0,                       True,  None, None,  None),
                          ('a4',     0.0,                       True,  None, None,  None))
      result_continuum = continuum_model.fit(flux_continuum, wave=wave_continuum, params=parameters, weight=1.0/error_continuum)
      #print(result_continuum.fit_report())
      
      self.pair['CONTINUUM_FIT'] = 1
      self.pair['A1'] = result_continuum.best_values['a1']
      self.pair['A2'] = result_continuum.best_values['a2']
      self.pair['A3'] = result_continuum.best_values['a3']
      self.pair['A4'] = result_continuum.best_values['a4']
                          
   def setquality(self, quality):
   
      if quality == '?':
         self.pair['QUALITY_MGII'] = -999
      
      if quality == 'accept':
         self.pair['QUALITY_MGII'] = 1
         
      if quality == 'multiple':
         self.pair['QUALITY_MGII'] = 2
      
      if quality == 'reject':
         self.pair['QUALITY_MGII'] = 0
         
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
      
      # Continuum left flag
      self.continuum_left0 = 0.0
      self.continuum_left1 = 0.0
      
      self.continuum_right0 = 0.0
      self.continuum_right1 = 0.0
      
      self.boundary0_2796 = 0.0
      self.boundary1_2796 = 0.0
      
      # Right continuum flag
      
      
      # Get the pairs
      self.pairs = Table.read('../catalogs/pairs_MgII.fits')
      self.N = len(self.pairs)
      self.index = index
      
      self.pair = self.pairs[self.index]
      self.filename = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_BG'], self.pair['MJD_BG'], self.pair['FIBERID_BG'])
      self.z = self.pair['REDSHIFT_FG']
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
      self.plot_spec.scene().sigMouseMoved.connect(self.mouseMoved_spec)
      

      # Add plot_spec to the layout
      self.layout.addWidget(self.plot_spec, 0, 0)
      
      self.features = Table.read('MgII.csv', format='ascii')
      
      # Go!
      # Get the spectrum
      self.set_spec()
      
      self.draw()
      self.widget.show()
      self.app.exec_()
      
      

checker = measuremgii(args.i, args.xsize, args.ysize)