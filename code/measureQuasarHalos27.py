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

c_kms = 299792.0

# Define some helper functions
def getspec(filename):
   
   spec = Table.read(filename)
   spec.add_column(Column(np.zeros(len(spec)), 'error'))
   spec.add_column(Column(np.zeros(len(spec)), 'wave'))
   spec['error'] = 1/np.sqrt(spec['ivar'])
   spec['wave'] = 10.0**spec['loglam']
   spec = spec[spec['ivar'] > 0.0]
   
   return spec
   
# Measure equivalent width
def measureWr_full(wave_obs, flux, error, z):
   
   # Go into the rest-frame
   wave_rest = wave_obs/(1 + z)
   Wr = simps(1-flux, wave_rest)
   dW = np.median(np.roll(wave_rest, -1) - wave_rest)
   WrErr = np.sqrt(simps(error**2, wave_rest))*dW
   
   return Wr, WrErr
   
# Quadratic function
def polynomial(wave, a0, a1, a2, a3, a4, a5, a6, degree):
   
   if degree == 0:
      print('Degree 0')
      return a0*wave**0
   elif degree == 1:
      print('Degree 1')
      return a0 + a1*wave
   elif degree == 2:
      print('Degree 2')
      return a0 + a1*wave + a2*wave**2
   elif degree == 3:
      print('Degree 3')
      return a0 + a1*wave + a2*wave**2 + a3*wave**3
   elif degree == 4:
      print('Degree 3')
      return a0 + a1*wave + a2*wave**2 + a3*wave**3 + a4*wave**4
   elif degree == 5:
      print('Degree 3')
      return a0 + a1*wave + a2*wave**2 + a3*wave**3 + a4*wave**4  + a5*wave**5
   elif degree == 6:
      print('Degree 3')
      return a0 + a1*wave + a2*wave**2 + a3*wave**3 + a4*wave**4 + a5*wave**5 + a6*wave**6

   
      
# Gaussian function
def gaussian(wave, amplitude, centroid, sigma):
   """1-d gaussian: gaussian(x, amp, cen, wid)"""
   return 1 + amplitude*np.exp(-0.5*((wave-centroid)/sigma)**2)

   


# Define the GUI class for verifying and measuring spectra.
class measurehalos:
   """Defining a GUI for verifying quasar classification / redshifts and measuring
      absorption from the gas around the foreground quasar in the background quasar
      spectrum"""
   
   def mouseMoved_fg(self, pos):
       """Keep track of where the mouse is and update the title with current
          mouse position."""
       self.plot_fg_current = True
       self.plot_bg_current = False
       self.plot_abs_current = False
       
       self.mouse_x_fg = self.plot_fg.mapToView(pos).x()
       self.mouse_y_fg = self.plot_fg.mapToView(pos).y()
       self.setTitle_fg()

       
   def mouseMoved_bg(self, pos):
       """Keep track of where the mouse is and update the title with current
          mouse position."""
       
       self.plot_fg_current = False
       self.plot_bg_current = True
       self.plot_abs_current = False
       
       self.mouse_x_bg = self.plot_bg.mapToView(pos).x()
       self.mouse_y_bg = self.plot_bg.mapToView(pos).y()
       self.setTitle_bg()
       
   def mouseMoved_abs(self, pos):
       """Keep track of where the mouse is and update the title with current
          mouse position."""
       
       self.plot_fg_current = False
       self.plot_bg_current = False
       self.plot_abs_current = True
       
       self.mouse_x_abs = self.plot_abs.mapToView(pos).x()
       self.mouse_y_abs = self.plot_abs.mapToView(pos).y()
       self.setTitle_abs()

   
   def keypress_fg(self,event):
      """Handle all keypress events in the spectrum plot."""
      
      if self.plot_fg_current:
      
         #print('FG keypress: {}'.format(event.text()))
         
         # Quit without saving
         
         if (event.text() == '?') | (event.text() == '/'):
         
            self.help()
            
         # clear
         if (event.text() == '\\') | (event.text() == '|'):
         
            self.clear()
         
         if (event.text() == 'Q'):
            self.widget.close()
            
         # Quit without saving
         if (event.text() == 'q'):
            self.save()
            self.widget.close()
         
         if (event.text() == 'w'):
            self.resetRange_fg()
            
         if (event.text() == 'W'):
            self.plot_fg.autoRange()
            
         if (event.text() == '=') | (event.text() == '+'):
            self.smoothing_fg = self.smoothing_fg + 2
            
            if self.smoothing_fg == 3:
               self.smoothing_fg = 5
            
            #print(self.smoothing_fg)      
            #print(self.smooth_spec_fg()
            
         if (event.text() == '-') | (event.text() == '_'):
            self.smoothing_fg = self.smoothing_fg - 2
            
            if self.smoothing_fg < 5:
               self.smoothing_fg = 1
         
            self.smooth_spec_fg()
            
         if (event.text() == 'n') | (event.text() == 'N'):
            self.advance(1)
          
         if (event.text() == 'b') | (event.text() == 'B') | (event.text() == 'p') | (event.text() == 'P'):
            self.advance(-1)
            
         if (event.text() == 't') | (event.text() == 'T'):
            self.quality_fg(1)
         
         if (event.text() == 'r') | (event.text() == 'R'):
            self.quality_fg(0)
   
   def keypress_bg(self,event):
      """Handle all keypress events in the spectrum plot."""
      
      if self.plot_bg_current:
      
         #print('BG keypress: {}'.format(event.text()))
         
         if (event.text() == '?') | (event.text() == '/'):
         
            self.help()
            
         # clear
         if (event.text() == '\\') | (event.text() == '|'):
         
            self.clear()
         
         # Quit without saving
         if (event.text() == 'Q'):
            self.widget.close()
            
         # Quit without saving
         if (event.text() == 'q'):
            self.save()
            self.widget.close()
         
         if (event.text() == 'w'):
            self.resetRange_bg()
            
         if (event.text() == 'W'):
            self.plot_bg.autoRange()
            
         if (event.text() == '=') | (event.text() == '+'):
            self.smoothing_bg = self.smoothing_bg + 2
            
            if self.smoothing_bg == 3:
               self.smoothing_bg = 5
                  
            self.smooth_spec_bg()
            
         if (event.text() == '-') | (event.text() == '_'):
            self.smoothing_bg = self.smoothing_bg - 2
            
            if self.smoothing_bg < 5:
               self.smoothing_bg = 1
         
            self.smooth_spec_bg()
            
         if (event.text() == 'n') | (event.text() == 'N'):
            self.advance(1)
          
         if (event.text() == 'b') | (event.text() == 'B') | (event.text() == 'p') | (event.text() == 'P'):
            self.advance(-1)
            
         
         if (event.text() == 't') | (event.text() == 'T'):
            self.quality_bg(1)
         
         if (event.text() == 'r') | (event.text() == 'R'):
            self.quality_bg(0)
   
   def clear(self):
      
      print('Clearing')
      
      self.pair['QUALITY_BG'] = -1
      self.pair['QUALITY_FG'] = -1
      self.pair['QUALITY_MGII'] = -1
      self.pair['REDSHIFT_MGII'] = self.pair['REDSHIFT_FG']
      self.pair['CONTINUUM_LEFT0'] = 0.0
      self.pair['CONTINUUM_LEFT1'] = 0.0
      self.pair['CONTINUUM_RIGHT0'] = 0.0
      self.pair['CONTINUUM_RIGHT1'] = 0.0
      self.pair['CONTINUUM_FIT'] = 0.0
      self.pair['DEGREE'] = 1
      self.pair['A0'] = 0.0
      self.pair['A1'] = 0.0
      self.pair['A2'] = 0.0
      self.pair['A3'] = 0.0
      self.pair['FIT_2796'] = 0.0
      self.pair['BOUNDARY0_2796'] = 0.0
      self.pair['BOUNDARY1_2796'] = 0.0
      self.pair['WR_2796'] = 0.0
      self.pair['WRERR_2796'] = 0.0
      self.pair['GAUSS_CENTROID_2796'] = 0.0
      self.pair['GAUSS_CENTROIDERR_2796'] = 0.0
      self.pair['GAUSS_AMPLITUDE_2796'] = 0.0
      self.pair['GAUSS_AMPLITUDEERR_2796'] = 0.0
      self.pair['GAUSS_SIGMA_2796'] = 0.0
      self.pair['GAUSS_SIGMAERR_2796'] = 0.0
      self.pair['FIT_2803'] = 0.0
      self.pair['BOUNDARY0_2803'] = 0.0
      self.pair['BOUNDARY1_2803'] = 0.0
      self.pair['WR_2803'] = 0.0
      self.pair['WRERR_2803'] = 0.0
      self.pair['GAUSS_CENTROID_2803'] = 0.0
      self.pair['GAUSS_CENTROIDERR_2803'] = 0.0
      self.pair['GAUSS_AMPLITUDE_2803'] = 0.0
      self.pair['GAUSS_AMPLITUDEERR_2803'] = 0.0
      self.pair['GAUSS_SIGMA_2803'] = 0.0
      self.pair['GAUSS_SIGMAERR_2803'] = 0.0
         
      self.set_spec()
      self.draw()
   
   def keypress_abs(self,event):
      """Handle all keypress events in the spectrum plot."""
      
      if self.plot_abs_current:
      
         #print('ABS keypress: {}'.format(event.text()))
         
         if (event.text() == '?') | (event.text() == '/'):
         
            self.help()
         
         # Quit without saving
         if (event.text() == 'Q'):
            self.widget.close()
            
         # Quit without saving
         if (event.text() == 'q'):
            self.save()
            self.widget.close()
         
         if (event.text() == 'w'):
            self.resetRange_abs()
            
         if (event.text() == 'W'):
            self.plot_abs.autoRange()
            
         if (event.text() == 'n') | (event.text() == 'N'):
            self.advance(1)
          
         if (event.text() == 'b') | (event.text() == 'B') | (event.text() == 'p') | (event.text() == 'P'):
            self.advance(-1)
            
                  
         if (event.text() == 'r') | (event.text() == 'R'):
            print('Rejecting')
            self.quality_abs(0) # Rejected, bad data
         
         if (event.text() == 'u') | (event.text() == 'U'):
            print('undetected')
            self.quality_abs(1) # undetected
            
         if (event.text() == 'y') | (event.text() == 'Y'):
            print('single-line')
            self.quality_abs(2) # single-line detection
            
         if (event.text() == 't') | (event.text() == 'T'):
            print('high confidence detection')
            self.quality_abs(3) # top level confidence
         
         # Mark the MgII 2796 redshift
         if (event.text() == 'm') | (event.text() == 'M'):
         
            self.setMgIIredshift()
            
         # clear
         if (event.text() == '\\') | (event.text() == '|'):
         
            self.clear()
      
         # Mark continuum left 0
         if (event.text() == 'a') | (event.text() == 'A'):
         
            self.continuum_left0 = self.mouse_x_abs
            self.pair['CONTINUUM_LEFT0'] = self.continuum_left0
            self.draw()
         
         if (event.text() == 's') | (event.text() == 'S'):
         
            self.continuum_left1 = self.mouse_x_abs
            self.pair['CONTINUUM_LEFT1'] = self.continuum_left1
            self.draw()
         
         if (event.text() == 'd') | (event.text() == 'D'):
         
            self.continuum_right0 = self.mouse_x_abs
            self.pair['CONTINUUM_RIGHT0'] = self.continuum_right0
            self.draw()
      
         if (event.text() == 'f') | (event.text() == 'F'):
         
            self.continuum_right1 = self.mouse_x_abs
            self.pair['CONTINUUM_RIGHT1'] = self.continuum_right1
            self.draw()
            
         
         # Mark 2796 left 0
         if (event.text() == 'z') | (event.text() == 'Z'):
         
            self.boundary0_2796 = self.mouse_x_abs
            self.pair['BOUNDARY0_2796'] = self.boundary0_2796
            self.draw()
         
         # Mark 2796 left 1
         if (event.text() == 'x') | (event.text() == 'X'):
         
            self.boundary1_2796 = self.mouse_x_abs
            self.pair['BOUNDARY1_2796'] = self.boundary1_2796
            self.draw()
         
         # Mark 2803 left 0
         if (event.text() == 'c') | (event.text() == 'C'):
         
            self.boundary0_2803 = self.mouse_x_abs
            self.pair['BOUNDARY0_2803'] = self.boundary0_2803
            self.draw()
         
         # Mark 2803 left 1
         if (event.text() == 'v') | (event.text() == 'V'):
         
            self.boundary1_2803 = self.mouse_x_abs
            self.pair['BOUNDARY1_2803'] = self.boundary1_2803
            self.draw()
            
         
         if (event.text() == '0') | (event.text() == ')'):
            self.degree = 0
            self.pair['DEGREE'] = 0
            self.draw()
            
         if (event.text() == '1') | (event.text() == '!'):
            self.degree = 1
            self.pair['DEGREE'] = 1
            self.draw()
            
         if (event.text() == '2') | (event.text() == '@'):
            self.degree = 2
            self.pair['DEGREE'] = 2
            self.draw()
            
         if (event.text() == '3') | (event.text() == '#'):
            self.degree = 3
            self.pair['DEGREE'] = 3
            self.draw()
         
         if (event.text() == '4') | (event.text() == '$'):
            self.degree = 4
            self.pair['DEGREE'] = 4
            self.draw()
         
         if (event.text() == '5') | (event.text() == '%'):
            self.degree = 5
            self.pair['DEGREE'] = 5
            self.draw()
            
         if (event.text() == '6') | (event.text() == '^'):
            self.degree = 6
            self.pair['DEGREE'] = 6
            self.draw()
   
   def help(self):
      
      # Print out help message:
      
      print('Top plot: visually inspect foreground quasar spectrum, fit, and verify redshift')
      print('  Keystrokes:')
      print('     ?: help menu')
      print('     n: next -- save and move to next quasar pair')
      print('     b or p: save and move to previous quasar pair')
      print('     q: quit -- save and quit')
      print('     Q: quit no save -- quit without saving')
      print('     w: show original x and y range')
      print('     W: show full x and y range in data')
      print('     +: increase smoothing')
      print('     -: decrease smoothing')
      print('     t: accept -- accept as good spectrum, classification, or redshift')
      print('     r: reject -- reject as bad spectrum, classification, redshift')
      
      print('')
      print('Middle plot: visually inspect background quasar spectrum, fit, and verify redshift')
      print('  Keystrokes: same as top')
      
      
      print('')
      print('Bottom plot: Estimate quasar continuum and measure absorption from the foreground quasar')
      print('  General Keystrokes: ')
      print('     ?: help menu')
      print('     n: next -- save and move to next quasar pair')
      print('     b or p: save and move to previous quasar pair')
      print('     q: quit -- save and quit')
      print('     Q: quit no save -- quit without saving')
      print('     w: show original x and y range')
      print('     W: show full x and y range in data')
      print('  Measurement Keystrokes: ')
      print('     m: Set nominal redshift for MgII 2796. GUI will show line at expected position of doublets. Use to judge whether this is MgII')      
      print('     a: Set the left-most boundary of the left continuum. Will draw vertical line at location')      
      print('     s: Set the right-most boundary of the left continuum. Will draw vertical line at location')      
      print('     d: Set the left-most boundary of the right continuum. Will draw vertical line at location')      
      print('     f: Set the right-most boundary of the right continuum. Will draw vertical line at location and fit continuum if all regions are set')      
      print('     0: order of continuum fit = 0 -- constant function')
      print('     1: order of continuum fit = 1 -- linear function')
      print('     2: order of continuum fit = 2 -- quadratic function')
      print('     3: order of continuum fit = 3 -- cubic function')
      print('     z: Set the left-most boundary of MgII 2796 feature.')      
      print('     x: Set the right-most boundary of MgII 2796 feature. Will fit a Gaussian if continuum is fit and both 2796 boundaries are set')
      print('     c: Set the left-most boundary of MgII 2803 feature.')      
      print('     v: Set the right-most boundary of MgII 2803 feature. Will fit a Gaussian if continuum is fit and both 2796 boundaries are set')
      print('     t: classify absorber as high confidence (two-line) detection')
      print('     y: classify absorber as potential (single-line) detection')
      print('     u: classify absorber as undetected')
      print('     r: reject due to e.g. bad data')
      print('     \: clear all and start over for this absorber')
      
         
      print('')
      print('Bottom textbox: place to comment on the system if e.g. there is a lot of absorption spread over a large wavelength range. Will save when focus changes to something else.')
      
   
   def setMgIIredshift(self):
      
      wave_MgII_obs = self.mouse_x_abs
      self.z_abs = wave_MgII_obs/2796.35 - 1.0
      self.pair['REDSHIFT_MGII'] = self.z_abs
      self.draw()
      
   def fitcontinuum(self):
      
      print('Fitting continuum')
      
      
      index = np.where(((self.spec_abs['wave'] >= self.pair['CONTINUUM_LEFT0']) & (self.spec_abs['wave'] <= self.pair['CONTINUUM_LEFT1'])) | ((self.spec_abs['wave'] >= self.pair['CONTINUUM_RIGHT0']) & (self.spec_abs['wave'] <= self.pair['CONTINUUM_RIGHT1'])))
      
      wave_continuum = self.spec_abs[index]['wave']
      flux_continuum = self.spec_abs[index]['flux']
      error_continuum = self.spec_abs[index]['error']
      
      continuum_model = lmfit.Model(polynomial)
      parameters = lmfit.Parameters()
      
      if self.degree == 0:
         vary0 = True
         vary1 = False
         vary2 = False
         vary3 = False
         vary4 = False
         vary5 = False
         vary6 = False
      elif self.degree == 1:
         vary0 = True
         vary1 = True
         vary2 = False
         vary3 = False
         vary4 = False
         vary5 = False
         vary6 = False
      elif self.degree == 2:
         vary0 = True
         vary1 = True
         vary2 = True
         vary3 = False
         vary4 = False
         vary5 = False
         vary6 = False
      elif self.degree == 3:
         vary0 = True
         vary1 = True
         vary2 = True
         vary3 = True
         vary4 = False
         vary5 = False
         vary6 = False
      elif self.degree == 4:
         vary0 = True
         vary1 = True
         vary2 = True
         vary3 = True
         vary4 = True
         vary5 = False
         vary6 = False
      elif self.degree == 5:
         vary0 = True
         vary1 = True
         vary2 = True
         vary3 = True
         vary4 = True
         vary5 = True
         vary6 = False
      elif self.degree == 6:
         vary0 = True
         vary1 = True
         vary2 = True
         vary3 = True
         vary4 = True
         vary5 = True
         vary6 = True
         
      parameters.add_many(('a0',     np.median(flux_continuum), vary0,  None, None,  None),
                          ('a1',     0.0,                       vary1,  None, None,  None),
                          ('a2',     0.0,                       vary2,  None, None,  None),
                          ('a3',     0.0,                       vary3,  None, None,  None),
                          ('a4',     0.0,                       vary4,  None, None,  None),
                          ('a5',     0.0,                       vary5,  None, None,  None),
                          ('a6',     0.0,                       vary6,  None, None,  None),
                          ('degree', self.degree,               False, None, None,  None))
      result_continuum = continuum_model.fit(flux_continuum, wave=wave_continuum, params=parameters, weight=1.0/error_continuum)
      
      self.pair['CONTINUUM_FIT'] = 1
      self.pair['A0'] = result_continuum.best_values['a0']
      self.pair['A1'] = result_continuum.best_values['a1']
      self.pair['A2'] = result_continuum.best_values['a2']
      self.pair['A3'] = result_continuum.best_values['a3']
      self.pair['A4'] = result_continuum.best_values['a4']
      self.pair['A5'] = result_continuum.best_values['a5']
      self.pair['A6'] = result_continuum.best_values['a6']
      self.pair['DEGREE'] = self.degree
      
         
   def advance(self, dIndex):
      
      self.save()
      print(self.index)
      self.index = self.index + dIndex
      if self.index < 0:
         self.index = len(self.pairs)-1
      if self.index > len(self.pairs)-1:
         self.index = 0
         
      self.pair = self.pairs[self.index]
      
      
      
      self.set_spec()
      self.smooth_spec_bg()
      self.smooth_spec_fg()
      
   
   def quality_fg(self, qualityKey):
      
      self.pair['QUALITY_FG'] = qualityKey
      self.draw()
      
      
   def quality_bg(self, qualityKey):
      
      self.pair['QUALITY_BG'] = qualityKey
      self.draw()
      
   
   def quality_abs(self, qualityKey):
      
      self.pair['QUALITY_MGII'] = qualityKey
      self.draw()
   
   def smooth_spec_bg(self):
      """Smooth the spectrum using Savitzky-Golay filter."""
      
      if self.smoothing_bg > 1:
         self.flux_bg_toplot = savgol_filter(self.spec_bg['flux'], self.smoothing_bg, 3)
         self.error_bg_toplot = savgol_filter(self.spec_bg['error'],
                                              self.smoothing_bg, 3)/np.sqrt(self.smoothing_bg)
         self.model_bg_toplot = savgol_filter(self.spec_bg['model'], self.smoothing_bg, 3)
      if self.smoothing_bg == 1:
         self.flux_bg_toplot = self.spec_bg['flux']
         self.error_bg_toplot = self.spec_bg['error']
         self.model_bg_toplot = self.spec_bg['model']
         
      self.draw()

   def smooth_spec_fg(self):
      """Smooth the spectrum using Savitzky-Golay filter."""
      
      if self.smoothing_fg > 1:
         self.flux_fg_toplot = savgol_filter(self.spec_fg['flux'], self.smoothing_fg, 3)
         self.error_fg_toplot = savgol_filter(self.spec_fg['error'],
                                              self.smoothing_fg, 3)/np.sqrt(self.smoothing_fg)
         self.model_fg_toplot = savgol_filter(self.spec_fg['model'], self.smoothing_fg, 3)
      if self.smoothing_fg == 1:
         self.flux_fg_toplot = self.spec_fg['flux']
         self.error_fg_toplot = self.spec_fg['error']
         self.model_fg_toplot = self.spec_fg['model']
         
      self.draw()
      
   def resetRange_bg(self):
      
      self.plot_bg.setXRange(np.min( self.wave_bg_toplot),
                             np.max( self.wave_bg_toplot))      
      self.plot_bg.setYRange(np.max( self.model_bg_toplot)*-0.1,
                             np.max( self.model_bg_toplot)*1.3)
   
   def resetRange_fg(self):
      
      self.plot_fg.setXRange(np.min( self.wave_fg_toplot),
                             np.max( self.wave_fg_toplot))      
      self.plot_fg.setYRange(np.max( self.model_fg_toplot)*-0.1,
                             np.max( self.model_fg_toplot)*1.3)
                             
   def resetRange_abs(self):
      
      w_2800 = 2796.35*(1 + self.z_fg)
      dw_plus = 4000.0/c_kms*w_2800
      dw_minus = 3000.0/c_kms*w_2800
      wMin = w_2800 - dw_minus
      wMax = w_2800 + dw_plus
      self.plot_abs.setXRange(wMin, wMax)
      index = np.where((self.wave_abs_toplot > wMin) & (self.wave_abs_toplot < wMax))
      tempMax = np.max(self.model_abs_toplot[index])
      self.plot_abs.setYRange(tempMax*-0.1, tempMax*1.3)
   
   def set_spec(self):
      """Set spectrum"""
      
      # Set the foreground and background redshift
      self.z_fg = self.pair['REDSHIFT_FG']
      self.z_bg = self.pair['REDSHIFT_BG']
      self.z_abs = self.pair['REDSHIFT_MGII']
      self.degree = self.pair['DEGREE']
      self.comment_text.setText(self.pair['COMMENT'])
      
      # Foreground spectrum
      self.filename_fg = '../spectra/fg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_FG'],
                                                                   self.pair['MJD_FG'],
                                                                   self.pair['FIBERID_FG'])
      self.spec_fg = getspec(self.filename_fg)
      
      self.wave_fg_toplot = self.spec_fg['wave']
      self.flux_fg_toplot = self.spec_fg['flux']
      self.error_fg_toplot = self.spec_fg['error']
      self.model_fg_toplot = self.spec_fg['model']
      self.plot_fg.setXRange(np.min( self.wave_fg_toplot),
                             np.max( self.wave_fg_toplot))      
      self.plot_fg.setYRange(np.max( self.model_fg_toplot)*-0.1,
                             np.max( self.model_fg_toplot)*1.3)
                             
                             
                             
      # Background spectrum
      self.filename_bg = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_BG'],
                                                                   self.pair['MJD_BG'],
                                                                   self.pair['FIBERID_BG'])
      self.spec_bg = getspec(self.filename_bg)
      
      self.wave_bg_toplot = self.spec_bg['wave']
      self.flux_bg_toplot = self.spec_bg['flux']
      self.error_bg_toplot = self.spec_bg['error']
      self.model_bg_toplot = self.spec_bg['model']
      self.plot_bg.setXRange(np.min( self.wave_bg_toplot),
                             np.max( self.wave_bg_toplot))      
      self.plot_bg.setYRange(np.max( self.model_bg_toplot)*-0.1,
                             np.max( self.model_bg_toplot)*1.3)
                             
      
      # Absorption spectrum. Independent smoothing, etc.
      spec_abs = getspec(self.filename_bg)
      #spec_abs['wave'] = spec_abs['wave']/(1 + self.z_fg)
      #velocity = (spec_abs['wave'] - 2796.35)/2796.35*c_kms
      
      
      self.spec_abs = spec_abs
      self.wave_abs_toplot = self.spec_abs['wave']
      self.flux_abs_toplot = self.spec_abs['flux']
      self.error_abs_toplot = self.spec_abs['error']
      self.model_abs_toplot = self.spec_abs['model']
      
      w_2800 = 2796.35*(1 + self.z_fg)
      dw_plus = 4000.0/c_kms*w_2800
      dw_minus = 3000.0/c_kms*w_2800
      wMin = w_2800 - dw_minus
      wMax = w_2800 + dw_plus
      self.plot_abs.setXRange(wMin, wMax)
      index = np.where((self.wave_abs_toplot > wMin) & (self.wave_abs_toplot < wMax))
      tempMax = np.max(self.model_abs_toplot[index])
      self.plot_abs.setYRange(tempMax*-0.1, tempMax*1.3)
      
      
   
   
   def save(self):
      print('Saved')    
      print('Saved {}/{}.'.format(self.index+1, self.N))
      self.pairs.write('../catalogs/pairs_MgII_{}.fits'.format(args.u), overwrite=True)
      
      
   
   def setTitle_fg(self):
      
      if self.pair['QUALITY_FG'] == -1:
         quality_fg = '?'
      elif self.pair['QUALITY_FG'] == 0:
         quality_fg = 'rejected'
      elif self.pair['QUALITY_FG'] == 1:
         quality_fg = 'accepted'
         
      name_fg = self.pair['NAME_FG']
      
      
      self.plot_fg.setTitle('{}/{}   foreground    J{} {}  z={:0.4f} x,y={:0.2f}, {:0.2f}'.format(self.index+1, self.N, name_fg, quality_fg, self.z_fg, self.mouse_x_fg, self.mouse_y_fg))

   def setTitle_bg(self):
      
      if self.pair['QUALITY_BG'] == -1:
         quality_bg = '?'
      elif self.pair['QUALITY_BG'] == 0:
         quality_bg = 'rejected'
      elif self.pair['QUALITY_BG'] == 1:
         quality_bg = 'accepted'
         
      name_bg = self.pair['NAME_BG']
      
      
      self.plot_bg.setTitle('{}/{}   background    J{} {}  z={:0.4f}'.format(self.index+1, self.N, name_bg, quality_bg, self.z_bg))
      

   def setTitle_abs(self):
      
      if self.pair['QUALITY_MGII'] == -1:
         quality_abs = '?'
      elif self.pair['QUALITY_MGII'] == 0:
         quality_abs = 'rejected'
      elif self.pair['QUALITY_MGII'] == 1:
         quality_abs = 'undetected'
      elif self.pair['QUALITY_MGII'] == 2:
         quality_abs = 'single-line'
      elif self.pair['QUALITY_MGII'] == 3:
         quality_abs = 'confident'
         
      name_bg = self.pair['NAME_BG']
      
      
      self.plot_abs.setTitle('{}/{}   quality={}   d={:0.1f}   z={:0.4f} Wr(2796) = {:0.2f} +/- {:0.2f}   Wr(2803) = {:0.2f} +/- {:0.2f}'.format(self.index+1, self.N, quality_abs, self.pair['D'], self.z_abs, self.pair['WR_2796'], self.pair['WRERR_2796'], self.pair['WR_2803'], self.pair['WRERR_2803']))
   
   
   
         
   def draw(self):
      
      # Clear plots
      self.plot_fg.clear()
      self.plot_bg.clear()
      self.plot_abs.clear()
      
      # Plot foreground
      self.plot_fg.plot(self.wave_fg_toplot, self.flux_fg_toplot,
                        pen=pg.mkPen('w', width=1), clear=True)
      self.plot_fg.plot(self.wave_fg_toplot, self.model_fg_toplot,
                        pen=pg.mkPen('r', width=2))
      self.plot_fg.plot(self.wave_fg_toplot, self.error_fg_toplot,
                        pen=pg.mkPen('b', width=1))
                        
      xRange = self.plot_fg.getViewBox().state['viewRange'][0]
      
      for feature in self.features_quasar:
         

         if (feature['wave']*(1.0 + self.z_fg) >= xRange[0]) & (feature['wave']*(1.0 + self.z_fg) <= xRange[1]):
            self.plot_fg.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z_fg),
                                pen=pg.mkPen('w', width=2, style=QtCore.Qt.DotLine),
                                label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                        
      
      if self.pair['QUALITY_FG'] == 1:
      
         # Plot background              
         self.plot_bg.plot(self.wave_bg_toplot, self.flux_bg_toplot,
                           pen=pg.mkPen('w', width=1), clear=True)
         self.plot_bg.plot(self.wave_bg_toplot, self.model_bg_toplot,
                           pen=pg.mkPen('r', width=2))
         self.plot_bg.plot(self.wave_bg_toplot, self.error_bg_toplot,
                           pen=pg.mkPen('b', width=1))
         
         for feature in self.features_quasar:
         
            if (feature['wave']*(1.0 + self.z_bg) >= xRange[0]) & (feature['wave']*(1.0 + self.z_bg) <= xRange[1]):
               self.plot_bg.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z_bg),
                                   pen=pg.mkPen('w', width=2, style=QtCore.Qt.DotLine),
                                   label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                   labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                   
         
         if args.q:                 
            for feature in self.features_qsoals:
            
               if (feature['wave']*(1.0 + self.z_abs) >= xRange[0]) & (feature['wave']*(1.0 + self.z_abs) <= xRange[1]):
                  self.plot_bg.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z_abs),
                                       pen=pg.mkPen('y', width=2, style=QtCore.Qt.DotLine),
                                       label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                       labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                   
      
      if (self.pair['QUALITY_FG'] == 1) & (self.pair['QUALITY_BG'] == 1):
                                
         # Plot absorption
         self.plot_abs.plot(self.wave_abs_toplot, self.flux_abs_toplot,
                            pen=None, symbol='o', symbolSize=3, clear=True)
         err = pg.ErrorBarItem(x=self.wave_abs_toplot, y=self.flux_abs_toplot, 
                               top=self.error_abs_toplot, bottom=self.error_abs_toplot)
         self.plot_abs.addItem(err)
         
         
         xRange = self.plot_abs.getViewBox().state['viewRange'][0]
         for feature in self.features:
         
            if (feature['wave']*(1.0 + self.pair['REDSHIFT_MGII']) >= xRange[0]) & (feature['wave']*(1.0 + self.pair['REDSHIFT_MGII']) <= xRange[1]):
               self.plot_abs.addItem(pg.InfiniteLine(feature['wave']*(1 + self.pair['REDSHIFT_MGII']),
                                   pen=pg.mkPen('y', width=2, style=QtCore.Qt.DotLine),
                                   label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                   labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))

         
         # Draw continuum lines
         if self.pair['CONTINUUM_LEFT0'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['CONTINUUM_LEFT0'],
                                pen=pg.mkPen('g', width=2, style=QtCore.Qt.DotLine),
                                label='left0', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                
         if self.pair['CONTINUUM_LEFT1'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['CONTINUUM_LEFT1'],
                                pen=pg.mkPen('g', width=2, style=QtCore.Qt.DotLine),
                                label='left1', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                
         if self.pair['CONTINUUM_RIGHT0'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['CONTINUUM_RIGHT0'],
                                pen=pg.mkPen('g', width=2, style=QtCore.Qt.DotLine),
                                label='right0', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                
         if self.pair['CONTINUUM_RIGHT1'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['CONTINUUM_RIGHT1'],
                                pen=pg.mkPen('g', width=2, style=QtCore.Qt.DotLine),
                                label='right1', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                
         if self.pair['BOUNDARY0_2796'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['BOUNDARY0_2796'],
                                pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine),
                                label='2796 left', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
         
         if self.pair['BOUNDARY1_2796'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['BOUNDARY1_2796'],
                                pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine),
                                label='2796 right', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))  
                                
         
         if self.pair['BOUNDARY0_2803'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['BOUNDARY0_2803'],
                                pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine),
                                label='2803 left', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
         
         if self.pair['BOUNDARY1_2803'] > 0:
            self.plot_abs.addItem(pg.InfiniteLine(self.pair['BOUNDARY1_2803'],
                                pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine),
                                label='2803 right', labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                
         
         # If all the continuum regions are set then fit the continuum and draw
         if (self.pair['CONTINUUM_LEFT0'] > 0) & (self.pair['CONTINUUM_LEFT1'] > 0) & (self.pair['CONTINUUM_RIGHT0'] > 0) & (self.pair['CONTINUUM_RIGHT1'] > 0):
            
            self.fitcontinuum()
         
         if self.pair['CONTINUUM_FIT'] == 1:
            
            print(len(self.spec_abs['wave']))
            continuum = polynomial(self.spec_abs['wave'], self.pair['A0'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'], self.pair['A5'], self.pair['A6'], self.pair['DEGREE'])
            print(continuum.shape)
            self.plot_abs.plot(self.spec_abs['wave'], continuum,
                                pen=pg.mkPen('g', width=2))
         
         
         if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['BOUNDARY0_2796'] > 0.0) & (self.pair['BOUNDARY1_2796'] > 0.0):
            
            self.fit2796()
            
         if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['FIT_2796'] == 1):
            
            continuum = polynomial(self.spec_abs['wave'], self.pair['A0'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'], self.pair['A5'], self.pair['A6'], self.pair['DEGREE'])
            model_2796 = continuum*gaussian(self.spec_abs['wave'], self.pair['GAUSS_AMPLITUDE_2796'], self.pair['GAUSS_CENTROID_2796'], self.pair['GAUSS_SIGMA_2796'])
            self.plot_abs.plot(self.spec_abs['wave'], model_2796,
                                pen=pg.mkPen('r', width=2))
         
         if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['BOUNDARY0_2803'] > 0.0) & (self.pair['BOUNDARY1_2803'] > 0.0):
            
            self.fit2803()
                                
         if (self.pair['CONTINUUM_FIT'] == 1) & (self.pair['FIT_2803'] == 1):
            
            continuum = polynomial(self.spec_abs['wave'], self.pair['A0'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'], self.pair['A5'], self.pair['A6'], self.pair['DEGREE'])
            model_2803 = continuum*gaussian(self.spec_abs['wave'], self.pair['GAUSS_AMPLITUDE_2803'], self.pair['GAUSS_CENTROID_2803'], self.pair['GAUSS_SIGMA_2803'])
            self.plot_abs.plot(self.spec_abs['wave'], model_2803,
                                pen=pg.mkPen('r', width=2))
         
      self.setTitle_fg()
      self.setTitle_bg()
      self.setTitle_abs()
      
      
      
      
      
   def fit2796(self):
      
      index = np.where((self.spec_abs['wave'] >= self.pair['BOUNDARY0_2796']) & (self.spec_abs['wave'] <= self.pair['BOUNDARY1_2796']))
      wave_2796 = self.spec_abs[index]['wave']
      flux_2796 = self.spec_abs[index]['flux']
      error_2796 = self.spec_abs[index]['error']
      continuum_2796 = polynomial(wave_2796, self.pair['A0'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'], self.pair['A5'], self.pair['A6'], self.pair['DEGREE'])
      
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
      
      index = np.where((self.spec_abs['wave'] >= self.pair['BOUNDARY0_2803']) & (self.spec_abs['wave'] <= self.pair['BOUNDARY1_2803']))
      wave_2803 = self.spec_abs[index]['wave']
      flux_2803 = self.spec_abs[index]['flux']
      error_2803 = self.spec_abs[index]['error']
      continuum_2803 = polynomial(wave_2803, self.pair['A0'], self.pair['A1'], self.pair['A2'], self.pair['A3'], self.pair['A4'], self.pair['A5'], self.pair['A6'], self.pair['DEGREE'])
      
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
      
      print(np.min(wave_2803))
      print(np.max(wave_2803))
      
      Wr2803, WrErr2803 = measureWr_full(wave_2803, flux_normed_2803, error_normed_2803, self.pair['REDSHIFT_MGII'])
      print('Wr = {} +/- {}'.format(Wr2803, WrErr2803))
      print('z={}'.format(self.pair['REDSHIFT_MGII']))
      self.pair['WR_2803'] = Wr2803
      self.pair['WRERR_2803'] = WrErr2803
   
   def updateComment(self, event):
      
      print('Updating comment')
      self.pair['COMMENT'] = self.comment_text.text()
   
   def __init__(self, user='', index=0, xsize=1000, ysize=500):
      
      # Set initial parameters
      self.z_fg = 0.0
      self.z_bg = 0.0
      self.z_abs = 0.0
      self.xsize = xsize
      self.ysize = ysize
      self.smoothing_fg = 1
      self.smoothing_bg = 1
      self.mouse_x_bg = 0.0
      self.mouse_y_bg = 0.0
      self.mouse_x_fg = 0.0
      self.mouse_y_fg = 0.0
      self.mouse_x_abs = 0.0
      self.mouse_y_abs = 0.0
      
      self.degree = 0
      
      # Read in list of features
      self.features = Table.read('MgII.csv', format='ascii')
      
      
      
      # Read in the pair table
      print('Reading pair table')
      self.pairs = Table.read('../catalogs/pairs_MgII_{}.fits'.format(args.u))
      self.N = len(self.pairs)
      self.index = index
      
      # Get the current pair
      self.pair = self.pairs[self.index]
      self.filename_bg = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_BG'],
                                                                   self.pair['MJD_BG'],
                                                                   self.pair['FIBERID_BG'])
      self.z_bg = self.pair['REDSHIFT_BG']
      
      self.filename_fg = '../spectra/bg/spec-{}-{}-{}.fits'.format(self.pair['PLATE_FG'],
                                                                   self.pair['MJD_FG'],
                                                                   self.pair['FIBERID_FG'])
      self.z_fg = self.pair['REDSHIFT_FG']
      
      # Get the quasar line list
      self.features_quasar = Table.read('quasar.csv', format='ascii')
      
      self.features_qsoals = Table.read('QSOALS.csv', format='ascii')
      
      
      self.app = QtGui.QApplication([])       # Always start by initializing Qt
      self.widget = QtGui.QWidget()       # Define a top-level widget
      
      # Set the widget size
      self.widget.resize(self.xsize, self.ysize)
      
      # Set the background plotting widget
      print('Create bg quasar plot region')
      self.plot_bg = pg.PlotWidget()
      self.plot_bg.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_bg.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_bg.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_bg.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_bg.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_bg.getAxis('top').setStyle(tickLength=-15)
      self.plot_bg.getAxis('left').setStyle(tickLength=-15)
      self.plot_bg.getAxis('right').setStyle(tickLength=-15)
      self.plot_bg.showAxis('right')
      self.plot_bg.showAxis('top')
      self.plot_bg.setLabel('bottom', 'Wavelength [&#8491;]')
      self.plot_bg.setLabel('left', 'Flux')
      self.setTitle_bg()
      
      
      
      # Set the foreground plotting widget
      print('Create fg quasar plot region')
      self.plot_fg = pg.PlotWidget()
      self.plot_fg.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_fg.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_fg.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_fg.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_fg.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_fg.getAxis('top').setStyle(tickLength=-15)
      self.plot_fg.getAxis('left').setStyle(tickLength=-15)
      self.plot_fg.getAxis('right').setStyle(tickLength=-15)
      self.plot_fg.showAxis('right')
      self.plot_fg.showAxis('top')
      self.plot_fg.setLabel('bottom', 'Wavelength [&#8491;]')
      self.plot_fg.setLabel('left', 'Flux')
      self.setTitle_fg()
      
      
      # Create the absorption spectrum plotting widget.
      print('Create absorption plot region')
      self.plot_abs = pg.PlotWidget()
      self.plot_abs.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_abs.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_abs.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_abs.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_abs.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_abs.getAxis('top').setStyle(tickLength=-15)
      self.plot_abs.getAxis('left').setStyle(tickLength=-15)
      self.plot_abs.getAxis('right').setStyle(tickLength=-15)
      self.plot_abs.showAxis('right')
      self.plot_abs.showAxis('top')
      self.plot_abs.setLabel('bottom', 'Wavelength [&#8491;]')
      self.plot_abs.setLabel('left', 'Flux')
      self.setTitle_abs()
      
      
      # Add comment bar
      self.comment_text = QtGui.QLineEdit('comments here')
      
      # Create layout.
      print('Layout')
      self.layout = QtGui.QGridLayout()
      self.widget.setLayout(self.layout)
      
      # Add mouse listeners
      self.plot_fg.scene().sigMouseMoved.connect(self.mouseMoved_fg)
      self.plot_bg.scene().sigMouseMoved.connect(self.mouseMoved_bg)
      self.plot_abs.scene().sigMouseMoved.connect(self.mouseMoved_abs)
      
      # Add keystroke listeners
      self.plot_bg.keyPressEvent = self.keypress_bg
      self.plot_fg.keyPressEvent = self.keypress_fg
      self.plot_abs.keyPressEvent = self.keypress_abs
      
      # Add text change event to the comment bar
      self.comment_text.focusOutEvent = self.updateComment
      
      # Add plot_bg to the layout
      self.layout.addWidget(self.plot_fg, 0, 0) # background on the bottom
      self.layout.addWidget(self.plot_bg, 1, 0) # foreground on the top
      self.layout.addWidget(self.plot_abs, 2, 0)
      self.layout.addWidget(self.comment_text, 3, 0)
      
      # Get the spectrum
      print('Reading in first spectra')
      self.set_spec()
      
      # Draw
      self.draw()
      
      # Show
      print('Done')
      print('For help, press ?')
      self.widget.show()
      self.app.exec_()
      
print('Starting app')
# Begin by parsing user name, etc
parser = argparse.ArgumentParser(description='Verify foreground/background quasars and measure absorbing gas around the foreground in the background quasar spectrum')
parser.add_argument('-u', metavar='user', type=str, help='user name for file storage', required=True)
parser.add_argument('-q', metavar='qsoals', type=bool, help='plot all qsoals lines', default=False)
parser.add_argument('-i', metavar='i', type=int, help='index to start with', default=0)
parser.add_argument('-xsize', metavar='xsize', type=int, help='xsize in pixels', default=1000)
parser.add_argument('-ysize', metavar='ysize', type=int, help='ysize in pixels', default=1500)
args = parser.parse_args()

# Instantiation
measure = measurehalos(args.u, args.i, args.xsize, args.ysize)