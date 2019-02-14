import numpy as np
from astropy.io import fits
from astropy.table import Table, Column


pairs = Table.read('../catalogs/pairs.fits')

wgetfile = open('../spectra/bg/download.sh', 'w')
mvfile = open('../spectra/bg/mv.sh', 'w')


for pair in pairs:
   
   wgetstring = 'wget https://dr14.sdss.org/optical/spectrum/view/data/format\=fits/spec\=lite\?plateid\={}\&mjd\={}\&fiberid\={}\n'.format(pair['PLATE_BG'], pair['MJD_BG'], pair['FIBERID_BG'])
   mvstring = 'mv "spec=lite?plateid={}&mjd={}&fiberid={}" spec-{}-{}-{}.fits\n'.format(pair['PLATE_BG'], pair['MJD_BG'], pair['FIBERID_BG'], pair['PLATE_BG'], pair['MJD_BG'], pair['FIBERID_BG'])
   wgetfile.write(wgetstring)
   mvfile.write(mvstring)

wgetfile.close()
mvfile.close()