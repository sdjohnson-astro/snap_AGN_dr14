import numpy as np
from astropy.io import fits
from astropy.table import Table, Column


pairs = Table.read('../catalogs/pairs.fits')

wgetfile = open('../spectra/fg/download.sh', 'w')
mvfile = open('../spectra/fg/mv.sh', 'w')


for pair in pairs:
   
   wgetstring = 'wget https://dr14.sdss.org/optical/spectrum/view/data/format\=fits/spec\=lite\?plateid\={}\&mjd\={}\&fiberid\={}\n'.format(pair['PLATE_FG'], pair['MJD_FG'], pair['FIBERID_FG'])
   mvstring = 'mv "spec=lite?plateid={}&mjd={}&fiberid={}" spec-{}-{}-{}.fits\n'.format(pair['PLATE_FG'], pair['MJD_FG'], pair['FIBERID_FG'], pair['PLATE_FG'], pair['MJD_FG'], pair['FIBERID_FG'])
   wgetfile.write(wgetstring)
   mvfile.write(mvstring)

wgetfile.close()
mvfile.close()