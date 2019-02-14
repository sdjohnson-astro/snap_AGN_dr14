from astropy.table import Table

pairs = Table.read('../catalogs/pairs.fits')
pairs.write('../catalogs/pairs.txt', format='ascii.fixed_width', overwrite=True)