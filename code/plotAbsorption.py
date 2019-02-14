import numpy as np
from astropy.table import Table, Column
from matplotlib import pyplot as plt

pairs = Table.read('../catalogs/pairs_MgII.fits')
pairs.add_column(Column(np.zeros(len(pairs)), name='WR2796_LIMIT'))

# Compute the limits
# 70 km/s pixel at 2800 Ang is
pix_A = 70.0/3e5*2800.0
nPix = 3.0
pairs['WR2796_LIMIT'] = pix_A*np.sqrt(3.0)/pairs['SN_MGII']*3.0



#pairs = pairs[((pairs['QUALITY_MGII'] == 1) | (pairs['QUALITY_MGII'] == 2)) & (pairs['WR2796_LIMIT'] <= 0.3)]
pairs = pairs[((pairs['QUALITY_MGII'] == 1) | (pairs['QUALITY_MGII'] == 2))]

print(pairs)
print(pairs.colnames)


fig, ax = plt.subplots(1, figsize=(5,5), sharex=True, sharey=True)

index = np.where((pairs['FIT_2796'] == 0) | (pairs['WR_2796']/pairs['WRERR_2796'] <= 3.0))
ax.scatter(pairs[index]['D']*3261.56e-6, pairs[index]['WR2796_LIMIT'], 10, marker='v',
           color='grey', alpha=0.5, label=r'$\rm non{-}' + 'detections\ with\ upper\ limits\ N={}$'.format(len(index[0])))


index = np.where((pairs['FIT_2796'] == 1) & (pairs['WR_2796']/pairs['WRERR_2796'] > 3.0))
ax.scatter(pairs[index]['D']*3261.56e-6, pairs[index]['WR_2796'], 20, 
           color='red', edgecolor='black', alpha=1, label=r'$\rm detections\ with\ measurements\ N={}$'.format(len(index[0])))
ax.errorbar(pairs[index]['D']*3261.56e-6, pairs[index]['WR_2796'],
            pairs[index]['WRERR_2796'],
            fmt='none', color='red', label=None)


ax.legend()
ax.minorticks_on()
ax.set_xlabel(r'$\rm projected\ distance\ [million\ light\ years]$')
ax.set_ylabel(r'$\rm absorption\ strength\ [equivalent\ width\ \AA]$')

ax.set_xlim([0.0, 1])
ax.set_ylim([0.1, 10])
ax.set_yscale('log')

fig.tight_layout()
plt.savefig('../plots/pairs_measurements.pdf')




# make a covering fraction plot.
covering = Table.read('../catalogs/coveringfractions.dat', format='ascii')
limit = 0.5

for cov in covering:
   
   thisPairs = pairs[(pairs[:]['D']*3261.56e-6 >= cov['d_Mly_min']) & (pairs[:]['D']*3261.56e-6 <= cov['d_Mly_max'])]
   index_nondet = np.where((thisPairs['FIT_2796'] == 0) | (thisPairs['WR_2796']/thisPairs['WRERR_2796'] <= 3.0) & (thisPairs['WR2796_LIMIT'] <= 0.5))
   index_det = np.where((thisPairs['FIT_2796'] == 1) & (thisPairs['WR_2796']/thisPairs['WRERR_2796'] > 3.0) & (thisPairs['WR_2796'] > 0.5))
   
   cov['n_det'] = len(index_det[0])*1.0
   cov['n_nondet'] = len(index_nondet[0])*1.0
   cov['n_tot'] = cov['n_det'] + cov['n_nondet']
   
covering['cov'] = covering['n_det']/covering['n_tot']
covering['err'] = np.sqrt(covering['n_det'])/covering['n_tot'] 
   
covering.write('../catalogs/coveringfractions_measured.dat', format='ascii.fixed_width', overwrite=True)



fig, ax = plt.subplots(1, figsize=(5,5), sharex=True, sharey=True)

ax.scatter(covering['d_Mly'], covering['cov'], 30, 
           color='red', edgecolor='black', alpha=1, label=None)
ax.errorbar(covering['d_Mly'], covering['cov'],
            covering['err'],
            fmt='none', color='red', label=None)

ax.legend()
ax.minorticks_on()
ax.set_xlabel(r'$\rm projected\ distance\ [million\ light\ years]$')
ax.set_ylabel(r'$\rm fraction\ of\ quasars\ with\ absorption >{:0.1f} \AA$'.format(limit))

ax.set_xlim([0.0, 1])
ax.set_ylim([0.0, 1.0])

fig.tight_layout()
plt.savefig('../plots/pairs_fraction.pdf')
