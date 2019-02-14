import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

pairs = Table.read('../catalogs/pairs.fits')
print(pairs)
print(pairs.colnames)

quasars_fg = pairs[pairs['QUALITY_FG'] != -1]
quasars_bg = pairs[pairs['QUALITY_BG'] != -1]

print(quasars_fg)
print('')
print(quasars_bg)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5), sharex=True, sharey=True)

ax1.scatter(quasars_fg['REDSHIFT_FG'], quasars_fg['SNMEDIAN_R_FG'], 7.0, 
           color='black', alpha=0.5, label=r'$\rm accepted\ N={}$'.format((np.sum(quasars_fg['QUALITY_FG'] == 1))))
index = np.where(quasars_fg['QUALITY_FG'] != 1)
           
ax1.scatter(quasars_fg[index]['REDSHIFT_FG'], quasars_fg[index]['SNMEDIAN_R_FG'], 10.0, 
           color='red', alpha=1, label=r'$\rm rejected\ N={}$'.format((np.sum(quasars_fg['QUALITY_FG'] != 1))))

ax1.set_xlim([0.0, 4.0])
ax1.set_ylim([0.1, 100])


ax1.legend()
ax1.set_title(r'$\rm foreground\ quasars$')
ax1.set_xlabel(r'$\rm redshift \propto distance\ from \ our \ Galaxy$')
ax1.set_ylabel(r'$\rm median\ signal\ to\ noise\ ratio$')
ax1.set_yscale('log')
ax1.minorticks_on()


ax2.scatter(quasars_bg['REDSHIFT_BG'], quasars_bg['SNMEDIAN_R_BG'], 7.0, 
           color='black', alpha=0.5, label=r'$\rm accepted\ N={}$'.format((np.sum(quasars_bg['QUALITY_BG'] == 1))))
index = np.where(quasars_bg['QUALITY_BG'] != 1)
           
ax2.scatter(quasars_bg[index]['REDSHIFT_BG'], quasars_bg[index]['SNMEDIAN_R_BG'], 10.0, 
           color='red', alpha=1, label=r'$\rm rejected\ N={}$'.format((np.sum(quasars_bg['QUALITY_BG'] != 1))))

ax2.set_xlim([0.0, 4.0])

ax2.legend(loc=1)
ax2.set_title(r'$\rm background\ quasars$')
ax2.set_xlabel(r'$\rm median\ redshift \propto distance\ from \ our \ Galaxy$')
#ax2.set_ylabel(r'$\rm signal\ to\ noise\ ratio$')
ax2.set_yscale('log')
ax2.minorticks_on()




fig.tight_layout()
plt.savefig('../plots/quasars_quality.pdf')

