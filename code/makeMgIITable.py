import numpy as np
from astropy.table import Table, Column
import argparse



pairs = Table.read('../catalogs/pairs.fits')



pairs = pairs[pairs['WAVE_MGII'] > pairs['WAVE_LYA_FOREST']] # Remove pairs outside of Lya forest
pairs = pairs[(pairs['WAVE_MGII'] > 3800) & (pairs['WAVE_MGII'] < 9200)]

pairs.add_column(Column(np.zeros(len(pairs)), name='SN_MGII'))

index = np.where((pairs['WAVE_MGII'] > 3600) & (pairs['WAVE_MGII'] <= 3800))
pairs['SN_MGII'][index] = pairs['SNMEDIAN_U_BG'][index]

index = np.where((pairs['WAVE_MGII'] > 3800) & (pairs['WAVE_MGII'] <= 5500))
pairs['SN_MGII'][index] = pairs['SNMEDIAN_G_BG'][index]

index = np.where((pairs['WAVE_MGII'] > 5500) & (pairs['WAVE_MGII'] <= 6800))
pairs['SN_MGII'][index] = pairs['SNMEDIAN_R_BG'][index]

index = np.where((pairs['WAVE_MGII'] > 6800) & (pairs['WAVE_MGII'] <= 8200))
pairs['SN_MGII'][index] = pairs['SNMEDIAN_I_BG'][index]

index = np.where((pairs['WAVE_MGII'] > 8200) & (pairs['WAVE_MGII'] <= 10000))
pairs['SN_MGII'][index] = pairs['SNMEDIAN_Z_BG'][index]

print(pairs['WAVE_MGII', 'SNMEDIAN_U_BG', 'SNMEDIAN_G_BG', 'SNMEDIAN_R_BG', 'SNMEDIAN_I_BG', 'SNMEDIAN_Z_BG', 'SN_MGII'])
print(pairs.colnames)



pairs = pairs[pairs['SN_MGII'] >= 10.0]
pairs.sort('D')


pairs.add_column(Column(np.arange(len(pairs)), name='INDEX_MGII'))
pairs.add_column(Column(np.zeros(len(pairs))-1, name='QUALITY_MGII'))
pairs.add_column(Column(np.zeros(len(pairs)), name='REDSHIFT_MGII'))
pairs['REDSHIFT_MGII'] = pairs['REDSHIFT_FG']

# Continuum fitting
pairs.add_column(Column(np.zeros(len(pairs)), name='CONTINUUM_LEFT0'))
pairs.add_column(Column(np.zeros(len(pairs)), name='CONTINUUM_LEFT1'))
pairs.add_column(Column(np.zeros(len(pairs)), name='CONTINUUM_RIGHT0'))
pairs.add_column(Column(np.zeros(len(pairs)), name='CONTINUUM_RIGHT1'))
pairs.add_column(Column(np.zeros(len(pairs)), name='CONTINUUM_FIT'))
pairs.add_column(Column(np.zeros(len(pairs))+1, name='DEGREE'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A0'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A1'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A2'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A3'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A4'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A5'))
pairs.add_column(Column(np.zeros(len(pairs)), name='A6'))


# 2796
pairs.add_column(Column(np.zeros(len(pairs)), name='FIT_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='BOUNDARY0_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='BOUNDARY1_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='WR_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='WRERR_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_CENTROID_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_CENTROIDERR_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_AMPLITUDE_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_AMPLITUDEERR_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_SIGMA_2796'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_SIGMAERR_2796'))

# 2803
pairs.add_column(Column(np.zeros(len(pairs)), name='FIT_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='BOUNDARY0_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='BOUNDARY1_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='WR_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='WRERR_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_CENTROID_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_CENTROIDERR_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_AMPLITUDE_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_AMPLITUDEERR_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_SIGMA_2803'))
pairs.add_column(Column(np.zeros(len(pairs)), name='GAUSS_SIGMAERR_2803'))
pairs.add_column(Column(np.chararray(len(pairs), 100), name='COMMENT'))
pairs['COMMENT'] = ''



pairs.write('../catalogs/pairs_MgII.fits', overwrite=True)