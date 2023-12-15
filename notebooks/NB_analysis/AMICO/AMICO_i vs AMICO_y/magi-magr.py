import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, hstack, vstack
from clevar.catalog import ClCatalog
from numpy.ma import masked
import sys
ami_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/Catalog_members.fits'

amy_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_y/Catalog_members.fits'
ami_mb = Table.read(ami_inpath)
amy_mb = Table.read(amy_inpath)
z_min, z_max = .5, .8
ami_mb_cdt = ami_mb[(ami_mb['z']>z_min)*(ami_mb['z']<z_max)]
amy_mb_cdt = amy_mb[(amy_mb['z']>z_min)*(amy_mb['z']<z_max)]
ami_mb_delt = ami_mb_cdt['mag_r'] - ami_mb_cdt['mag_i']
amy_mb_delt = amy_mb_cdt['mag_r'] - amy_mb_cdt['mag_i']
x_bins = np.linspace(16,28,60)
y_bins = np.linspace(-1,2,60)
#print(amy_mb_cdt['mag_i'], amy_mb_delt)
#print(len(amy_mb_cdt['mag_i']), len(amy_mb_delt))
ami_hist = np.histogram2d(ami_mb_cdt['mag_i'], ami_mb_delt, bins = (x_bins,y_bins), weights = ami_mb_cdt['pmem'])
amy_hist = np.histogram2d(amy_mb_cdt['mag_i'], amy_mb_delt, bins = (x_bins,y_bins), weights = amy_mb_cdt['pmem'])
ami_hist = ami_hist[0]
amy_hist = amy_hist[0]
ami_hist[ami_hist==0] = np.nan
amy_hist[amy_hist==0] = np.nan
x, y = np.meshgrid(x_bins, y_bins)
fig, (ax1, ax2) = plt.subplots(1,2, figsize = (14,5))
#print(x)
#print(y)
c = ax1.pcolormesh(x, y, ami_hist.T, cmap='jet', vmin=0, vmax = 30000)
ax1.set_xlabel('mag$_i$', fontsize = 13)
ax1.set_ylabel('mag$_r$ - mag$_i$', fontsize = 13)
ax1.set_title('AMICO$_i$' + str(z_min) + '<z<' + str(z_max), fontsize = 13)
c1 = ax2.pcolormesh(x, y, amy_hist.T, cmap='jet', vmin=0, vmax = 35000)
ax2.set_xlabel('mag$_i$', fontsize = 13)
ax2.set_ylabel('mag$_r$ - mag$_i$', fontsize = 13)
ax2.set_title('AMICO$_y$ ' + str(z_min) + '<z<' + str(z_max), fontsize = 13)
ax1.set_xlim([16,28])
ax1.set_ylim([-.1,2])
ax2.set_xlim([16,28])
ax2.set_ylim([-.1,2])
fig.colorbar(c, ax=ax1, label = 'Nb of members')
fig.colorbar(c1, ax=ax2, label = 'Nb of members')
plt.savefig('/pbs/home/n/namourou/test_jupyter/cluster_challenge/newtypeoffig_58.png')