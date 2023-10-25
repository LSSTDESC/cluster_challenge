#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, hstack, vstack
from clevar.catalog import ClCatalog
from astropy.io import ascii

inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/'
c1 = Table.read(inpath + 'Catalog.fits')
c1_mb = Table.read(inpath + 'Catalog_members.fits')

cosmo_values = ascii.read('/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/magstar_files/istar.asc')
cosmo_values.rename_column('col1', 'z')
cosmo_values.rename_column('col2', 'mag_i')

zeros1 = np.zeros(13)
c1_mb_c = Table(zeros1, names = ('id', 'id_cluster', 'ra', 'dec', 'z', 'pmem', 'lambstar', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y', 'ind_cl'))

c1_mb['id'] = c1_mb['id'].astype(int)
c1_mb['id_cluster'] = c1_mb['id_cluster'].astype(int)

for i in range(len(cosmo_values['z'])-1):
    z_min = cosmo_values['z'][i]
    z_max = cosmo_values['z'][i+1]
    mag_i_max = cosmo_values['mag_i'][i] + 1.5
    cdt1 = c1_mb[(c1_mb['z']<z_max)*(c1_mb['z']>z_min)*(c1_mb['mag_i']<mag_i_max)]
    c1_mb_c = vstack([c1_mb_c, cdt1])
    
nbins_z = 25
nbins_mag = 25
z_bins=np.linspace(0,max(cosmo_values['z']),25)
mag_bins=np.linspace(13,28,25)
n_gal, z, m = np.histogram2d(c1_mb_c['z'], c1_mb_c['mag_i'], bins=[nbins_z,nbins_mag], range=[[z_bins[0],z_bins[nbins_z-1]],[mag_bins[0],mag_bins[nbins_mag-1]]])
density = n_gal/439.78986

rich = []
for id in c1['id']:
    rich.append(np.sum(c1_mb_c[c1_mb_c['id_cluster']==id]['pmem']))
bins = np.linspace(0,200, num=101)
plt.hist(c1['mass'], alpha = .5, bins=bins, label = "LAMBSTAR")

plt.hist(rich, alpha = .5, bins=bins, label = "By hand")
plt.legend()
plt.savefig("/pbs/home/n/namourou/fig1.png")
