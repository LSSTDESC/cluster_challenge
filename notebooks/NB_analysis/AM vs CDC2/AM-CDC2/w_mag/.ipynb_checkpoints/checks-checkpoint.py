#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, hstack, vstack
from clevar.catalog import ClCatalog
from numpy.ma import masked
import sys

inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/amico_cosmoDC2/mag_i/w_rich/"
outpath = inpath
#am_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/Catalog_members.fits' 
#cdc_inpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4_image/Catalog_members.fits"

print('reading tables')
c1 = ClCatalog.read(inpath + 'c1_p.fits', 'c1', full = True)
c1 = c1[c1['mt_cross']!= None]
c2 = ClCatalog.read(inpath + 'c2_p.fits', 'c2', full = True)
c2 = c2[c2['mt_cross']!= None]
c1.read_members(inpath + 'c1_p_members.fits', full = True)
c2.read_members(inpath + 'c2_p_members.fits', full = True)
#c2.members = c2.members[c2.members['mag_i']<=25.3]
print("Done, now it's just plotting", max(c1.members['mag_i']), max(c2.members['mag_i']), len(c1), len(c2))

bin1 = [0.2,0.5,0.8,1.0,1.2,1.5,1.8]
bin2 = [13,13.5]
for i in range(0,len(bin1)-1):
    for j in range(0,len(bin2)-1):
        z_min = bin1[i]
        z_max = bin1[i+1]
        m_min = bin2[j]
        m_max = bin2[j+1]
        
        c1_b = c1[(c1['z']>=0z_min)*(c1['z']<=z_max)]
        am_mb_cdt = c1_b.members
        am_mb_delt = am_mb_cdt['mag_r'] - am_mb_cdt['mag_i']
        
        c2_b = c2[(c2['z']>=z_min)*(c2['z']<=z_max)*(c2['log_mass']<=m_max)*(c2['log_mass']>=m_min)]
        cdc_mb_cdt = c2_b.members
        cdc_mb_delt = cdc_mb_cdt['mag_r'] - cdc_mb_cdt['mag_i']

        x_bins = np.linspace(16,26,100)
        y_bins = np.linspace(-0.5,2,60)

        am_hist = np.histogram2d(am_mb_cdt['mag_i'], am_mb_delt, bins = (x_bins,y_bins), weights = am_mb_cdt['pmem'])
        cdc_hist = np.histogram2d(cdc_mb_cdt['mag_i'], cdc_mb_delt, bins = (x_bins,y_bins))
        am_hist = am_hist[0]
        cdc_hist = cdc_hist[0]
        vmax1 = np.max(am_hist)
        vmax2 = np.max(cdc_hist)
        vmax = (vmax1 + 3*vmax2)/4
        print(vmax1, vmax2, vmax)
        am_hist[am_hist==0] = np.nan
        cdc_hist[cdc_hist==0] = np.nan

        x, y = np.meshgrid(x_bins, y_bins)
        fig, (ax1, ax2) = plt.subplots(1,2, figsize = (14,5))
        print('Plotting...')
        c = ax1.pcolormesh(x, y, am_hist.T, cmap='jet', vmin=0, vmax = vmax)
        ax1.set_xlabel('mag$_i$', fontsize = 13)
        ax1.set_ylabel('mag$_r$ - mag$_i$', fontsize = 13)
        ax1.set_title('AMICO ' + str(z_min) + '<z<' + str(z_max), fontsize = 13)
        cc1 = ax2.pcolormesh(x, y, cdc_hist.T, cmap='jet', vmin=0, vmax = vmax)
        ax2.set_xlabel('mag$_i$', fontsize = 13)
        ax2.set_ylabel('mag$_r$ - mag$_i$', fontsize = 13)
        ax2.set_title('CosmoDC2 ' + str(z_min) + '<z<' + str(z_max) + ';' + str(m_min) + '<m<' + str(m_max), fontsize = 13)
        ax1.set_xlim([16,26])
        ax1.set_ylim([-.5,2])
        ax2.set_xlim([16,26])
        ax2.set_ylim([-.5,2])
        fig.colorbar(c, ax=ax1, label = 'Nb of members')
        fig.colorbar(cc1, ax=ax2, label = 'Nb of members')
        plt.savefig('/pbs/home/n/namourou/test_jupyter/cluster_challenge/plots/amico_cosmoDC2/mag_i/p_matching/dmagvsmag/w_rich/' + str(z_min) + '-' + str(z_max) + '_' + str(m_min) + '-' + str(m_max) + '.png')
