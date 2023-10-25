#!/usr/bin/env python
# coding: utf-8
print('Reading AMICO files')
###import
import GCRCatalogs
import numpy as np
from astropy.table import Table
import sys
import os
import pandas as pd

#cut = 3.5
#print("performing cut" + str(cut))
# ---Paths---#

#CosmoDC2
cl_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/cosmoDC2_photoz_flexzboost_v1_yband/map_detections_refined_noBuffer_all_noDoubles_handFix.fits'
mb_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_y/all_maps.fits'
outpath = '/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/amico/cosmoDC2.fzb/magy/'

#DC2
#cl_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/DC2_v0_yband/map_detections_refined_noBuffer_all_noDoubles.fits'
#mb_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/DC2/mag_y/all_maps.fits'
#outpath = '/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/amico/DC2.fzb.magy/sn' + str(cut) + '/'
#print('Reading clusters at', cl_inpath, '\nReading members at', mb_inpath, '\nSaving table at', outpath)


# ---Reading tables---#

amico_cl_data = Table.read(cl_inpath)['ID', 'Xphys', 'Yphys', 'Zphys','LAMBSTAR', 'SN_NO_CLUSTER', 'SN', 'UID', 'TILE']    
c1 = Table([amico_cl_data['UID'],amico_cl_data['Xphys'],amico_cl_data['Yphys'],amico_cl_data['Zphys'],amico_cl_data['LAMBSTAR'], amico_cl_data['SN']],names=('id_cl','ra_cl','dec_cl','z_cl','mass','snr_cl'))
print("Number of AMICO clusters = ", len(amico_cl_data))

mb_amico_data = Table.read(mb_inpath)
c1_members = Table([mb_amico_data['mb_id'],mb_amico_data['uid'].astype(int),mb_amico_data['ra'],mb_amico_data['dec'],mb_amico_data['z'],mb_amico_data['prob'], mb_amico_data['f_prob'], mb_amico_data['lambstar'], mb_amico_data['mag_g'], mb_amico_data['mag_i'], mb_amico_data['mag_r'], mb_amico_data['mag_z'], mb_amico_data['mag_y'],mb_amico_data['SN']], names=('id_mb','clid_mb','ra_mb','dec_mb','z_mb','pmem', 'f_prob', 'lambstar', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y', 'snr_cl'))


# ---Cut---#


#z_max = 2
#sn_min = cut
#ls_min = 38
#c1 = c1[(c1['sn']>=sn_min)]
#c1_members = c1_members[(c1_members['sn']>=sn_min)]


print('First 5 elements of cluster list and member list :\n', c1[:5])
print(c1_members[:5])


# ---Saving---#

c1.write(outpath + 'Catalog.fits', overwrite=True)
c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()
