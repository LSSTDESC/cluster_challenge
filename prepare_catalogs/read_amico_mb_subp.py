#!/usr/bin/env python
# coding: utf-8
print('Reading amico files before merging')
###import
import numpy as np
from astropy.table import Table, vstack
import sys
import os
import time

start = time.time()
tile = str(sys.argv[1])
print('Reading tile = ', str(sys.argv[1]))

#---Paths---#

cl_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/cosmoDC2_photoz_flexzboost_v1_iband/map_detections_refined_noBuffer_all_noDoubles.fits'

cat_path = ('/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/' + tile + '_map_associations_w_mag' + '.fits')

outpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/'
if not os.path.exists(outpath):
    os.makedirs(outpath)
print('Reading clusters at', cl_inpath, '\nReading members at', cat_path, '\nSaving table at', outpath)  

#---Reading tables---#


amico_cl_data = Table.read(cl_inpath)['ID', 'Xphys', 'Yphys', 'Zphys','LAMBSTAR', 'SN_NO_CLUSTER', 'SN', 'UID', 'TILE']  

table = Table.read(cat_path)

print('Number of clusters in tract = ', len(amico_cl_data[amico_cl_data['TILE'] == int(tile)]))


#---Settings before creating the table---#

tiles = [10066,10067,10068,10069,10070,10071,10072,10073,10074,10193,10194,10195,10196,10197,10198,10199,10200,10201,10202,10321,10322,10323,10324,10325,10326,10327,10328,10329,10444,10445,10446,10447,10448,10449,10450,10451,10452,8786,8787,8788,8789,8790,8791,8792,8793,8794,8913,8914,8915,8916,8917,8918,8919,8920,8921,9042,9043,9044,9045,9046,9047,9048,9049,9050,9169,9170,9171,9172,9173,9174,9175,9176,9177,9178,9298,9299,9300,9301,9302,9303,9304,9305,9306,9425,9426,9427,9428,9429,9430,9431,9432,9433,9434,9554,9555,9556,9557,9558,9559,9560,9561,9562,9681,9682,9683,9684,9685,9686,9687,9688,9689,9690,9810,9811,9812,9813,9814,9815,9816,9817,9818,9937,9938,9939,9940,9941,9942,9943,9944,9945,9946]

new_table = {'mb_id' : [], 'cl_id' : [], 'prob' : [], 'f_prob' : [], 'uid' : [], 
             'ra' :[], 'dec':[],'z':[], 'lambstar' : [], 'SN' : [], 
             'mag_g' : [], 'mag_i' : [], 'mag_r' : [], 'mag_z' : [], 'mag_y' : []}


#---Creating the table---#

for cluster_id in amico_cl_data[amico_cl_data['TILE']==int(tile)]['ID']: 
    if cluster_id in table['ASSOC_ID']:
        table_cdt = table[(table['ASSOC_ID']==cluster_id)]
        mb_id, cl_id, prob, f_prob = table_cdt['GALID'], table_cdt['ASSOC_ID'], table_cdt['ASSOC_PROB'], table_cdt['FIELD_PROB']
        mag_g, mag_i, mag_r, mag_z, mag_y = table_cdt['mag_g'], table_cdt['mag_i'], table_cdt['mag_r'], table_cdt['mag_z'], table_cdt['mag_y']
        amico_cl = amico_cl_data[(amico_cl_data['TILE'] == int(tile)) * (amico_cl_data['ID'] == cluster_id)]
        lenght = np.ones(len(table_cdt['ASSOC_ID']))
        uid, ra, dec, z, lambstar, sn = int(amico_cl['UID'][0])*lenght, amico_cl['Xphys'][0]*lenght, amico_cl['Yphys'][0]*lenght, amico_cl['Zphys'][0]*lenght, amico_cl['LAMBSTAR'][0]*lenght, amico_cl['SN'][0]*lenght
        new_table['mb_id'] += list(mb_id)
        new_table['cl_id'] += list(cl_id)
        new_table['prob'] += list(prob) 
        new_table['f_prob'] += list(f_prob)
        new_table['uid'] += list(uid)
        new_table['ra'] += list(ra)
        new_table['dec'] += list(dec)
        new_table['z'] += list(z)
        new_table['lambstar'] += list(lambstar)
        new_table['SN'] += list(sn)
        new_table['mag_g'] += list(mag_g)
        new_table['mag_i'] += list(mag_i)
        new_table['mag_r'] += list(mag_r)
        new_table['mag_z'] += list(mag_z)
        new_table['mag_y'] += list(mag_y)
        
end = time.time()
print('Catalog created in', (end-start), 'seconds')

print('created table lenght = ', len(new_table['mb_id']))


#---Saving---#

Table(new_table).write(outpath + tile + '_map_associations.fits', overwrite=True)
print(new_table['mb_id'][:5], new_table['uid'][:5])
print('--------done--------')

sys.exit()
