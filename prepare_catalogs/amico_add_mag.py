print('Adding magnitudes in AMICO catalog')
#!/usr/bin/env python
# coding: utf-8
import numpy as np
from astropy.table import Table, vstack
import healpy as hp
import sys
import os
import time

start = time.time()


#---Paths---#

inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/cosmoDC2_photoz_flexzboost_v1_iband/DETECTIONS_DERIVED/'
outpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/'
print('\nCatalog loaded at', inpath, 'and will be saved in', outpath)

curr_tile = str(sys.argv[1])
print('\nprocessing', curr_tile)

if not os.path.exists(outpath):
    os.makedirs(outpath)


#---Reading tables---#

mbn = Table.read(inpath + str(curr_tile) + '_map_associations_noBuffer.fits')

tile_t = Table.read('/sps/lsst/groups/clusters/amico_validation_project/catalogs/CosmoDC2/cosmodc2_neighbours.fits')

tile_l = tile_t[tile_t['tile']==int(curr_tile)]['list_of_neighbour_tiles'][0].split(',') #Because neighbours listed w/ comas in a long string, need to convert it to a list
for i in range(len(tile_l)):
    tile_l[i] = int(tile_l[i])
print('Healpix list', ':', tile_l, '\nNow merging all input healpix files')

healpath = '/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/'
large_heal = Table.read(healpath + str(tile_l[0]) + '/galaxies.fits')
for tile in tile_l[1:]:
    heal = Table.read(healpath + str(tile) + '/galaxies.fits')
    large_heal = vstack([large_heal, heal])
print('Done. \nIncorporation of magnitudes')

  
#---Add Magnitudes---#

heal_dict = {}
for i, ID in enumerate(large_heal['galaxy_id']):
    heal_dict[ID] = heal_dict.get(ID, [])+[i]
mbn['matched'] = np.array([i in heal_dict for i in mbn['GALID']])
mbn['mt_ids'] = None
for i, ID in enumerate(mbn['GALID']):
    mbn['mt_ids'][i] = heal_dict.get(ID, [])

mbn['mag_g'] = 0.0
mbn['mag_i'] = 0.0
mbn['mag_r'] = 0.0
mbn['mag_z'] = 0.0
mbn['mag_y'] = 0.0
mags = ['mag_g','mag_i','mag_r','mag_z','mag_y']
per_list = [.1,.2,.4,.6,.8,1.0]
k = 0
mbn['mt_id_final'] = None

for i, amg in enumerate(mbn):
    per = round(i/len(mbn),2)
    if per in per_list and k != per:
        print(i, per)
        k = per
    if len(amg['mt_ids'])==1:
        mbn['mt_id_final'][i] = amg['mt_ids'][0]
        healsh = large_heal[amg['mt_ids']]
        for mag in mags:
            mbn[mag][i] =  healsh[mag]
    elif len(amg['mt_ids'])>1:
        healsh = large_heal[amg['mt_ids']]
        j = np.arange(healsh.size, dtype=float)[0]
        #print(dc2h['mt_ids'], mtmask, j)
        mbn['mt_id_final'][i] = amg['mt_ids'][j]
        for mag in mags:
            mbn[mag][i] =  healsh[mag][j]
            
            
print(mbn[:5])
print('There are :', len(mbn[mbn['mag_g']==0]['mag_g']), 'galaxies which are not found inside input catalog')

end = time.time()
print('Program ran in :', end-start)


#---Saving---#

forsave = mbn['GALID','FIELD_PROB','ASSOC_ID','ASSOC_PROB', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y']
forsave.write(outpath + str(curr_tile) + '_map_associations_w_mag'  + '.fits', overwrite = True)

print('Done.')

sys.exit()