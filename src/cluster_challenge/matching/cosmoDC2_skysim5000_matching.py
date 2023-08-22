#!/usr/bin/env python
# coding: utf-8

###Original notebook from M. Aguena, all credit to him!
###Source: http://dev.linea.org.br/~aguena/share/desc/get_halos_full.html
###Adapted by T .Guillemin

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.table import vstack
from matplotlib.ticker import MultipleLocator
import healpy as hp
import GCRCatalogs
from GCR import GCRQuery
import sys
import os
import shutil

outpath = 'full_halos/'
if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.makedirs(outpath)

#load catalogs
#full_with_photozs_v1pzdat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_v1')
#full_with_photozs_v1pzb = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_flexzboost_v1')
#full_with_photozs_v1pzdat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_v1')
full_with_photozs_v1pzdat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image_with_photozs_v1')

full_pixels32 = np.loadtxt('full_pix32.txt', dtype=int)
from astropy.table import vstack

#not used
#cols_halo = [
#    'halo_id', 'halo_mass', 
#    'ra_true', 'dec_true', 'ra', 'dec',
#    'redshift_true', 'galaxyID', 'galaxy_id', 'redshift']
#    #+['photoz_mask', 'photoz_mean', 'photoz_median', 'photoz_mode']

halo_cols =[
    'halo_id', 'halo_mass',
    'baseDC2/target_halo_fof_halo_id',
    'ra_true', 'dec_true', 'ra', 'dec', 'redshift_true',
    'baseDC2/target_halo_redshift', 'galaxyID', 'redshift'
]

hcats = []
for p in full_pixels32:
      
    #print(p)

    #if(p!=10200):
    #    continue

    data = full_with_photozs_v1pzdat.get_quantities(
        halo_cols, native_filters=[f'healpix_pixel == {p}'],
        filters=['halo_mass > 10**(12.885)', 'is_central==True'])
    print(f'    {len(data["galaxyID"]):,}')
    
    temp_cat = Table()
    for col in halo_cols:
        temp_cat[col] = data[col]
        if col in ['halo_id', 'galaxyID']:
            temp_cat[col] = temp_cat[col].astype(int)
        
    hcats.append(temp_cat)
    del temp_cat

dc2_halos = vstack(hcats)
del hcats
dc2_halos.write(outpath+'halos_m12.885.fits')

#load skysim
#skysim = GCRCatalogs.load_catalog('skysim5000_v1.1.1_small')
skysim = GCRCatalogs.load_catalog('skysim5000_v1.1.1')
print('skysim loaded')
ss_cols = [
    'halo_id', 'halo_mass',
    'baseDC2/target_halo_fof_halo_id',
    'ra_true', 'dec_true', 'redshift_true',
    'ra', 'dec', 'redshift',
    'baseDC2/sod_halo_mass', 'baseDC2/sod_halo_radius', 'baseDC2/target_halo_redshift',
    ]

skysim_dat = Table(skysim.get_quantities(
    ss_cols, filters=['halo_mass > 10**(12.88)', 'is_central==True']))

for c in skysim_dat.colnames:
    if 'halo_mass' in c:
        skysim_dat[c].info.format = '.4g'
skysim_dat['m200c'] = skysim_dat['baseDC2/sod_halo_mass']/skysim.cosmology.h
skysim_dat['r200c'] = skysim_dat['baseDC2/sod_halo_radius']/skysim.cosmology.h
print(f'{len(skysim_dat):,} -> {(skysim_dat["m200c"]>1e13).sum():,}')
skysim_dat.write(outpath+'skysim5000_halos_m12.88.fits')
skysim_dat[skysim_dat["m200c"]>1e13].write(outpath+'skysim5000_halos_m200c_13.0.fits')

#Match halos on ID
#There can be repeated IDs, so it is better to add an ID list to the dictionary instead of the rownumber
cosmoDC2_dict = {}
for i, ID in enumerate(dc2_halos['baseDC2/target_halo_fof_halo_id']):
        cosmoDC2_dict[ID] = cosmoDC2_dict.get(ID, [])+[i]
skysim_dat['matched'] = np.array([i in cosmoDC2_dict for i in skysim_dat['baseDC2/target_halo_fof_halo_id']])
skysim_dat['mt_ids'] = None
for i, ID in enumerate(skysim_dat['baseDC2/target_halo_fof_halo_id']):
    skysim_dat['mt_ids'][i] = cosmoDC2_dict.get(ID, [])

print(f"""
Unique cosmoDC2 ids: {len(cosmoDC2_dict):,}
SkySim halos in cosmoDC2: {(skysim_dat['matched']).sum():,}
""")

skysim_dict = {}
for i, ID in enumerate(skysim_dat['baseDC2/target_halo_fof_halo_id']):
    skysim_dict[ID] = skysim_dict.get(ID, [])+[i]
dc2_halos['matched'] = np.array([i in skysim_dict for i in dc2_halos['baseDC2/target_halo_fof_halo_id']])
dc2_halos['mt_ids'] = None
for i, ID in enumerate(dc2_halos['baseDC2/target_halo_fof_halo_id']):
    dc2_halos['mt_ids'][i] = skysim_dict.get(ID, [])

print(f"""
SkySim catalog size: {len(skysim_dat):,}
Unique SkySim ids: {len(skysim_dict):,}
cosmoDC2 halos in SkySim: {(dc2_halos['matched']).sum():,}
""")

dc2_halos['mt_id_final'] = None
for i, dc2h in enumerate(dc2_halos):
    if len(dc2h['mt_ids'])==1:
        dc2_halos['mt_id_final'][i] = dc2h['mt_ids'][0]
    elif len(dc2h['mt_ids'])>1:
        sksh = skysim_dat[dc2h['mt_ids']]
        mtmask = (dc2h['halo_mass']==sksh['halo_mass'])*\
                 (abs(dc2h['ra_true']-sksh['ra_true'])<.1)*\
                 (abs(dc2h['dec_true']-sksh['dec_true'])<.1)
        if mtmask.sum()==1:
            j = np.arange(mtmask.size, dtype=int)[mtmask][0]
            #print(dc2h['mt_ids'], mtmask, j)
            dc2_halos['mt_id_final'][i] = dc2h['mt_ids'][j]
        else:
            print(mtmask)

print(len(dc2_halos))
print((dc2_halos['mt_id_final']!=None).sum())

print(dc2_halos[dc2_halos['mt_id_final']==None][:10])

from astropy.coordinates import SkyCoord
from astropy import units as u
skcoord = lambda cat: SkyCoord(cat['ra_true']*u.deg, cat['dec_true']*u.deg, frame='icrs')

dc2_halos['SkyCoord'] = skcoord(dc2_halos)

skysim_dat['SkyCoord'] = skcoord(skysim_dat)

dc2_halos['dist_deg'] = 999.
dc2_halos_mt = dc2_halos[dc2_halos['mt_id_final']!=None]
dc2_halos['dist_deg'][dc2_halos['mt_id_final']!=None] = dc2_halos_mt['SkyCoord'].separation(
    skysim_dat['SkyCoord'][list(dc2_halos_mt['mt_id_final'])]).value
del dc2_halos_mt

dc2_halos['dist_deg'][dc2_halos['mt_id_final']!=None].max()

len(dc2_halos), (dc2_halos['mt_id_final']!=None).sum(), (dc2_halos['dist_deg']<1).sum()

dc2test = dc2_halos[dc2_halos['dist_deg']>1][
    dc2_halos['halo_mass'][dc2_halos['dist_deg']>1].argmax()]
dc2test

skysim_dat[
    (dc2test['SkyCoord'].separation(skysim_dat['SkyCoord']).value<.1)*\
    (abs(dc2test['redshift_true']-skysim_dat['redshift_true'])<.1)
    ]

dc2_halos2 = dc2_halos[dc2_halos['halo_mass']>1e13]
print(len(dc2_halos2), (dc2_halos2['mt_id_final']!=None).sum(), (dc2_halos2['dist_deg']<1).sum())
del dc2_halos2

dc2_halos['m200c'] = -1.
dc2_halos['r200c'] = -1.

dc2_halos['m200c'][dc2_halos['dist_deg']<1] = skysim_dat['m200c'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
dc2_halos['r200c'][dc2_halos['dist_deg']<1] = skysim_dat['r200c'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])] 

dc2_halos[dc2_halos['dist_deg']>1]

closest_pair = []
for h in dc2_halos[dc2_halos['dist_deg']>1]:
    mass_diff = abs(np.log10(h['halo_mass'])-np.log10(dc2_halos['halo_mass']))
    z_diff = abs(h['redshift_true']-dc2_halos['redshift_true'])
    dist = 5*mass_diff+z_diff
    dist[h['halo_id']==dc2_halos['halo_id']] = 9999.
closest_pair.append(dist.argmin())

samp_orig, samp_closest = dc2_halos[dc2_halos['dist_deg']>1], dc2_halos[closest_pair]

for col in ('m200c', 'r200c'):
    dc2_halos[col][dc2_halos['dist_deg']>1] = dc2_halos[col][closest_pair]

dc2_halos['skysim_ind'] = -1
dc2_halos['skysim_ind'][dc2_halos['mt_id_final']!=None] = dc2_halos['mt_id_final'][dc2_halos['mt_id_final']!=None]

dc2_halos['skysim_inds'] = [','.join(np.array(inds, dtype=str)) for inds in dc2_halos['mt_ids']]

dc2_halos['skysim_halo_id'] = -1
dc2_halos['skysim_halo_id'][dc2_halos['dist_deg']<1] = skysim_dat['halo_id'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
#add other variables from skysim
dc2_halos['skysim_ra'] = -1
dc2_halos['skysim_dec'] = -1
dc2_halos['skysim_redshift'] = -1
dc2_halos['skysim_ra_true'] = -1
dc2_halos['skysim_dec_true'] = -1
dc2_halos['skysim_redshift_true'] = -1
dc2_halos['skysim_ra'][dc2_halos['dist_deg']<1] = skysim_dat['ra'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
dc2_halos['skysim_dec'][dc2_halos['dist_deg']<1] = skysim_dat['dec'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
dc2_halos['skysim_redshift'][dc2_halos['dist_deg']<1] = skysim_dat['redshift'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
dc2_halos['skysim_ra_true'][dc2_halos['dist_deg']<1] = skysim_dat['ra_true'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
dc2_halos['skysim_dec_true'][dc2_halos['dist_deg']<1] = skysim_dat['dec_true'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]
dc2_halos['skysim_redshift_true'][dc2_halos['dist_deg']<1] = skysim_dat['redshift_true'][list(dc2_halos['mt_id_final'][dc2_halos['dist_deg']<1])]

print(dc2_halos[:5])

dc2_halos.colnames

dc2_halos[
    ['halo_id',
     'halo_mass',
     'baseDC2/target_halo_fof_halo_id',
     'ra_true',
     'dec_true',
     'ra',
     'dec',
     'redshift_true',
     'baseDC2/target_halo_redshift',
     'galaxyID',
     'skysim_halo_id',
     'skysim_ind',
     'skysim_inds',
     'dist_deg',
     #'SkyCoord',
     'm200c',
     'r200c',
     ]
    ].write(outpath+'halos_m12.885_mtskysim.fits')

dc2_halos['mass_fof'] = dc2_halos['halo_mass']

dc2_halos[
    ['halo_id',
     #'ra_true',
     #'dec_true',
     'ra',
     'dec',
     #'redshift_true',
     'redshift',
     'mass_fof',
     'm200c',
     #'r200c',
     'skysim_halo_id',
     ]
    ][dc2_halos['m200c']>1e13].write(outpath+'halos_m200c_13.0.fits')


print(len(dc2_halos))
print((dc2_halos['m200c']>1e13).sum())
print(dc2_halos.colnames)
print('Halos whith m200c>1e13 and not matched:', ((dc2_halos['m200c']>1e13)*(dc2_halos['dist_deg']>1)).sum())

#reduce table to a few variables
dc2_halos_reduced = dc2_halos['halo_id','ra','dec','redshift','ra_true','dec_true','redshift_true','mass_fof','m200c','r200c','skysim_halo_id'][dc2_halos['m200c']>1e13]# 'skysim_ra','skysim_dec','skysim_redshift','skysim_ra_true','skysim_dec_true','skysim_redshift_true']
dc2_halos_reduced.write(outpath+'halos_m200c_13.0_reduced.fits')
print(dc2_halos_reduced.colnames)
print(dc2_halos_reduced[:5])

#extra step: renaming columns to be in synch with the redmapper validation project
#file_mricci = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/debug/matched_pairs_Mfofcut.fits'
#data = Table.read(file_mricci)
#print(data.colnames)
#'cat2_id', 'cat2_z', 'cat2_ra', 'cat2_dec', 'cat2_mass_fof', 'cat2_M200c' 

#files written
#halos_m12.885.fits
#skysim5000_halos_m12.88.fits
#skysim5000_halos_m200c_13.0.fits
#halos_m12.885_mtskysim.fits
#halos_m200c_13.0.fits
#halos_m200c_13.0_reduced.fits

sys.exit()
