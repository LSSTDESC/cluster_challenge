

import numpy as np
from astropy import units as u
from astropy.io import fits
import sys
import os
import json
import shutil

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import MembershipMatch
from clevar.cosmology import AstroPyCosmology
from IPython.display import display



requested_cats = [sys.argv[1].split('_'), sys.argv[2].split('_')]
matching_method = sys.argv[3]
matching_preference = sys.argv[4]

inpath  = '/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/'
outpath = '/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/'



## COLLECT THE CATALOGS TO MATCH.
cats = []
if np.any([c[0] == 'cosmoDC2' for c in requested_cats]) :
	matched_to_cosmoDC2 = True
else :
	matched_to_cosmoDC2 = False


for c in requested_cats :
	if c[0] == 'cosmoDC2' :
		cltags = {
			'id':'id_cl',
			'ra':'ra_cl',
			'dec':'dec_cl',
			'z':'z_cl',
			'mass':'mass',
			'log_mass':'log10_mass'}
		mbtags = {
			'id':'id_mb',
			'id_cluster':'clid_mb',
			'ra':'ra_mb',
			'dec':'dec_mb',
			'z':'z_mb',
			'pmem':'pmem',
			'mag_g':'mag_g',
			'mag_r':'mag_r',
			'mag_i':'mag_i',
			'mag_z':'mag_z',
			'mag_y':'mag_y'}

		print(f"loading from {inpath}halos/cosmoDC2/{c[1]}/...")
		if os.path.exists(inpath + f'halos/cosmoDC2/{c[1]}/') :
			cat = ClCatalog.read(inpath + f'halos/cosmoDC2/{c[1]}/Catalog.fits', f'cosmoDC2.{c[1]}', tags=cltags)
			#cat.read_members(inpath + f'halos/cosmoDC2/{c[1]}/Catalog_members.fits', tags=mbtags)
			cats.append(cat)
		else :
			sys.exit(f'The version {c[1]} of cosmoDC2 does not exist.')
		print(f"{'.'.join(c)} catalog is loaded...")
	else :
		cltags = {
			'id':'id_cl',
			'ra':'ra_cl',
			'dec':'dec_cl',
			'z':'z_cl',
			'mass':'mass',
			#'mass_err':'mass_err',
			#'snr':'snr_cl'}
        }
		mbtags = {
			'id':'id_mb',
			'id_cluster':'clid_mb',
			'ra':'ra_mb',
			'dec':'dec_mb',
			'z':'z_mb',
			'pmem':'pmem',
			'mag_g':'mag_g',
			'mag_r':'mag_r',
			'mag_i':'mag_i',
			'mag_z':'mag_z',
			'mag_y':'mag_y'}
		
		print(f"loading from {inpath}{c[0]}/{c[1]}/{c[2]}/...")
		if os.path.exists(inpath + f'{c[0]}/{c[1]}/{c[2]}/') :
			cat = ClCatalog.read(inpath + f'{c[0]}/{c[1]}/{c[2]}/Catalog.fits', f'{c[0]}.{c[1]}.{c[2]}', tags=cltags)
			#cat.read_members(inpath + f'{c[0]}/{c[1]}/{c[2]}/Catalog_members.fits', tags=mbtags)
			cats.append(cat)
		else :
			sys.exit(f"The catalog {'_'.join(c)} does not exist.")
		print(f"{'.'.join(c)} catalog is loaded...")

## SET OUTPUT PATH	
## OUTPATH DIRECTORY STRUCTURE:
##	after_matching/
##		|- cosmoDC2_wazp.cosmoDC2.truez/
##			|- v0_v0/
##				|- member_matching/
##				|- proximity_matching/
##					|- deltaz_0.05_matchradius_1.0mpc_pref_angular_proximity
##					|- deltaz_0.05_matchradius_1.0mpc_pref_more_massive
##		|- cosmoDC2_wazp.cosmoDC2.fzb/
##		|- cosmoDC2_wazp.DC2/


## SET ORDERING OF CATALOG NAMBES TO BE CONSISTENT IN DIRECTORIES.
if np.any([c[0] == 'cosmoDC2' for c in requested_cats]) :
        first = np.argwhere(np.array([c[0] == 'cosmoDC2' for c in requested_cats]))[0][0]
        last  = np.argwhere(np.array([c[0] != 'cosmoDC2' for c in requested_cats]))[0][0]
else :
        (first, last) = np.argsort([c[0] for c in requested_cats])
outpath += f"{'.'.join(requested_cats[first][:-1])}_{'.'.join(requested_cats[last][:-1])}/"
outpath += f"{requested_cats[first][-1]}_{requested_cats[last][-1]}/"



## PERFORM PROXIMITY MATCHING
if matching_method == 'proximity' :
	mt = ProximityMatch()

	delta_z = float(sys.argv[5])#0.05
	match_radius =  float(sys.argv[6])#1 ## in Mpc
	match_config = {
		'type':'cross',				## OPTIONS: cross, cat1, cat2
		'which_radius':'max',			## OPTIONS: cat1, cat2, min, max
		'preference':matching_preference,	## OPTIONS: more_massive, angular_proximity, redshift_proximity
		'catalog1':{'delta_z':delta_z, 'match_radius':f'{match_radius} mpc'},
		'catalog2':{'delta_z':delta_z, 'match_radius':f'{match_radius} mpc'}
		}

	outpath += f'proximity_matching/deltaz_{delta_z}_matchradius_{match_radius}mpc_pref_{matching_preference}/'

	if not os.path.exists(outpath) :
		os.makedirs(outpath)
	print(f'OUTPATH = {outpath}')

	cats_raw = [cats[0].raw(), cats[1].raw()]
	cosmo = AstroPyCosmology()

	print('MATCHING...')
	mt.match_from_config(cats_raw[0], cats_raw[1], match_config, cosmo=cosmo)
	cats_raw[0].cross_match()
	cats_raw[1].cross_match()

	print('WRITING MATCHED FILES...')
	cats_raw[0].write(outpath + f'{cats[0].name}.fits', overwrite=True)
	cats_raw[1].write(outpath + f'{cats[1].name}.fits', overwrite=True)
	mt.save_matches(cats_raw[0], cats_raw[1], out_dir=outpath, overwrite=True)

	#to print summary
	mt.load_matches(cats_raw[0], cats_raw[1], out_dir=outpath)
	#display(cats_raw[0])
	#display(cats_raw[1])

## PERFORM MEMBER MATCHING
if matching_method == 'member' :
	mt = MembershipMatch()
	
	minimum_share_fraction = float(sys.argv[5])
	match_config = {
	  	'type':'cross',				## OPTIONS: cross, cat1, cat2
	  	'preference':matching_preference,	## OPTIONS: shared_member_fraction, more_massive, angular_proximity, redshift_proximity
	  	'minimum_share_fraction':minimum_share_fraction,
	  	'match_members_kwargs': {'method':'id'},
	  	}
	
	outpath += f'member_matching/fshare_{minimum_share_fraction}_pref_{matching_preference}/'
	
	if not os.path.exists(outpath) :
		os.makedirs(outpath)
	print(f'OUTPATH = {outpath}')
	
	print('MATCHING...')
	mt.match_from_config(cats[0], cats[1], match_config)
	
	## TO INVESTIGATE MATCHING ISSUES
	mt.fill_shared_members(cats[0], cats[1])
	mt.save_shared_members(cats[0], cats[1], fileprefix=outpath+'mem_share')
	
	cats[0].cross_match()
	cats[1].cross_match()
	
	print('WRITING...')
	cats[0].write(outpath + f'{cats[0].name}.fits', overwrite=True)
	cats[1].write(outpath + f'{cats[1].name}.fits', overwrite=True)
	mt.save_matches(cats[0], cats[1], out_dir=outpath, overwrite=True)
	
	## TO PRINT SUMMARY
	mt.load_matches(cats[0], cats[1], out_dir=outpath)
	#display(cats[0])
	#display(cats[1])

sys.exit()
