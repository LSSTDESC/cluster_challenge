

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


inpath  = '/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/'
outpath = '/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/'


## SELECT THE CATALOGS TO MATCH.
available_runons = ['cosmoDC2', 'DC2']
available_algos  = ['wazp', 'redmapper', 'amico']
available_matching_methods = ['member', 'proximity']
try :
	algo	= sys.argv[1]
	runon	= sys.argv[2]
	algo_runon_v	= sys.argv[3]
	matching_method = sys.argv[4]
except ValueError :
	print(f'Invalid argument selection. Please check.')

if not (np.isin(runon, available_runons) & np.isin(algo, available_algos) & np.isin(matching_method, available_matching_methods)) :
	sys.exit('Invalid argument selection. Please check.')



## COLLECT THE CATALOGS TO MATCH.
cats = []


## COSMODC2 (IT IS ASSUMED THAT ALL MATCHES WILL BE MADE TO FULL COSMODC2 i.e. v1)
if True :
	cosmoDC2_v = 'v1'	## FOR COSMODC2.SMALL, CHANGE TO v0
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
		'pmem':'pmem'}

	if os.path.exists(inpath + f'cosmoDC2/halos/{cosmoDC2_v}/') :
		cat = ClCatalog.read(inpath + f'cosmoDC2/halos/{cosmoDC2_v}/Catalog.fits', f'cosmoDC2_{cosmoDC2_v}', tags=cltags)
		cat.read_members(inpath + f'cosmoDC2/halos/{cosmoDC2_v}/Catalog_members.fits', tags=mbtags)
		cats.append(cat)
	else :
		sys.exit(f'The version {cosmoDC2_v} of cosmoDC2 does not exist. Please see {inpath} for available versions.')
	
## PREPARE CLFINDER CATALOG
cltags = {
	'id':'id_cl',
	'ra':'ra_cl',
	'dec':'dec_cl',
	'z':'z_cl',
	'mass':'mass'}
mbtags = {
	'id':'id_mb',
	'id_cluster':'clid_mb',
	'ra':'ra_mb',
	'dec':'dec_mb',
	'z':'z_mb',
	'pmem':'pmem'}

if os.path.exists(inpath + f'{runon}/{algo}/{algo_runon_v}/') :
	cat = ClCatalog.read(inpath + f'{runon}/{algo}/{algo_runon_v}/Catalog.fits', f'{algo}.{runon}.{algo_runon_v}', tags=cltags)
	cat.read_members(inpath + f'{runon}/{algo}/{algo_runon_v}/Catalog_members.fits', tags=mbtags)
	cats.append(cat)
else :
	sys.exit(f'The version {algo_runon_v} of {algo}.{runon} does not exist. Please see {inpath} for available versions.')






## OUTPATH DIRECTORY STRUCTURE: (ASSUMING THERE ARE NOT COSMODC2 - DC2 COMPARISONS)
##	after_matching/
##		|- cosmoDC2/
##			|
##			|- cosmoDC2_wazp.cosmoDC2/
##				|- v0_v0/
##				|- v1_v0/
##			|- cosmoDC2_redmapper.cosmoDC2/
##			|- cosmoDC2_amico.cosmoDC2/
##		|- DC2/
##			|- cosmoDC2_wazp.DC2/
##				|- v0_v0/
##			|- cosmoDC2_redmapper.DC2/
##			|- cosmoDC2_amico.DC2/

## PREPARE THE OUTPATH DIRECTORY.
outpath += f'{runon}/cosmoDC2_{algo}.{runon}/{cosmoDC2_v}_{algo_runon_v}/'



## PERFORM PROXIMITY MATCHING
if matching_method == 'proximity' :
	mt = ProximityMatch()

	delta_z = 0.05
	match_radius =  1 ## in Mpc
	match_config = {
		'type':'cross',				## OPTIONS: cross, cat1, cat2
		'which_radius':'max',			## OPTIONS: cat1, cat2, min, max
		'preference':'angular_proximity',	## OPTIONS: more_massive, angular_proximity, redshift_proximity
		'catalog1':{'delta_z':delta_z, 'match_radius':f'{match_radius} mpc'},
		'catalog2':{'delta_z':delta_z, 'match_radius':f'{match_radius} mpc'}
		}

	outpath += f'proximity_matching_deltaz_{delta_z}_matchradius_{match_radius}mpc/'

	if os.path.exists(outpath) :
		shutil.rmtree(outpath)
	os.makedirs(outpath)
	print(f'OUTPATH = {outpath}')

	cosmo = AstroPyCosmology()
	mt.match_from_config(cats[0], cats[1], match_config, cosmo=cosmo)
	cats[0].cross_match()
	cats[1].cross_match()

	cats[0].write(outpath + f'{cats[0].name}.fits', overwrite=True)
	cats[1].write(outpath + f'{cats[1].name}.fits', overwrite=True)
	mt.save_matches(cats[0], cats[1], out_dir=outpath, overwrite=True)

	#to print summary
	mt.load_matches(cats[0], cats[1], out_dir=outpath)
	display(cats[0])
	display(cats[1])

## PERFORM MEMBER MATCHING
if matching_method == 'member' :
	mt = MembershipMatch()
	
	minimum_share_fraction = 0.0
	match_config = {
	  	'type':'cross',				## OPTIONS: cross, cat1, cat2
	  	'preference':'shared_member_fraction',	## OPTIONS: shared_member_fraction, more_massive, angular_proximity, redshift_proximity
	  	'minimum_share_fraction':minimum_share_fraction,
	  	'match_members_kwargs': {'method':'id'},
	  	}
	mt.match_from_config(cats[0], cats[1], match_config)
	
	outpath += f'member_matching_fshare_{minimum_share_fraction}/'
	
	if os.path.exists(outpath) :
		shutil.rmtree(outpath)
	os.makedirs(outpath)
	print(f'OUTPATH = {outpath}')
	
	## TO INVESTIGATE MATCHING ISSUES
	mt.fill_shared_members(cats[0], cats[1])
	mt.save_shared_members(cats[0], cats[1], fileprefix=outpath+'mem_share')
	
	cats[0].cross_match()
	cats[1].cross_match()
	
	cats[0].write(outpath + f'{cats[0].name}.fits', overwrite=True)
	cats[1].write(outpath + f'{cats[1].name}.fits', overwrite=True)
	mt.save_matches(cats[0], cats[1], out_dir=outpath, overwrite=True)
	
	## TO PRINT SUMMARY
	mt.load_matches(cats[0], cats[1], out_dir=outpath)
	#display(cats[0])
	#display(cats[1])

sys.exit()
