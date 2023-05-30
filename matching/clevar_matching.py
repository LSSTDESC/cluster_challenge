

import numpy as np
from astropy import units as u
from astropy.io import fits
import sys
import os
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
available_catalogs = ['cosmoDC2', 'wazp_cosmoDC2', 'wazp_DC2', 'redmapper_cosmoDC2', 'redmapper_DC2']
try :
	## STRUCTURE: [[cat_name, version], [cat_name, version]]
	## cat_name OPTIONS: cosmoDC2, wazp_cosmoDC2, wazp_DC2, redmapp_cosmoDC2, redmapper_DC2
	## version OPTIONS: v0, v1, ... (DEPENDING ON CATALOG AVAILABILITY)
	catalogs = np.array([[sys.argv[1], sys.argv[2]], [sys.argv[3], sys.argv[4]]])
except ValueError :
	print(f'Must choose a catalog name: {available_catalogs}. \n \
		See {inpath} for available versions.')




## MATCH BY PROXIMITY OR MEMBER
match_by = sys.argv[5]
if match_by == 'proximity_matching' :
	proximity_matching = True
	member_matching = False
elif match_by == 'member_matching' :
	proximity_matching = False
	member_matching = True
else :
	sys.exit(f'{match_by} is not a valid matching method.')


## COLLECT THE CATALOGS TO MATCH.
cats = []


## COSMODC2
if np.any(catalogs[:,0] == 'cosmoDC2') :
	version = catalogs[:,1][catalogs[:,0] == 'cosmoDC2'][0]
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

	if os.path.exists(inpath + f'cosmoDC2/halos/{version}/') :
		cat = ClCatalog.read(inpath + f'cosmoDC2/halos/{version}/Catalog.fits', f'cosmoDC2_{version}', tags=cltags)
		cat.read_members(inpath + f'cosmoDC2/halos/{version}/Catalog_members.fits', tags=mbtags)
		cats.append(cat)
	else :
		sys.exit(f'The version {version} of {catalog} does not exist. Please see {inpath} for available versions.')
	


## WAZP RUN ON COSMODC2
if np.any(catalogs[:,0] == 'wazp_cosmoDC2') :
	algo   = 'wazp'
	runon  = 'cosmoDC2'
	version = catalogs[:,1][catalogs[:,0] == 'wazp_cosmoDC2'][0]
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

	if os.path.exists(inpath + f'{runon}/{algo}/{version}/') :
		cat = ClCatalog.read(inpath + f'{runon}/{algo}/{version}/Catalog.fits', f'{algo}_{runon}_{version}', tags=cltags)
		cat.read_members(inpath + f'{runon}/{algo}/{version}/Catalog_members.fits', tags=mbtags)
		cats.append(cat)
	else :
		sys.exit(f'The version {version} of {catalog} does not exist. Please see {inpath} for available versions.')



## WAZP RUN ON DC2
if np.any(catalogs[:,0] == 'wazp_DC2') :
	algo   = 'wazp'
	runon  = 'DC2'
	version = catalogs[:,1][catalogs[:,0] == 'wazp_DC2'][0]
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

	if os.path.exists(inpath + f'{runon}/{algo}/{version}/') :
		cat = ClCatalog.read(inpath + f'{runon}/{algo}/{version}/Catalog.fits', f'{algo}_{runon}_{version}', tags=cltags)
		cat.read_members(inpath + f'{runon}/{algo}/{version}/Catalog_members.fits', tags=mbtags)
		cats.append(cat)
	else :
		sys.exit(f'The version {version} of {catalog} does not exist. Please see {inpath} for available versions.')



## REDMAPPER RUN ON COSMODC2
if np.any(catalogs[:,0] == 'redmapper_cosmoDC2') :
	algo   = 'redmapper'
	runon  = 'cosmoDC2'
	version = catalogs[:,1][catalogs[:,0] == 'redmapper_cosmoDC2'][0]
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

	if os.path.exists(inpath + f'{runon}/{algo}/{version}/') :
		cat = ClCatalog.read(inpath + f'{runon}/{algo}/{version}/Catalog.fits', f'{algo}_{runon}_{version}', tags=cltags)
		cat.read_members(inpath + f'{runon}/{algo}/{version}/Catalog_members.fits', tags=mbtags)
		cats.append(cat)
	else :
		sys.exit(f'The version {version} of {catalog} does not exist. Please see {inpath} for available versions.')



## REDMAPPER RUN ON DC2
if np.any(catalogs[:,0] == 'redmapper_DC2') :
	algo   = 'redmapper'
	runon  = 'DC2'
	version = catalogs[:,1][catalogs[:,0] == 'redmapper_DC2'][0]
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

	if os.path.exists(inpath + f'{runon}/{algo}/{version}/') :
		cat = ClCatalog.read(inpath + f'{runon}/{algo}/{version}/Catalog.fits', f'{algo}_{runon}_{version}', tags=cltags)
		cat.read_members(inpath + f'{runon}/{algo}/{version}/Catalog_members.fits', tags=mbtags)
		cats.append(cat)
	else :
		sys.exit(f'The version {version} of {catalog} does not exist. Please see {inpath} for available versions.')



## OUTPATH DIRECTORY STRUCTURE: (ASSUMING THERE ARE NOT COSMODC2 - DC2 COMPARISONS)
##	after_matching/
##		|- cosmoDC2/
##			|- cosmoDC2_clfinder/
##				|- cosmoDC2_wazp.cosmoDC2/
##					|- v0_v0/
##					|- v1_v0/
##				|- cosmoDC2_redmapper.cosmoDC2/
##			|- clfinder_clfinder/
##				|- redmapper.cosmoDC2_wazp.cosmoDC2/
##					|- v0_v0/
##				|- amico.cosmoDC2_wazp.cosmoDC2/
##				|- amico.cosmoDC2_redmapper.cosmoDC2/
##		|- DC2/
##			|- clfinder_clfinder/
##				|- redmapper.DC2_wazp.DC2/
##					|- v0_v0/
##				|- amico.DC2_wazp.DC2/
##				|- amico.DC2_redmapper.DC2/

## PREPARE THE OUTPATH DIRECTORY.
if np.any(['cosmoDC2' in c for c in catalogs[:,0]]) :
	outpath += 'cosmoDC2/'
	if np.any(catalogs[:,0] == 'cosmoDC2') :
		outpath += 'cosmoDC2_clfinder/'
		algo = catalogs[:,0][catalogs[:,0] != 'cosmoDC2'][0].split('_', 1)[0]
		vA = catalogs[:,1][catalogs[:,0] == 'cosmoDC2'][0]
		vB = catalogs[:,1][catalogs[:,0] != 'cosmoDC2'][0]
		outpath += f'cosmoDC2_{algo}.cosmoDC2/{vA}_{vB}/'
	else :
		outpath += 'clfinder_clfinder/'
		clfinders = sorted([c.split('_',1)[0] for c in catalogs[:,0]])
		vA = catalogs[:,1][[clfinders[0] in c for c in catalogs[:,0]]][0]
		vB = catalogs[:,1][[clfinders[1] in c for c in catalogs[:,0]]][0]
		outpath += f'{clfinders[0]}.cosmoDC2_{clfinders[1]}.cosmoDC2/{vA}_{vB}/'
else :
	outpath += 'DC2/clfinder_clfinder/'
	clfinders = sorted([c.split('_',1)[0] for c in catalogs[:,0]])
	vA = catalogs[:,1][[clfinders[0] in c for c in catalogs[:,0]]][0]
	vB = catalogs[:,1][[clfinders[1] in c for c in catalogs[:,0]]][0]
	outpath += f'{clfinders[0]}.DC2_{clfinders[1]}.DC2/{vA}_{vB}/'
	



## PERFORM PROXIMITY MATCHING
if proximity_matching :
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
if member_matching :
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
