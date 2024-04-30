

import numpy as np
from astropy import units as u
from astropy.io import fits
import os, sys, yaml
import json
import shutil

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import output_catalog_with_matching
from clevar.match import output_matched_catalog
from clevar.cosmology import AstroPyCosmology
from IPython.display import display


config = sys.argv[1]
with open(config) as fstream :
	cfg = yaml.safe_load(fstream)


#requested_cats = [cfg['cats']['cat1'].split('.'), cfg['cats']['cat2'].split('.')]
mt_method = cfg['matching']['method']
mt_preference = cfg['matching']['pref']

inpath  = cfg['paths']['matching']['in']
outpath = cfg['paths']['matching']['out']



## COLLECT THE CATALOGS TO MATCH.
cats = []
cat_locs = []

#for c in requested_cats :
for c in ['cat1', 'cat2'] :
	print(f"loading from {inpath[c]}...")
	if os.path.exists(inpath[c]) :
		cltags = cfg['cat_keys'][c]['cluster']
		mbtags = cfg['cat_keys'][c]['member']

		cat_locs.append({inpath[c]})
		cat = ClCatalog.read(os.path.join(inpath[c], 'Catalog.fits'), cfg['cats'][c], tags=cltags, full=True)
		cat.read_members(os.path.join(inpath[c], 'Catalog_members.fits'), tags=mbtags, full=True)
		cats.append(cat)
	else :
		sys.exit(f"The {cfg['cats'][c]} catalog does not exist.")
	print(f"{cfg['cats'][c]} catalog is loaded...")
	


### apply cuts
#if not (cfg['matching']['cuts'] is None) :
#	for key in cfg['matching']['cuts']['cat1'].keys() :
		

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


#outpath += f"{'.'.join(requested_cats[0][:-1])}_{'.'.join(requested_cats[1][:-1])}/"
#outpath += f"{requested_cats[0][-1]}_{requested_cats[1][-1]}/"



## PERFORM PROXIMITY MATCHING
if mt_method == 'proximity' :
	from clevar.match import ProximityMatch
	mt = ProximityMatch()
	
	delta_z = float(cfg['matching']['delta_z'])
	match_radius = float(cfg['matching']['match_radius'])	## in Mpc
	
	outpath += f'proximity_matching/deltaz_{delta_z}_matchradius_{match_radius}mpc_pref_{mt_preference}/'

	match_config = {
		'type':'cross',			## OPTIONS: cross, cat1, cat2
		'which_radius':'max',		## OPTIONS: cat1, cat2, min, max
		'preference':mt_preference,	## OPTIONS: more_massive, angular_proximity, redshift_proximity
		'catalog1':{'delta_z':delta_z, 'match_radius':f'{match_radius} mpc'},
		'catalog2':{'delta_z':delta_z, 'match_radius':f'{match_radius} mpc'},
		}


	if os.path.exists(outpath) :
		shutil.rmtree(outpath)
	os.makedirs(outpath)
	print(f'OUTPATH = {outpath}')

	#cats_raw = [cats[0].raw(), cats[1].raw()]
	cosmo = AstroPyCosmology()

	print('MATCHING...')
	mt.match_from_config(cats[0], cats[1], match_config, cosmo=cosmo)

	print('WRITING...')
	cats[0].write(f"{outpath}{cats[0].name}.fits", overwrite=True)
	cats[1].write(f"{outpath}{cats[1].name}.fits", overwrite=True)




## PERFORM MEMBER MATCHING
if mt_method == 'member' :
	from clevar.match import MembershipMatch
	mt = MembershipMatch()
	
	minimum_share_fraction = float(cfg['matching']['minimum_share_fraction'])

	outpath += f'member_matching/fshare_{minimum_share_fraction}_pref_{mt_preference}/'

	match_config = {
	  	'type':'cross',			## OPTIONS: cross, cat1, cat2
	  	'preference':mt_preference,	## OPTIONS: shared_member_fraction, more_massive, angular_proximity, redshift_proximity
	  	'minimum_share_fraction':minimum_share_fraction,
	  	'match_members_kwargs': {'method':'id'},
	  	}
	
	
	if os.path.exists(outpath) :
		shutil.rmtree(outpath)
	os.makedirs(outpath)
	print(f'OUTPATH = {outpath}')
	
	print('MATCHING...')
	mt.match_from_config(cats[0], cats[1], match_config)
	
	print('WRITING...')
	cats[0].write(f"{outpath}{cats[0].name}.fits", overwrite=True)
	cats[1].write(f"{outpath}{cats[1].name}.fits", overwrite=True)



sys.exit()
