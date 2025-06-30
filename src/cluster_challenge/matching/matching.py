

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import os, sys, yaml
import json
import shutil
import healpy as hp

###clevar
import clevar
print(clevar.version.__version__)
from clevar.catalog import ClCatalog
from clevar.match import output_catalog_with_matching
from clevar.match import output_matched_catalog
from clevar.cosmology import AstroPyCosmology


config = sys.argv[1]
with open(config) as fstream :
	cfg = yaml.safe_load(fstream)


mt_method = cfg['matching']['cl_method']
mt_type = cfg['matching']['type']
mt_preference = cfg['matching']['pref']

inpath  = cfg['paths']['matching']['in']
outpath = cfg['paths']['matching']['out']


def mt_footprint(config, cats) :
    src = __file__.replace('matching.py', '')

    ftp = [np.load(f"{src}../footprints/{config['footprints']['cat1']['pixFile']}"),
            np.load(f"{src}../footprints/{config['footprints']['cat2']['pixFile']}")]

    cat_hpix = [
            hp.ang2pix(
                config['footprints']['cat2']['NSIDE'],
                cats[0][config['cat_keys']['cat1']['cluster']['ra']],
                cats[0][config['cat_keys']['cat1']['cluster']['dec']],
                nest=False,
                lonlat=True),
            hp.ang2pix(
                config['footprints']['cat1']['NSIDE'],
                cats[1][config['cat_keys']['cat2']['cluster']['ra']],
                cats[1][config['cat_keys']['cat2']['cluster']['dec']],
                nest=False,
                lonlat=True)]
    mask = [np.isin(cat_hpix[0], ftp[1]),
            np.isin(cat_hpix[1], ftp[0])]

    return mask



## COLLECT THE CATALOGS TO MATCH.
def read_members(cat, cfg, cat12) :
    mb_tags = cfg['cat_keys'][cat12]['member']
    mb_cat = Table.read(inpath[cat12]['member'])
    cat.add_members(
            id = mb_cat[mb_tags['id']].astype(str),
            id_cluster = mb_cat[mb_tags['id_cluster']].astype(type(cat['id'][0])),
            ra = mb_cat[mb_tags['ra']],
            dec = mb_cat[mb_tags['dec']],
            z = mb_cat[mb_tags['z']],
            mag_u = mb_cat[mb_tags['mag_u']],
            mag_g = mb_cat[mb_tags['mag_g']],
            mag_r = mb_cat[mb_tags['mag_r']],
            mag_i = mb_cat[mb_tags['mag_i']],
            mag_z = mb_cat[mb_tags['mag_z']],
            mag_y = mb_cat[mb_tags['mag_y']]
            )
    del mb_cat


cats = []
for c in ['cat1', 'cat2'] :
    print(f"loading from\n\t{inpath[c]['cluster']}...")
    if os.path.exists(inpath[c]['cluster']) & os.path.exists(inpath[c]['member']) :
        cltags = cfg['cat_keys'][c]['cluster']
        mbtags = cfg['cat_keys'][c]['member']
        
        cat = ClCatalog.read(inpath[c]['cluster'], c, tags=cltags, full=False)
        if ('snr' in cfg['matching'].keys()) :
            if (c in cfg['matching']['snr'].keys()) :
                cat = cat[cat[cfg['cat_keys'][c]['cluster']['snr']] > cfg['matching']['snr'][c]]

        read_members(cat, cfg, c)
        if not (cat.leftover_members is None) :
            print(f"{cfg['cats'][c]} has {len(cat.leftover_members['id']) / (len(cat.leftover_members['id']) + len(cat.members['id'])) * 100:.1f}% of members in leftover_members")
        else :
            print(f"{cfg['cats'][c]} has no leftover_members and {len(cat.members['id'])} members")
        cats.append(cat)
    else :
    	sys.exit(f"The {cfg['cats'][c]} catalog does not exist.")
    print(f"{cfg['cats'][c]} catalog is loaded...")



## APPLY FOOTPRINT MATCHING
print("matching footprints...")
mask = mt_footprint(cfg, cats)
for i in range(2) :
    print(f"\tcat{i+1}", end='...')
    cats[i] = cats[i][mask[i]]
    print('DONE')


## PERFORM PROXIMITY MATCHING
if mt_method == 'proximity' :
	from clevar.match import ProximityMatch
	mt = ProximityMatch()
	
	delta_z = float(cfg['matching']['delta_z'])
	match_radius = float(cfg['matching']['match_radius'])	## in Mpc
	
	outpath += f'proximity_matching/deltaz_{delta_z}_matchradius_{match_radius}mpc_pref_{mt_preference}/'

	match_config = {
		'type':mt_type,			## OPTIONS: cross, cat1, cat2
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
    
    minimum_share_fraction1 = float(cfg['matching']['minimum_share_fraction1'])
    minimum_share_fraction2 = float(cfg['matching']['minimum_share_fraction2'])
    delta_z = float(cfg['matching']['delta_z'])
    if 'mb_method' in cfg['matching'].keys() :
        mt_method_mb = cfg['matching']['mb_method']
    else :
        mt_method_mb = 'id'
    
    outpath += f'member_matching/fshare_{minimum_share_fraction2}_pref_{mt_preference}_mbmt_{mt_method_mb}/'
    
    match_config = {
      	'type':mt_type,			## OPTIONS: cross, cat1, cat2
      	'preference':mt_preference,	## OPTIONS: shared_member_fraction, more_massive, angular_proximity, redshift_proximity
      	'minimum_share_fraction1':minimum_share_fraction1,
      	'minimum_share_fraction2':minimum_share_fraction2,
        'match_members_kwargs': {'method':mt_method_mb, 'delta_z':delta_z},
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
    
    for cat in cats :
    	cat.members['mt_multi_self'] = cat.members['match']
    	del cat.members['match']
    
    cats[0].members.write(f"{outpath}{cats[0].name}_mbs.fits", overwrite=True)
    cats[1].members.write(f"{outpath}{cats[1].name}_mbs.fits", overwrite=True)

sys.exit()
