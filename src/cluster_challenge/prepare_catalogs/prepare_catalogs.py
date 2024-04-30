
## IMPORTS
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
import sys, os, yaml, shutil

from src.opening_catalogs_functions import *

with open(sys.argv[1]) as fstream :
	cfg = yaml.safe_load(fstream)


print(f"Producing reduced {cfg['name']} catalog")


## outpath
outpath = cfg['outpath']

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print(f'outpath = {outpath}')

## GET DATA
cl_table, mb_table = get_cl_and_mb_tables(cfg)

## EDIT SELECT COLUMN NAMES
[cl_table.rename_column(old, new) for (old, new) in zip(cfg['change_cols']['cluster']['old'], cfg['change_cols']['cluster']['new'])]
[mb_table.rename_column(old, new) for (old, new) in zip(cfg['change_cols']['member']['old'], cfg['change_cols']['member']['new'])]


### HALOS CATALOG DOES NOT HAVE PMEM VALUES
#if cfg['algo'] == 'halos' :
#	mb_table['pmem'] = np.ones_like(mb_table['z_mb'])
if 'pmem' not in mb_table.colnames :
	mb_table['pmem'] = np.ones_like(mb_table['z_mb'])


## APPLY CUTS
#cl_table['n200-n500kpc'] = cl_table['mass'] - cl_table['n500kpc']
if 'cuts' in cfg.keys() :
	for key in cfg['cuts'].keys() :
		if key == 'n200-n500kpc' :
			cut = np.array((cl_table['mass'] - cl_table['n500kpc']) > 0)
		else :
			cut = np.array([eval(cfg['cuts'][key]) for _ in cl_table[key]])
		cl_table = cl_table[cut]
	mb_table = mb_table[np.isin(mb_table['clid_mb'], cl_table['id_cl'])]

#if cfg['algo'] != 'halos' :
#	## REMOVE NON-CLUSTER MEMBERS
#	mb_table = mb_table[np.isin(mb_table['clid_mb'], cl_table['id_cl'])]

#if cfg['algo'] == 'redmapper' :	## REDMAPPER DOES NOT HAVE SNR COMPUTED FOR THE CLUSTERS
#	cl_table['snr_cl'] = cl_table['mass'] / cl_table['mass_err']


## WRITE TABLES TO FILES
cl_table.write(outpath + 'Catalog.fits', overwrite=True)
mb_table.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()

