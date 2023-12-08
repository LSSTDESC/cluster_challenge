
## THIS SCRIPT MAKES A REDUCED TABLE FOR WAZP RUN ON COSMODC2 FROM THE GCR CATALOG
## AND SAVES IT IN THE before_matching SUBDIRECTORY.


## IMPORTS
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
import sys
import os
import shutil

from opening_catalogs_functions import *


algo = 'pywazp'
runon = 'cosmoDC2.truez'


## KEEP TRACK OF VERSIONS (FEEL FREE TO ADD additional_comments TO THE VERSIONS IF NEEDED)
versions = [
        {'v':'v0',
	'cat_name':'pywazp.cosmoDC2.truez',
	'min_richness':0,
	'description':'WaZP-Python run on cosmoDC2 with true redshifts',},
	]


## SET THE VERSION TO WORK WITH (DEFAULT: v0)
try :
	index = np.argwhere([versions[i]['v'] == str(sys.argv[1]) for i in range(len(versions))])[0][0]
	version = versions[index]
except :
	version = versions[0]
print(f"Producing reduced {algo}.{runon} catalog:  {version['v']} \t {version['cat_name']}")


## OUTPATH
outpath = f"/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/{algo}/{runon}/{version['v']}/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print(f'outpath = {outpath}')

## GET DATA
cl_data = Table(fits.open(
        f"/sps/lsst/users/rsolomon/wazp/catalogs/cosmoDC2_truez/wazp_clusters.fits")[1].data)
mb_data = Table(fits.open(
        f"/sps/lsst/users/rsolomon/wazp/catalogs/cosmoDC2_truez/wazp_members.fits")[1].data)


## MAKE TABLES
cl_table = Table([
        cl_data['id'],          ## CLUSTER ID: id_cl
        cl_data['ra'],          ## CLUSTER RA: ra_cl
        cl_data['dec'],         ## CLUSTER DEC: dec_cl
        cl_data['zp'],          ## CLUSTER REDSHIFT: z_cl
        cl_data['n200'],        ## CLUSTER RICHNESS: mass (FOR SORTING PURPOSES)
        cl_data['n200_err'],    ## CLUSTER RICHNESS ERROR: mass_err
	cl_data['snr']],	## CLUSTER SNR: snr_cl
        names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass', 'mass_err', 'snr_cl'))

mb_table = Table([
        mb_data['id_gal'],	## MEMBER ID: id_mb
        mb_data['id_cl'],	## CORRESPONDING CLUSTER ID: clid_mb
        mb_data['ra'],          ## MEMBER RA: ra_mb
        mb_data['dec'],         ## MEMBER DEC: dec_mb
        mb_data['zp'],          ## MEMBER REDSHIFT: z_mb
        mb_data['pmem']],        ## MEMBER PMEM: pmem
        names=('id_mb', 'clid_mb', 'ra_mb', 'dec_mb', 'z_mb', 'pmem',))

## APPLY min_richness
cl_table = cl_table[cl_table['mass'] > version['min_richness']]
mb_table = mb_table[np.isin(mb_table['clid_mb'], cl_table['id_cl'])]


## WRITE TABLES TO FILES
cl_table.write(outpath + 'Catalog.fits', overwrite=True)
mb_table.write(outpath + 'Catalog_members.fits', overwrite=True)


## ADD CURRENT ITERATION OF VERSION DOCUMENTATION TO versions.txt
from datetime import datetime
with open(f'/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/{algo}/{runon}/versions.txt', 'a+') as vf :
	vf.write(f"\n{datetime.now()}\twriting: {version['v']}\n\t")
	for i in range(len(versions)) :
		keys = list(versions[i].keys())
		details = [f"{key}: {versions[i][key]}" for key in keys]
		vf.write('\t'.join(details) + '\n\t')


sys.exit()

