
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


algo = 'wazp'
runon = 'cosmoDC2'


## KEEP TRACK OF VERSIONS (FEEL FREE TO ADD additional_comments TO THE VERSIONS IF NEEDED)
versions = [
        {'v':'v0',
	'cat_name':'cosmoDC2_v1.1.4_wazp_v1.0_truez',
	'min_richness':0,
	'description':'WaZP run on cosmoDC2 with true redshifts as input',},
        {'v':'v1',
	'cat_name':'cosmoDC2_v1.1.4_wazp_v1.0_truez',
	'min_richness':20,
	'description':'WaZP run on cosmoDC2 with true redshifts as input',},
	]


## SET THE VERSION TO WORK WITH (DEFAULT: v0)
try :
	index = np.argwhere([versions[i]['v'] == str(sys.argv[1]) for i in range(len(versions))])[0][0]
	version = versions[index]
except :
	version = versions[0]
print(f"Producing reduced {algo}.{runon}catalog:  {version['v']} \t {version['cat_name']}")


## OUTPATH
outpath = f"/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/{runon}/{algo}/{version['v']}/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print(f'outpath = {outpath}')

## GET DATA
## POSSIBLE FILTERS: min_richness, min_z_cl, max_z_cl
cl_data, mb_data = wazp_cosmoDC2_cat_open(version['cat_name'], version['min_richness'],)

## MAKE TABLES
cl_table = Table([
        cl_data['cluster_id'],		## CLUSTER ID: id_cl
        cl_data['cluster_ra'],		## CLUSTER RA: ra_cl
        cl_data['cluster_dec'],		## CLUSTER DEC: dec_cl
        cl_data['cluster_z'],		## CLUSTER REDSHIFT: z_cl
        cl_data['cluster_ngals']],	## CLUSTER RICHNESS: mass (FOR SORTING PURPOSES)
        names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass'))

mb_table = Table([
        mb_data['member_id'],		## MEMBER ID: id_mb
        mb_data['member_id_cluster'],	## CORRESPONDING CLUSTER ID: clid_mb
        mb_data['member_ra'],		## MEMBER RA: ra_mb
        mb_data['member_dec'],		## MEMBER DEC: dec_mb
        mb_data['member_z'],		## MEMBER REDSHIFT: z_mb
	mb_data['member_pmem']],	## MEMBER PMEM: pmem
        names=('id_mb', 'clid_mb', 'ra_mb', 'dec_mb', 'z_mb', 'pmem'))

### COMPUTE RICHNESS
#mb_sorted_by_cl = np.argsort(mb_data['member_id_cluster'])[::-1]
#mb_pmem_sorted = mb_data['member_pmem'][mb_sorted_by_cl]
#indices = np.insert(np.cumsum(cl_data['cluster_nmem']), 0, 0)
#richness = np.array([np.sum(mb_pmem_sorted[indices[i]:indices[i+1]]) for i in range(len(indices)-1)])
#cl_table['mass'] = richness

## WRITE TABLES TO FILES
cl_table.write(outpath + 'Catalog.fits', overwrite=True)
mb_table.write(outpath + 'Catalog_members.fits', overwrite=True)


## ADD CURRENT ITERATION OF VERSION DOCUMENTATION TO versions.txt
from datetime import datetime
with open(f'/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/{runon}/{algo}/versions.txt', 'a+') as vf :
	vf.write(f"\n{datetime.now()}\twriting: {version['v']}\n\t")
	for i in range(len(versions)) :
		keys = list(versions[i].keys())
		details = [f"{key}: {versions[i][key]}" for key in keys]
		vf.write('\t'.join(details) + '\n\t')


sys.exit()

