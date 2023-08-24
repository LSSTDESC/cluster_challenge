## THIS SCRIPT MAKES A REDUCED TABLE FOR WAZP RUN ON DC2 FROM FILES STORED AT
## /sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/wazp/DC2/
## AND SAVES IT IN THE before_matching SUBDIRECTORY.


## IMPORTS
import numpy as np
from astropy.table import Table
from astropy.io import fits
import sys
import os
import shutil

from opening_catalogs_functions import *


algo = 'wazp'
runon = 'DC2'


## KEEP TRACK OF VERSIONS (FEEL FREE TO ADD additional_comments TO THE VERSIONS IF NEEDED)
versions = [
        {'v':'v0',
        'cat_name':'7081',
	'min_richness':0,
        'description':'WaZP run on DC2 from the LineA 7081 run',},
        {'v':'v1',
        'cat_name':'7081',
	'min_richness':20,
        'description':'WaZP run on DC2 from the LineA 7081 run',},
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
cl_data = Table(fits.open(
	f"/sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/wazp/DC2/{version['cat_name']}/wazp_cluster.fits")[1].data)
mb_data = Table(fits.open(
	f"/sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/wazp/DC2/{version['cat_name']}/wazp_membership.fits")[1].data)


## MAKE TABLES
cl_table = Table([
        cl_data['ID'],		## CLUSTER ID: id_cl
        cl_data['RA'],          ## CLUSTER RA: ra_cl
        cl_data['DEC'],         ## CLUSTER DEC: dec_cl
        cl_data['zp'],          ## CLUSTER REDSHIFT: z_cl
        cl_data['NGALS']],      ## CLUSTER RICHNESS: mass (FOR SORTING PURPOSES)
        names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass'))

mb_table = Table([
        mb_data['ID_g'],        ## MEMBER ID: id_mb
        mb_data['ID_CLUSTER'],	## CORRESPONDING CLUSTER ID: clid_mb
        mb_data['RA'],          ## MEMBER RA: ra_mb
        mb_data['DEC'],         ## MEMBER DEC: dec_mb
        mb_data['ZP'],          ## MEMBER REDSHIFT: z_mb
        mb_data['PMEM']],       ## MEMBER PMEM: pmem
        names=('id_mb', 'clid_mb', 'ra_mb', 'dec_mb', 'z_mb', 'pmem'))

### COMPUTE RICHNESS
#mb_sorted_by_cl = np.argsort(mb_data['ID_CLUSTER'])[::-1]
#mb_pmem_sorted = mb_data['PMEM'][mb_sorted_by_cl]
#indices = np.insert(np.cumsum(cl_data['NMEM']), 0, 0)
#richness = np.array([np.sum(mb_pmem_sorted[indices[i]:indices[i+1]]) for i in range(len(indices)-1)])
#cl_table['mass'] = richness


## APPLY min_richness
cl_table = cl_table[cl_table['mass'] > version['min_richness']]
mb_table = mb_table[np.isin(mb_table['clid_mb'], cl_table['id_cl'])]

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

