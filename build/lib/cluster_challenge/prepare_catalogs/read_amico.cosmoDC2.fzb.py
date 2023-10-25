
## THIS SCRIPT MAKES A REDUCED TABLE FOR AMICO RUN ON COSMODC2 WITH FZB
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


algo = 'amico'
runon = 'cosmoDC2.fzb'


## KEEP TRACK OF VERSIONS (FEEL FREE TO ADD additional_comments TO THE VERSIONS IF NEEDED)
versions = [
        {'v':'v0',
	'cat_name':'amico.cosmoDC2.fzb',
	'min_richness':0,
	'description':'AMICO run on cosmoDC2 with redshifts provide by FlexZBoost photoz algorithm (using zmodes)',},
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
        f"/sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/{algo}/{runon}/{version['v']}/Catalog.fits")[1].data)
mb_data = Table(fits.open(
        f"/sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/{algo}/{runon}/{version['v']}/Catalog_members.fits")[1].data)



## MAKE TABLES
cl_table = Table([
        cl_data['id'],		## CLUSTER ID: id_cl
        cl_data['ra'],		## CLUSTER RA: ra_cl
        cl_data['dec'],		## CLUSTER DEC: dec_cl
        cl_data['z'],		## CLUSTER REDSHIFT: z_cl
        cl_data['mass']],	## CLUSTER MASS: mass
        names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass'))

mb_table = Table([
        mb_data['id'],		## MEMBER ID: id_mb
        mb_data['id_cluster'],	## CORRESPONDING CLUSTER ID: clid_mb
        mb_data['ra'],		## MEMBER RA: ra_mb
        mb_data['dec'],		## MEMBER DEC: dec_mb
        mb_data['z'],		## MEMBER REDSHIFT: z_mb
	mb_data['pmem']],	## MEMBER PMEM: pmem
        names=('id_mb', 'clid_mb', 'ra_mb', 'dec_mb', 'z_mb', 'pmem'))

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

