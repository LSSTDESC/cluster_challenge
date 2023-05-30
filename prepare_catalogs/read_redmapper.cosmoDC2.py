
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


algo  = 'redmapper'
runon = 'cosmoDC2'

## KEEP TRACK OF VERSIONS (just add new versions, do not delete previous versions from dict)
versions = [
	{'v':'v0',
	'cat_name':'cosmoDC2_v1.1.4_redmapper_v0.8.1',
	'min_richness':20,
	'description':'redMaPPer run on cosmoDC2 galaxy cluster catalog',}
	,]


## SET THE VERSION TO WORK WITH (DEFAULT: v0)
try :
        index = np.argwhere([versions[i]['v'] == str(sys.argv[1]) for i in range(len(versions))])[0][0]
        version = versions[index]
except :
        version = versions[0]
print(f"Producing reduced {algo}.{runon} catalog:  {version['v']} \t {version['cat_name']}")


## OUTPATH
outpath = f"/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/{runon}/{algo}/{version['v']}/"

if os.path.exists(outpath) :
     shutil.rmtree(outpath)
os.makedirs(outpath)
print(f'outpath = {outpath}')

## GET DATA
## POSSIBLE FILTERS: min_richness, min_z_cl, max_z_cl
cl_data, mb_data = redmapper_cat_open(
	cat_name=version['cat_name'],
	min_richness=version['min_richness'],)

## REDMAPPER GCR CATALOG DOES NOT HAVE MEMBER REDSHIFTS OR CLUSTER MASS
## SO WE GET THIS FROM ANOTHER CATALOG
zs = Table(fits.open('/sps/lsst/groups/desc/shared/xgal/cosmoDC2/addons/redmapper_v1.1.4/cosmoDC2_v1.1.4_mag_run_redmapper_v0.8.1_lgt05_catalog_members.fit')[1].data)
zs.keep_columns(['id','zred'])

if not np.all(mb_data['id_member'] == zs['id']) :
	sys.exit(f'The member id does not match those in the xgal catalog from which the member redshifts are obtained. Please check.')


## MAKE TABLES
cl_table = Table([
        cl_data['cluster_id'],          ## CLUSTER ID: id_cl
        cl_data['ra'],          	## CLUSTER RA: ra_cl
        cl_data['dec'],         	## CLUSTER DEC: dec_cl
        cl_data['redshift'],           	## CLUSTER REDSHIFT: z_cl
        cl_data['richness']],      	## CLUSTER RICHNESS: mass (FOR SORTING PURPOSES)
        names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass'))

mb_table = Table([
        mb_data['id_member'],           ## MEMBER ID: id_mb
        mb_data['cluster_id_member'],   ## CORRESPONDING CLUSTER ID: clid_mb
        mb_data['ra_member'],           ## MEMBER RA: ra_mb
        mb_data['dec_member'],          ## MEMBER DEC: dec_mb
	zs['zred'],			## MEMBER REDSHIFT: z_mb
        mb_data['p_member']],	        ## MEMBER PMEM: pmem
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

