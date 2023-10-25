
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
runon = 'DC2'

## KEEP TRACK OF VERSIONS (just add new versions, do not delete previous versions from dict)
versions = [
	{'v':'v0',
	'cat_name':'dc2_redmapper_run2.2i_dr6_wfd_v0.8.1',
	'min_richness':0,
	'description':'redMaPPer run on DC2 galaxy cluster catalog',},
	{'v':'v1',
	'cat_name':'dc2_redmapper_run2.2i_dr6_wfd_v0.8.1',
	'min_richness':20,
	'description':'redMaPPer run on DC2 galaxy cluster catalog',},
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

if os.path.exists(outpath) :
     shutil.rmtree(outpath)
os.makedirs(outpath)
print(f'outpath = {outpath}')

## GET DATA
## POSSIBLE FILTERS: min_richness, min_z_cl, max_z_cl
cl_data, mb_data = redmapper_cat_open(
	cat_name=version['cat_name'],
	min_richness=version['min_richness'],)


## MAKE TABLES
cl_table = Table([
        cl_data['cluster_id'],          ## CLUSTER ID: id_cl
        cl_data['ra'],          	## CLUSTER RA: ra_cl
        cl_data['dec'],         	## CLUSTER DEC: dec_cl
        cl_data['redshift'],		## CLUSTER REDSHIFT: z_cl
        cl_data['richness'],		## CLUSTER RICHNESS: mass (FOR SORTING PURPOSES)
        cl_data['richness_err']],	## CLUSTER RICHNESS ERROR: mass_err
        names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass', 'mass_err'))
cl_table['snr_cl'] = cl_data['richness'] / cl_data['richness_err']

mb_table = Table([
        mb_data['id_member'],           ## MEMBER ID: id_mb
        mb_data['cluster_id_member'],   ## CORRESPONDING CLUSTER ID: clid_mb
        mb_data['ra_member'],           ## MEMBER RA: ra_mb
        mb_data['dec_member'],          ## MEMBER DEC: dec_mb
        mb_data['redshift_true_member'],## MEMBER REDSHIFT: z_mb
        mb_data['p_member'],		## MEMBER PMEM: pmem
	mb_data['mag_g_lsst_member'],	## MEMBER MAGNITUDE g BAND: mag_g
	mb_data['mag_r_lsst_member'],	## MEMBER MAGNITUDE r BAND: mag_r
	mb_data['mag_i_lsst_member'],	## MEMBER MAGNITUDE i BAND: mag_i
	mb_data['mag_z_lsst_member'],	## MEMBER MAGNITUDE z BAND: mag_z
	mb_data['mag_y_lsst_member']],	## MEMBER MAGNITUDE y BAND: mag_y
        names=('id_mb', 'clid_mb', 'ra_mb', 'dec_mb', 'z_mb', 'pmem', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y'))

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

