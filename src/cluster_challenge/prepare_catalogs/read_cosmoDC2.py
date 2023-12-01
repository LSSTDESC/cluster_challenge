
## THIS SCRIPT MAKES A REDUCED TABLE FOR COSMODC2 FROM THE GCR CATALOG
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

np.bool = bool	## REQUIRED DUE TO DEPRECATION OF np.bool


## KEEP TRACK OF VERSIONS
versions = [
	{'v':'v0',
        'cat_name':'halos_m200c_13.0',
        'min_mass':1e13,
        'min_mass_units':'Msun',
        'description':'the cosmoDC2_v1.1.4 galaxy cluster catalog cluster of mass > 1e13Msun',},
	]


## SET THE VERSION TO WORK WITH (DEFAULT: v0)
try :
        index = np.argwhere([versions[i]['v'] == str(sys.argv[1]) for i in range(len(versions))])[0][0]
        version = versions[index]
except :
        version = versions[0]
print(f"Producing reduced cosmoDC2 catalog:  {version['v']} \t {version['cat_name']}")


## OUTPATH
outpath = f"/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/halos/cosmoDC2/{version['v']}/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print(f'outpath = {outpath}')

## MASS CUT
#truth_data = DC2_cat_open(version['cat_name'], version['min_mass'], cluster_only=False)

## MAKE TABLES
#cl_data = truth_data[truth_data['is_central']==True]
#mb_data = truth_data[truth_data['is_central']==False]
cl_data = Table(fits.open(
	f"/sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/halos/cosmoDC2/halos_m200c_13.0.fits")[1].data)
mb_data = Table(fits.open(
	f"/sps/lsst/groups/clusters/cluster_comparison_project/initial_catalogs/halos/cosmoDC2/halos_m200c_13.0_members.fits")[1].data)


cl_table = Table([
	cl_data['halo_id'],		## CLUSTER ID: id_cl
	cl_data['ra'],			## CLUSTER RA: ra_cl
	cl_data['dec'],			## CLUSTER DEC: dec_cl
	cl_data['redshift_true'],	## CLUSTER REDSHIFT: z_cl
	cl_data['m200c']],		## HALO MASS: mass
	names=('id_cl', 'ra_cl', 'dec_cl', 'z_cl', 'mass'))
cl_table['log10_mass'] = np.log10(cl_table['mass'])

mb_table = Table([
	mb_data['galaxyID'],		## MEMBER ID: id_mb
	mb_data['halo_id'],		## CORRESPONDING CLUSTER ID: clid_mb
	mb_data['ra'],			## MEMBER RA: ra_mb
	mb_data['dec'],			## MEMBER DEC: dec_mb
	mb_data['redshift_true'],	## MEMBER REDSHIFT: z_mb
	mb_data['mag_true_g'],		## MEMBER MAGNITUDE IN g BAND: mag_g
	mb_data['mag_true_r'],		## MEMBER MAGNITUDE IN r BAND: mag_r
	mb_data['mag_true_i'],		## MEMBER MAGNITUDE IN i BAND: mag_i
	mb_data['mag_true_z'],		## MEMBER MAGNITUDE IN z BAND: mag_z
	mb_data['mag_true_y']],		## MEMBER MAGNITUDE IN y BAND: mag_y
	names=('id_mb', 'clid_mb', 'ra_mb', 'dec_mb', 'z_mb', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y'))
mb_table['pmem'] = np.ones_like(mb_table['z_mb'])

## WRITE TABLES TO FILES
cl_table.write(outpath + 'Catalog.fits', overwrite=True)
mb_table.write(outpath + 'Catalog_members.fits', overwrite=True)


## ADD CURRENT ITERATION OF VERSION DOCUMENTATION TO versions.txt
from datetime import datetime
with open('/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/halos/cosmoDC2/versions.txt', 'a+') as vf :
        vf.write(f"\n{datetime.now()}\twriting: {version['v']}\n\t")
        for i in range(len(versions)) :
                keys = list(versions[i].keys())
                details = [f"{key}: {versions[i][key]}" for key in keys]
                vf.write('\t'.join(details) + '\n\t')


sys.exit()

