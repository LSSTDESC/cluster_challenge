
import yaml, sys
import numpy as np

from clevar.catalog import ClCatalog

import src.saving_figures as sfigs
import fit_func as ff
from src.pltFunc_redshift import *


## read config files as online arguments 
config = sys.argv[1]
with open(config) as fstream:
	cfg = yaml.safe_load(fstream)


## get the file path for the matched catalogs as determined in config file
cat1 = cfg['cats']['cat1']
cat2 = cfg['cats']['cat2']
mt_method = cfg['matching']['method']
mt_pref = cfg['matching']['pref']
if mt_method == 'member' :
        mt_params = [cfg['matching']['minimum_share_fraction']]
else :
        mt_params = [cfg['matching']['delta_z'], cfg['matching']['match_radius']]

#inpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['in'])
#outpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['out'], addon='redshift')
inpath = ossfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['in'])
outpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['out'], addon='redshift')



## open the merged catalog
cltags = cfg['merged_keys']
mtcats = ClCatalog.read(f"{inpath}merged_catalog.fits", 'merged', tags=cltags)



## redshift vs redshift
x = mtcats['z_halo']
y = mtcats['z_cl']
xlabel = cfg['latex']['redshift_1']
ylabel = cfg['latex']['redshift_2']
title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
saveas = 'redshift_redshift'
plot_redshift_redshift(x, y, Nhex=200, xlabel=xlabel, ylabel=ylabel, diagonal=True, title=title, outpath=outpath, saveas=saveas)

## binned offset and STD from true redshift
xlabel = cfg['latex']['redshift_1']
title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
saveas = 'redshift_metrics_z1'
plot_redshift_std_and_mean(x, y, Nbins=100, x_true=True, title=title, outpath=outpath, saveas=saveas)


