

import yaml, sys
import numpy as np
from astropy.table import Table

from clevar.catalog import ClCatalog
from clevar.match import get_matched_pairs, output_matched_catalog
from clevar.match_metrics import scaling

import matplotlib.pyplot as plt
from plot_style import set_style, figsize
import saving_figures as sfigs
import fit_func as ff
from plotting_functions import plot_redshift_redshift, plot_zscatter_redshift, plot_redshift_std_and_mean


## Initialize style for plots.
set_style()

# read config files as online arguments 
config = sys.argv[1]

# open config files
with open(config) as fstream:
	cfg = yaml.safe_load(fstream)

cat1 = cfg['cats']['cat1']
cat2 = cfg['cats']['cat2']
mt_method = cfg['matching']['method']
mt_pref = cfg['matching']['pref']
if mt_method == 'member' :
        mt_params = [cfg['matching']['minimum_share_fraction']]
else :
        mt_params = [cfg['matching']['delta_z'], cfg['matching']['match_radius']]

inpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['in'])
outpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['out'],
        addon='redshift')



## Create a merged catalog for the cross-matched pairs if it does not already exist.
cltags = cfg['merged_keys']
try :
	mtcats = ClCatalog.read(f"{inpath}{cat1}_{cat2}_merged.fits", 'merged', tags=cltags)
except :
	cats = {cat1: ClCatalog.read_full(f"{inpath}{cat1}.fits"),
		cat2: ClCatalog.read_full(f"{inpath}{cat2}.fits")}
	output_matched_catalog(f"{inpath}{cat1}.fits", f"{inpath}{cat2}.fits", f"{inpath}{cat1}_{cat2}_merged.fits",
		cats[cat1], cats[cat2], matching_type='cross', overwrite=True)
	mtcats = ClCatalog.read(f"{inpath}{cat1}_{cat2}_merged.fits", 'merged', tags=cltags)


## redshift vs redshift
x = mtcats['z_halo']
y = mtcats['z_cl']
xlabel = '$z_{halo}$'
ylabel = '$z_{cl}$'
title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
saveas = 'redshift_redshift'
plot_redshift_redshift(x, y, Nhex=200, xlabel=xlabel, ylabel=ylabel, diagonal=True, title=title, outpath=outpath, saveas=saveas)

## redshift scatter vs redshift
saveas = 'zscatter_redshift'
plot_zscatter_redshift(x, y, Nhex=200, cfg=cfg, title=title, outpath=outpath, saveas=saveas)

## binned offset and STD from true redshift
xlabel = "$z_{halo}$"
title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
saveas = 'redshift_metrics_zhalo'
plot_redshift_std_and_mean(x, y, Nbins=100, x_true=True, title=title, outpath=outpath, saveas=saveas)
saveas = 'redshift_metrics_zcl'
plot_redshift_std_and_mean(x, y, Nbins=100, x_true=False, title=title, outpath=outpath, saveas=saveas)


