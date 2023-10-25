#!/bin/bash
echo "DESC setup"
#source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh 
#conda activate /sps/lsst/users/tguillem/DESC/desc_may_2021
cd /sps/lsst/users/tguillem/DESC/desc_april_2022
source setup.sh 
#conda env list

#import sys
#import python
echo "GO1"

#cd /sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/prepare_catalogs/
#cd /sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/matching/
#cd /sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/performance/
#cd /sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/pre_processing/
cd /sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/validation/
echo "GO2"

#python galaxy_selection.py ${healpix_pixel}
#python Skysim_selection.py ${healpix_pixel}
#python read_cosmoDC2.py
#python read_cosmoDC2_matching_skysim.py
#python clevar_matching.py
#python selection_function_halos.py
#python selection_function_halos_global_fit.py
#python cosmoDC2_skysim5000_matching.py
python footprint.py
echo "GO3"

#queues
#qsub -q mc_highmem_long Skysim_selection_batch.py
