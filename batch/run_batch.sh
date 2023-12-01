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
#cd /sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/validation/
#cd /pbs/home/n/namourou/workspace/cluster_challenge/pre_processing/AMICO/
#cd /pbs/home/n/namourou/workspace/cluster_challenge/prepare_catalogs/
#cd /pbs/home/n/namourou/workspace/cluster_challenge/validation/
#cd /pbs/home/n/namourou/workspace/cluster_challenge/matching/
#cd /pbs/home/n/namourou/workspace/cluster_challenge/performance/
cd /pbs/home/n/namourou/workspace/side_codes/clusters/DC2/DC2_mask/production/
echo "GO2"

#python galaxy_selection.py ${healpix_pixel}
#python Skysim_selection.py ${healpix_pixel}
#python read_cosmoDC2.py
#python read_cosmoDC2_matching_skysim.py
#python clevar_matching.py cosmoDC2_v0 amico_cosmoDC2.fzb_magy proximity
#python selection_function_halos.py
#python selection_function_halos_global_fit.py
#python cosmoDC2_skysim5000_matching.py
#python footprint.py
#python amico_members_preproc.py ${healpix_pixel}
#python merge_amico_mb.py
#python read_amico.py
#python amico_verifications.py
#python completeness_m.py amico_DC2/mag_y/ p_matching
#python completeness_z.py amico_DC2/mag_y/ p_matching
#python completeness_2D.py amico_DC2/mag_y/ p_matching
#python purity.py amico_DC2/mag_y/ p_matching
#python purity_versus_completeness.py cosmoDC2_DC2.masked amico_DC2.fzb.magy_v0 proximity
python mask_dc2_holes.py ${healpix_pixel}

echo "GO3"

#queues
#qsub -q mc_highmem_long Skysim_selection_batch.py
