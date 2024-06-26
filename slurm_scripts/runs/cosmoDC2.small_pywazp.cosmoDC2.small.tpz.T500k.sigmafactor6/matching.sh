#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=matching
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/runs/cosmoDC2.small_pywazp.cosmoDC2.small.tpz.T500k.sigmafactor6/matching.out
#SBATCH --mem=8GB
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/matching/matching.py runConfigs/cosmoDC2.small_pywazp.cosmoDC2.small.tpz.T500k.sigmafactor6/member_matching/fshare_0.0_pref_more_massive.cfg
