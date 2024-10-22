#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=performance
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/runs/cosmoDC2.small_pywazp.cosmoDC2.small.sigz0.01.zband/performance.out
#SBATCH --mem=8GB
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/performance/performance.py runConfigs/cosmoDC2.small_pywazp.cosmoDC2.small.sigz0.01/zband/member_matching/fshare_0.0_pref_more_massive.cfg
