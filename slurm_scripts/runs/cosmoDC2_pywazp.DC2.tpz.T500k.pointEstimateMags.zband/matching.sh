#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=matching
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/runs/cosmoDC2_pywazp.DC2.tpz.T500k.pointEstimateMags.zband/matching.out
#SBATCH --mem=64GB
module load conda
conda activate /sps/lsst/users/rsolomon/conda_envs/clevar/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/matching/matching.py runConfigs/cosmoDC2_pywazp.DC2.tpz.T500k/zband/member_matching/fshare_0.0_pref_more_massive.cfg
