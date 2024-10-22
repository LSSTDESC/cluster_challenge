#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=matching
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge//../..//slurm_outputs/runs/EXAMPLE_RUN/matching.out
#SBATCH --mem=8GB
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/matching/matching.py runConfigs/example.cfg
