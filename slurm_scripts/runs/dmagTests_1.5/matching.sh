#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=matching
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/runs/dmagTests_1.5/matching.out
#SBATCH --mem=8GB
module load conda
conda activate /sps/lsst/users/rsolomon/conda_envs/clevar/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/matching/matching.py runConfigs/dmag_tests/dmag_1.5.cfg
