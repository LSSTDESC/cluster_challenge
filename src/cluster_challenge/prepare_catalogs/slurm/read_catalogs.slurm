#!/bin/bash
#SBATCH --job-name=reduce_cluster_catalogs
##SBATCH --output=slurm_output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=100:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=fail
##SBATCH --mail-type=end
##SBATCH --mail-user=solomon@lapp.in2p3.fr


source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/

if [ -n $SLURM_JOB_ID ]; then
        SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
        SCRIPT_PATH=$(realpath $0)
fi

SHARED_PATH=$(dirname  $SCRIPT_PATH)

exec > $SHARED_PATH/slurm_output.out

echo $(date -u +%Y%m%d-%H:%M:%S)


python -u $SHARED_PATH/../read_cosmoDC2.py v0		## SMALL
python -u $SHARED_PATH/../read_cosmoDC2.py v1		## FULL

python -u $SHARED_PATH/../read_wazp.cosmoDC2.py v0		## RICHNESS 0
python -u $SHARED_PATH/../read_wazp.cosmoDC2.py v1		## RICHNESS 20	

python -u $SHARED_PATH/../read_wazp.cosmoDC2.fzb.py v0	## RICHNESS 0
python -u $SHARED_PATH/../read_wazp.cosmoDC2.fzb.py v1	## RICHNESS 20	

python -u $SHARED_PATH/../read_wazp.DC2.py v0		## RICHNESS 0
python -u $SHARED_PATH/../read_wazp.DC2.py v1		## RICHNESS 20

python -u $SHARED_PATH/../read_redmapper.cosmoDC2.py v0	## RICHNESS 0
python -u $SHARED_PATH/../read_redmapper.cosmoDC2.py v1	## RICHNESS 20

python -u $SHARED_PATH/../read_redmapper.DC2.py v0		## RICHNESS 0
python -u $SHARED_PATH/../read_redmapper.DC2.py v1		## RICHNESS 20

echo $(date -u +%Y%m%d-%H:%M:%S)
