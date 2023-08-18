#sbatch -p flash --mem=8G run_batch.sh
#sbatch --mem=64G run_batch.sh
#sbatch -p flash run_batch.sh
sbatch -t 12:00:00 --cpus-per-task=1 --mem=64G  run_batch.sh
