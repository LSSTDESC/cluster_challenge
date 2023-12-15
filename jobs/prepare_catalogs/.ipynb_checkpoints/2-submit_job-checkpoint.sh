#!/bin/sh

sbatch -t 1:00:00 --cpus-per-task=1 --mem=64G run_py.sh
echo 'done'
