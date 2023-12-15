#!/bin/sh
#Different possibilities are different path for dfferent cats ex : amico_cosmoDC2/mag_i/0.0/
#param2 p_matching mb_matching
param1=amico_cosmoDC2/mag_i/
param2=p_matching
sbatch -t 02:00:00 --cpus-per-task=1 --mem=20G --export=param1=${param1},param2=${param2} run_py.sh
echo 'done'