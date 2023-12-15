#!/bin/sh
#Different possibilities are amico_cosmoDC2, amico_redmapper, redmapper_cosmoDC2, amico_amico
#param2 p_matching mb_matching
param1=amico_cosmoDC2
param2=p_matching
param3=_sn6
sbatch -t 4:00:00 --cpus-per-task=1 --mem=100G --export=param1=${param1},param2=${param2},param3=${param3} run_py.sh
echo 'done'
    