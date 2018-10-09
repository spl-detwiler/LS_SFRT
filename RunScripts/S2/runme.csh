#!/bin/bash
#$ -N Noiriel
#$ -q cee
#$ -m beas
#$ -pe openmp 24
#$ -cwd

module load MATLAB/r2017b

mcc -m ADRE.m -a ./../flowandtransportscripts -a ./../LevelSetAdvanceFunctions

cd /dfs3/pub/trevorj/repos/LSM_Newv2/Noiriel/

uptime
./ADRE
uptime

module unload MATLAB/r2017b
