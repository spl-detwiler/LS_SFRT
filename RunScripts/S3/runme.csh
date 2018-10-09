#!/bin/bash
#$ -N v5v3rstrt
#$ -q cee
#$ -m beas
#$ -ckpt blcr
#$ -pe openmp 32
#$ -cwd

module load MATLAB/r2017b

mcc -m ADRE.m -a ./../flowandtransportscripts -a ./../LevelSetAdvanceFunctions

uptime
./ADRE
uptime

module unload MATLAB/r2017b
