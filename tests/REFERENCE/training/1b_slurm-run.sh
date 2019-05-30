#!/bin/bash
#SBATCH --job-name=fiber_array
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=1
#SBATCH --array=00-40

## Settings for "fiber2" case:
##
#SBATCH --partition=HM-dev
#SBATCH --mem-per-cpu=1024
#SBATCH --time=01:00:00

## Settings for "fiber3" case:
##
##SBATCH --partition=HM
##SBATCH --mem-per-cpu=????
##SBATCH --time=????


########### Further details -> man sbatch ##########

export OMP_NUM_THREADS=1
printf -v ID "%02d\n" $SLURM_ARRAY_TASK_ID
TRAJECTORYPATH=$PWD/trajectory_$ID
cd $TRAJECTORYPATH
time python3 MainKratos.py

