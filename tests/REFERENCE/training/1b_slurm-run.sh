#!/bin/bash
#SBATCH --job-name=fiber_array
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=1
#SBATCH --partition=HM-dev
#SBATCH --array=00-40

##Optional - Required memory in MB per core. Defaults are 1GB per core.
##SBATCH --mem-per-cpu=3072

##Optional - Estimated execution time
##Acceptable time formats include  "minutes",   "minutes:seconds",
##"hours:minutes:seconds",   "days-hours",   "days-hours:minutes" ,"days-hours:minutes:seconds".
##SBATCH --time=

########### Further details -> man sbatch ##########

export OMP_NUM_THREADS=1
printf -v ID "%02d\n" $SLURM_ARRAY_TASK_ID
TRAJECTORYPATH=/comp-des-mat/mraschi/Cases/REFERENCE_2branches_damage/training/trajectory_$ID
cd $TRAJECTORYPATH
time python3 MainKratos.py

