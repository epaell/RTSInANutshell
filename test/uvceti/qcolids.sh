#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --partition=workq
#SBATCH --account=mwasci
#SBATCH --export=NONE
#SBATCH --array=0
# Flag bad channels before averaging over the band
aprun -n 1 sflag.py ${SLURM_ARRAY_TASK_ID} 1.8
# Average remaining channels over the band
aprun -n 1 collapseobsids.py ${SLURM_ARRAY_TASK_ID}
