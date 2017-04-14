#!/bin/bash -l
#SBATCH --nodes=25
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --partition=workq
#SBATCH --account=mwasci
#SBATCH --export=NONE
rm -f image.log
./image.py
