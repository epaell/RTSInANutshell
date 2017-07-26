#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:15:00
#SBATCH --partition=gpuq
#SBATCH --account=mwasci
#SBATCH --export=NONE
#
#
./launcher.py
