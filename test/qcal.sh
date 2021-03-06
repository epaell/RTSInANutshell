#!/bin/bash -l
#SBATCH --nodes=25
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --partition=gpuq
#SBATCH --account=mwasci
#SBATCH --export=NONE
#SBATCH --mem=30000
#SBATCH --gres=gpu:1
rm -f cal.log
./cal.py
