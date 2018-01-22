#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --partition=workq
#SBATCH --account=mwasci
#SBATCH --export=NONE
# Run RM synthesis on the cube generated for the obsid
cd 1165925976
# Work on RTS cube directly (i.e. images prefixed with "2"), search cube from phi=-60 to 60 in steps of 0.5 rad/m^2, don't clean or write out a cube.
rmsynth.py 2 -60.0 60.0 0.5 0 0
