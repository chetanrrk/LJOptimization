#!/bin/bash
#SBATCH --job-name=opt_run
#SBATCH --account=Drude
#SBATCH --partition=sball
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=96:00:00
#SBATCH -o opt.out
#SBATCH -e opt.err

../non_halo_code/optimizeLennardJones lj-sim-blues-working.py

