#!/bin/bash
#####################################################
# machine set up (users should change this part)
#####################################################
#
#SBATCH --time=6:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=zipser 
#SBATCH --partition=kingspeak
#SBATCH -J FiT
#SBATCH -o slurm-%j.out-%N 
#SBATCH -e slurm-%j.err-%N 
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=james.russell@utah.edu

stdbuf -i0 -o0 -e0 ./FiT_Object_analysis_basic_with_NetCDF4.exex namelist.dat | tee FiT.log

exit 0
