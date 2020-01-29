#!/bin/bash
#####################################################
# machine set up (users should change this part)
#####################################################
#
#SBATCH --time=12:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=zipser 
#SBATCH --partition=kingspeak
#SBATCH -J IMERGtoFiT
#SBATCH -o slurm-%j.out-%N 
#SBATCH -e slurm-%j.err-%N 
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=james.russell@utah.edu

stdbuf -i0 -o0 -e0 python create_FiT_input_files.py  | tee IMERGtoFiT.log

exit 0
