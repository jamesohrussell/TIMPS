#!/bin/bash
#####################################################
# machine set up (users should change this part)
#####################################################
#
#SBATCH --time=12:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --account=zipser 
#SBATCH --partition=kingspeak
#SBATCH -J addPFvars
#SBATCH -o slurm-%j.out-%N 
#SBATCH -e slurm-%j.err-%N 
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=james.russell@utah.edu

stdbuf -i0 -o0 -e0 python add_PF_vars_parrallel.py | tee addvars.log

exit 0
