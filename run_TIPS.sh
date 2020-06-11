#!/bin/bash
# Example script to run tracking, process objects into precipitation 
# systems, and add various variables quantifying the precipitation systems

# Create FiT input files
python -u /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/create_FiT_input_files.py > create.log 2>&1 &

# Wait for previous script to finish
wait

# Run FiT algorithm
/uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/exe_for_tracking_2020/FiT_Object_analysis_basic_with_NetCDF4.exex /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/exe_for_tracking_2020/namelist.dat > FiT.log 2>&1 &
wait

# Process FiT objects into individual files for precipitation systems
python /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/process_FiTobs.py > proc.log 2>&1 &
wait

# Add variables quantifying precipitation systems
python /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/add_vars.py > add.log 2>&1 &
wait

# Add ERA5 variables to precipitation system files
python /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/add_ERA5_vars.py > addE5.log 2>&1 &
