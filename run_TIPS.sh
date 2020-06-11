#!/bin/bash
# Example script to run tracking, process objects into precipitation 
# systems, and add various variables quantifying the precipitation systems
#
# Important - Namelists can be found at the start of each script:
#  * create_FiT_input_files.py
#  * namelist.dat (for FiT algorithm)
#  * process_FiTobs.py
#  * add_vars.py
#  * add_ERA5_vars.py
# Each namelist controls the inputs for each script. 
# Make sure to change all namelists before running this script.
#
# Run this script by typing: 
#  ./run_TIPS.sh &
#  disown "process number e.g. 296651"
# Last line will disconnect process from the terminal window so that if the 
# terminal closes the process continues running.
#
# Scripts must be run sequentially (each relies on the previous and you cannot
# open the same file with two scripts at the same time).
# wait command will ensure previous script ran before running next script.
#
# Below is an example/sample for running the scripts required.

# Python: Creates the thresholded files that are input to the FiT algorithm
# Ensure all output is piped to a log file and script is run in the background
# or this will get messy in your terminal.
python -u /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/create_FiT_input_files.py > create.log 2>&1 &
wait

# Executable: Run compiled FiT algorithm (C++) to track objects
/uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/exe_for_tracking_2020/FiT_Object_analysis_basic_with_NetCDF4.exex /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/exe_for_tracking_2020/namelist.dat > FiT.log 2>&1 &
wait

# Python: Process FiT objects into individual files, each describing one 
# precipitation system.
python /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/process_FiTobs.py > proc.log 2>&1 &
wait

# Add variables quantifying precipitation systems
python /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/add_vars.py > add.log 2>&1 &
wait

# Add ERA5 variables to precipitation system files
python /uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/add_ERA5_vars.py > addE5.log 2>&1 &
