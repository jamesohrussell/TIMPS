#!/bin/bash
# Example script to run tracking, process objects into precipitation 
#  systems, and add various variables quantifying the precipitation 
#  systems.
#
# Important - must edit directories below, and namelist_TIPS.py and 
#  namelist.dat. These control the python scripts and FiT algorithm 
#  respectively.
#
# Run this script by typing: 
#  time ./run_TIPS.sh &
#  disown "process number e.g. 296651"
# Last line will disconnect process from the terminal window so that 
#  if the terminal closes the process continues running.
#
# Scripts must be run sequentially (each relies on the previous and 
#  you cannot open the same file with two scripts at the same time). 
#  wait command will ensure previous script ran before running next 
#  script.
#
# Below is an example/sample for running the scripts required.

# Inputs
datadir="/uufs/chpc.utah.edu/common/home/u0816744/varb2/FiT_CPEX-AW/"
pycodedir="/uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/"
trcodedir="/uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/exe_for_tracking_2020/"

# Make required directories
mkdir "${datadir}/FiT_input_test"  # FiT input data directory
mkdir "${datadir}/FiT_output_test" # FiT output data directory
mkdir "${datadir}/TIPS_test"       # TIPS files data directory
wait

# Python: Creates the thresholded files that are input to the FiT 
#  algorithm. Ensure all output is piped to a log file and script 
#  is run in the background or this will get messy in your terminal.
python -u "${pycodedir}create_FiT_input_files.py" > create.log 2>&1 &
wait

# Executable: Run compiled FiT algorithm (C++) to track objects
"${trcodedir}FiT_Object_analysis_basic_with_NetCDF4.exex" \
"${trcodedir}namelist.dat" > FiT.log 2>&1 &
wait

# Python: Process FiT objects into individual files, each describing 
#  one precipitation system.
python "${pycodedir}process_FiTobs.py" > proc.log 2>&1 &
wait

# Add variables quantifying precipitation systems
python "${pycodedir}add_vars.py" > add.log 2>&1 &
wait

# Add ERA5 variables to precipitation system files
python "${pycodedir}add_ERA5_vars.py" > addE5.log 2>&1 &
wait

# Remove extras
rm -f "${pycodedir}__pycache__/"
wait
