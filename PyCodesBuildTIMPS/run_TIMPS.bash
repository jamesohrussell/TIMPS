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
datadir="/uufs/chpc.utah.edu/common/home/u0816744/zips2/TIMPS/"
infiledir="FiTin_2015/"
outfiledir="FiTout_2015/"
Trkfiledir="Tracking_2015/"
TIMPSfiledir="TIMPS_2015/"
pycodedir="/uufs/chpc.utah.edu/common/home/zipser-group2/TIMPS/scripts_2015/"
trcodedir="/uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/exe_for_tracking_2020_v3/"

# Make required directories (should line up with those in python namelist)
echo "Making directories for output"
mkdir $datadir             # TIPS data directory
mkdir $datadir$infiledir   # FiT input data directory
mkdir $datadir$outfiledir  # FiT output data directory
mkdir $datadir$Trkfiledir  # Tracking data directory
mkdir $datadir$TIMPSfiledir # TIPS files data directory
wait

# Python: Creates the thresholded files that are input to the FiT 
#  algorithm. Ensure all output is piped to a log file and script 
#  is run in the background or this will get messy in your terminal.
echo "Creating FiT input files"
python -u "${pycodedir}create_FiT_input_files.py" > create.log 2>&1 &
wait

# Executable: Run compiled FiT algorithm (C++) to track objects
echo "Running FiT tracking algorithm"
"${trcodedir}FiT_Object_analysis_basic_with_NetCDF4.exex" \
"${datadir}${infiledir}namelist_FiT.txt" > FiT.log 2>&1 &
wait

# Executable: Run compiled FiT algorithm (C++) to track objects
echo "Compiling FiT tracking data"
python -u "${pycodedir}write_tracking_output.py" \
> writetrack.log 2>&1 &
wait

# Python: Process FiT objects into individual files, each describing 
#  one precipitation system.
echo "Processing FiT objects into TIPS"
python -u "${pycodedir}process_FiTobs_from_trackfiles.py" \
> proc.log 2>&1 &
wait

# Add variables quantifying precipitation systems
echo "Adding variables to TIPS"
python -u "${pycodedir}add_vars.py" > add.log 2>&1 &
wait

# Add ERA5 variables to precipitation system files
#echo "Adding ERA5 variables to TIPS"
#python -u "${pycodedir}add_ERA5_vars.py" > addE5.log 2>&1 &
#wait

# Remove extras
rm -rf "${pycodedir}__pycache__/"
wait

# Finish message
echo "Tracking and processing complete"
