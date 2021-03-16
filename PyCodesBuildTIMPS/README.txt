Python codes to build tracked IMERG mesoscale precipitation systems (TIMPS)

run_TIMPS.bash - runs the below python scripts in sequence

1) create_FiT_input_files.py - Creates input files for FiT algorithm using a set of thresholds.
2) write_tracking_output.py - Builds tracking output files.
3) process_FiTobs_from_originalfiles.py - Takes original IMERG files, FiT algorithm input files, and FiT algorithm output files and combines them into individual files for each precipitation system.
4) process_FiTobs_from_trackfiles.py - Takes tracked output files and combines them into individual files for each precipitation system.
5) add_vars.py - Calculates variables to quantify precipitation systems and writes them to the precipitation system files.
6) add_ERA5_vars.py - Adds ERA5 data to the precipitation system files.
