# Python codes to build tracked IMERG precipitation systems
#
# 1) create_FiT_input_files - Creates input files for FiT algorithm using 
#     a set of thresholds
# 2) process_FiTobs - Takes original IMERG files, FiT algorithm input files,
      and FiT algorithm output files and combines them into individual files
      for each precipitation system
# 3) add_vars - Calculates variables to quantify precipitation systems and 
      writes them to the precipitation system files
# 4) add_ERA5_vars - Adds ERA5 data to the precipitation system files
