#!/bin/sh

mpicxx -O2 -Wall -Wno-unused-result -Wno-unknown-pragmas -o FiT_Object_analysis_basic_with_NetCDF4_and_MPI.exex CC_Object_analysis_FiT_basic_with_NetCDF4_and_MPI.cc -I/usr/include/hdf5/serial/ -lboost_serialization -lgd -lnetcdf 

#mpicxx -O2 -Wall -Wno-unused-result -Wno-unknown-pragmas -o CC_Object_analysis_FiT_basic_with_NetCDF4_and_MPI.exex CC_Object_analysis_FiT_basic_with_NetCDF4_and_MPI.cc -I${PATH_TO_PNGWRITER_HEADER_FILE} -L${PATH_TO_PNGWRITER_LIBRARY_FILE} -DNO_FREETYPE -lboost_serialization -lmpi -lgd -ljpeg -lnetcdf_c++ -lnetcdf -lpngwriter -lpng -lhdf5_hl -lhdf5 -lfreetype -lcurl -lz

#mpicxx -O2 -Wall -Wno-unused-result -o CC_Prepare_compressed_netcdf_files_MPI.exex CC_Prepare_compressed_netcdf_files_MPI.cc -I${PATH_TO_PNGWRITER_HEADER_FILE} -L${PATH_TO_PNGWRITER_LIBRARY_FILE} -DNO_FREETYPE -lboost_serialization -lmpi -lgd -ljpeg -lnetcdf_c++ -lnetcdf -lpngwriter -lpng -lhdf5_hl -lhdf5 -lfreetype -lcurl -lz

