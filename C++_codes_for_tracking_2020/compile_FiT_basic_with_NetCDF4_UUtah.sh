#!/bin/bash

set -e

export FILENAME=CC_Object_analysis_FiT_basic_with_NetCDF4
export EXE=FiT_Object_analysis_basic_with_NetCDF4.exex

export NETCDF_DIR=/uufs/chpc.utah.edu/sys/installdir/netcdf/4.1.3
export NETCDF_INCDIR=$NETCDF_DIR/include
export NETCDF_LIBDIR=$NETCDF_DIR/lib

export GD_DIR=/uufs/chpc.utah.edu/sys/installdir/libgd/2.2.5
export GD_INCDIR=$GD_DIR/include
export GD_LIBDIR=$GD_DIR/lib


echo "  STEP1: Compiling the source ... "
g++ -c -I$NETCDF_INCDIR -I$GD_INCDIR ${FILENAME}.cc

echo "  STEP2: Linking the object files into an executable "
g++ -o $EXE ${FILENAME}.o \
       -Wl,-rpath=$NETCDF_LIBDIR -L$NETCDF_LIBDIR -lnetcdf_c++ -lnetcdf \
       -Wl,-rpath=$GD_LIBDIR  -L$GD_LIBDIR -lgd   

echo "  STEP3: Cleanup "
rm -rf ${FILENAME}.o


echo "  DONE !"
