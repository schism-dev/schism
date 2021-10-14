#!/usr/bin/env bash

set -xeuo pipefail

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77
export F90=mpif90

# build and install schism
mkdir build
cd build

# check if needed on Linux
#CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/cross-linux.cmake")
#   ${CMAKE_PLATFORM_FLAGS[@]} \

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_Fortran_FLAGS_RELEASE_INIT="-O2 -ffree-line-length-none" \
    -DCMAKE_Fortran_FLAGS_DEBUG_INIT="-g -ffree-line-length-none" \
    -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT="-O2 -g -ffree-line-length-none" \
    -DC_PREPROCESS_FLAG="-cpp" \
    -DNetCDF_FORTRAN_DIR="$PREFIX/lib" \
    -DNetCDF_C_DIR="$PREFIX/lib" \
    -DNetCDF_LIBRARIES="$PREFIX/lib/libnetcdff.dylib;$PREFIX/lib/libnetcdf.dylib" \
    -DNetCDF_INCLUDE_DIR="$PREFIX/include" \
    -DTVD_LIM="VL" \
    ../src
#make -j${CPU_COUNT:-1} # It doesn't seem to work because of the build of parmetis
make

#make install
cp -r bin/* $PREFIX/bin/
#ln -s $PREFIX/bin/pschism_TVD-VL $PREFIX/bin/schism # optional to be discussed

#clean up
cd ..
rm -r build

cd $BUILD_PREFIX

