#!/usr/bin/env bash

set -xeuo pipefail

#export mpi=mpich
export CC=mpicc
export FC=mpif90
export F77=mpif77
export F90=mpif90

#Fix an error?
sed -i -e "s|\!bflux0,|bflux0,|g" src/Hydro/schism_step.F90

# build and install schism
mkdir build
cd build

CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/cross-linux.cmake")

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_Fortran_FLAGS_RELEASE_INIT="-O2 -ffree-line-length-none" \
    -DCMAKE_Fortran_FLAGS_DEBUG_INIT="-g -ffree-line-length-none" \
    -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT="-O2 -g -ffree-line-length-none" \
    -DC_PREPROCESS_FLAG="-cpp" \
    -DNetCDF_FORTRAN_DIR="$PREFIX/lib" \
    -DNetCDF_C_DIR="$PREFIX/lib" \
    -DNetCDF_INCLUDE_DIR="$PREFIX/include" \
    -DNetCDF_LIBRARIES="$PREFIX/lib/libnetcdff.so;$PREFIX/lib/libnetcdf.so" \
    ${CMAKE_PLATFORM_FLAGS[@]} \
    -DTVD_LIM="VL" \
    ../src
#make -j${CPU_COUNT:-1}
#   -DNetCDF_LIBRARIES="$PREFIX/lib/libnetcdff.dylib;$PREFIX/lib/libnetcdf.dylib" \
make
cp -r bin/* $PREFIX/bin/
cp -r include/* $PREFIX/include/
cp -r lib/* $PREFIX/lib/
ln -s $PREFIX/bin/pschism_TVD-VL $PREFIX/bin/schism
#cd $BUILD_PREFIX
