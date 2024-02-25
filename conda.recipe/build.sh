#!/usr/bin/env bash

set -xeuo pipefail

#export mpi=mpich
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77
export F90=mpif90

#Remove metis from cmake
sed -i -e "s|set(PARMETIS_VER|#set(PARMETIS_VER|g" src/CMakeLists.txt
sed -i -e "s|set(PARMETIS_DIR|#set(PARMETIS_DIR|g" src/CMakeLists.txt
sed -i -e "s|add_subdirectory( \${PARMETIS_VER} )|#add_subdirectory( \${PARMETIS_VER} )|g" src/CMakeLists.txt

# build and install schism
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_Fortran_FLAGS_RELEASE_INIT="-O2 -ffree-line-length-none" \
    -DTVD_LIM=$TVD_LIM \
    -DOLDIO=$OLDIO \
    -DPREC_EVAP=$PREC_EVAP \
    -DUSE_GOTM=$GOTM \
    -DUSE_HA=$HA \
    -DUSE_SED2D=$SED2D \
    -DUSE_MARSH=$MARSH \
    -DUSE_WWM=$WWM \
    -DUSE_WW3=$WW3 \
    -DUSE_ICE=$ICE \
    -DUSE_ICM=$ICM \
    -DUSE_GEN=$GEN \
    -DUSE_AGE=$AGE \
    -DUSE_ECO=$ECO \
    -DICM_PH=$PH \
    -DUSE_COSINE=$COSINE \
    -DUSE_FIB=$FIB \
    -DUSE_FABM=$FABM \
    -DUSE_SED=$SED \
    ../src

make

##make install
cp -r bin/* $PREFIX/bin/
#make a symlink for convenience
ln -s $PREFIX/bin/pschism_TVD-$TVD_LIM $PREFIX/bin/schism

#clean up
cd ..
rm -r build

cd $BUILD_PREFIX
