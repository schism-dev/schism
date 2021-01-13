#!/usr/bin/env bash

set -xeuo pipefail

export FC=mpif90
export F77=mpif77
export F90=mpif90
export CC=mpicc
export CC=mpicc

#Fix an error?
sed -i -e "s|\!bflux0,|bflux0,|g" src/Hydro/schism_step.F90

cd mk
# change key values
cp Make.defs.ubuntu.gnu.mpich Make.defs.local &&\
sed -i -e "s|ENV = ubuntu.gnu.mpich|ENV = conda|g" Make.defs.local && \
sed -i -e "s|FCS = gfortran|FCS = mpif90|g" Make.defs.local &&\
sed -i -e "s|FCP = /usr/bin/mpif90.mpich -fc=gfortran -ffree-line-length-none|FCP = mpif90 -ffree-line-length-none|g" Make.defs.local &&\
sed -i -e "s|CDFLIBS = -L/usr/lib|CDFLIBS = -L\$(PREFIX)/lib|g" Make.defs.local &&\
sed -i -e "s|CDFMOD = -I/usr/include/|CDFMOD = -I/\$(PREFIX)/include -I/usr/local/include -I/usr/local/include/malloc|g" Make.defs.local &&\
#sed -i -e "s|MTSLIBS = -L./ParMetis-3.1-Sep2010 -lparmetis -lmetis|MTSLIBS = -L\$(PREFIX)/lib -lparmetis -lmetis|g" Make.defs.local &&\
sed -i -e "s|include ../mk/include_modules|export TVD_LIM=VL|g" Make.defs.local && \
sed -i -e "s|#EXEC   := othername.ex|EXEC   := schism|g" Make.defs.local


# compile
cd ../src
make clean
make psc

cp schism $PREFIX/bin/
#cp libschism.a $PREFIX/lib/
