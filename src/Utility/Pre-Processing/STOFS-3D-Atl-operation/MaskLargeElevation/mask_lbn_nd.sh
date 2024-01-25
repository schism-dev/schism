#!/bin/bash

infile=out2d_1.nc
outfile=out2d_1.mask.nc

ncks -A -v idmask lbnd_node_mask.nc $infile
ncap2 -s 'where(idmask==1) elevation=float(-99999.)' $infile $outfile
ncatted -O -a missing_value,elevation,a,f,-99999. $outfile
