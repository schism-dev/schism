#!/usr/bin/env python3
#add a linear shift on elev2D based on lat
from pylib import *

#input
lats=[0, 27, 28, 32, 33, 90]
msl_shifts=[-0.25, -0.25, 0.0, 0.0, 0.56, 0.56]

#read grid and elev
if os.path.exists('../hgrid.npz'):
   gd=loadz('../hgrid.npz').hgrid
else:
   gd=read_schism_hgrid('../hgrid.gr3')
C=ReadNC('elev2D.th.nc')

#prepare for the shift
slat=gd.y[gd.iobn[0]]
elev_adjust=np.interp(slat,lats,msl_shifts)

#modify elev
C.time_series.val=C.time_series.val+elev_adjust[None,:,None,None]
WriteNC('elev2D.th.modify.nc',C)

#change link
bdir=os.path.abspath(os.path.curdir)
os.system('cd ../; rm elev2D.th.nc; ln -sf {}/elev2D.th.modify.nc elev2D.th.nc'.format(bdir))

