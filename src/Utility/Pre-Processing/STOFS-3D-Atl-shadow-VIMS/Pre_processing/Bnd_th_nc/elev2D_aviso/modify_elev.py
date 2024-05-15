#!/usr/bin/env python3
#add a linear shift on elev2D based on lat
from pylib import *

#input
lats=[0, 40, 40.5, 90]
msl_shifts=[-0.42, -0.42, -0.42, -0.42]

#read grid and elev
if os.path.exists('hgrid.npz'):
   gd=loadz('hgrid.npz').hgrid
else:
   gd=read_schism_hgrid('hgrid.gr3')
C=ReadNC('elev2D.th.nc.aviso')

#prepare for the shift for the first boundary
bnds = [0, 1]
ind1 = 0
ind2 = 0
for i in bnds:
    slat = gd.y[gd.iobn[i]]
    elev_adjust=np.interp(slat,lats,msl_shifts)
    ind1 = ind2
    ind2 = ind1 + gd.iobn[i].shape[0]

    #modify elev
    C.time_series.val[:, ind1:ind2, :, :] = C.time_series.val[:, ind1:ind2, :, :] +elev_adjust[None,:,None,None]

WriteNC('elev2D.th.modify.nc',C)

