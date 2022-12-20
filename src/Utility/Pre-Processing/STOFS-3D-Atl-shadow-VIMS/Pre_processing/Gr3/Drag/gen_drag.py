from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.fgrid import DragCoefficient

hgrid=Hgrid.open('./hgrid.gr3', crs='epsg:4326')
depth1=-1.0
depth2=-3.0
bfric_river=0.0025
bfric_land=0.025
fgrid=DragCoefficient.linear_with_depth(hgrid, depth1, depth2, bfric_river, bfric_land)
regions=['GoME+0.001.reg', 'Lake_Charles_0.reg']
values=[0.001, 0.0]
flags=[1, 0]
for reg, value, flag in zip(regions, values, flags):
    #fgrid.modify_by_region(hgrid, f'./{reg}', 0.001,depth1, 1)
    fgrid.modify_by_region(hgrid, f'./{reg}', value,depth1, flag)

fgrid.write('drag.gr3', overwrite=True)
