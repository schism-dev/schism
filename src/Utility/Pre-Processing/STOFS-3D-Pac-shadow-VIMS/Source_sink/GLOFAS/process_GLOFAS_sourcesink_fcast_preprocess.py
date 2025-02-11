import shapely as shl
import pylib as pl
import numpy as np
import rasterio as ras
import rasterio.features as rasfeat
import rasterio.plot as rasplot
import netCDF4 as nc
import scipy.ndimage as ndimage
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import os
import datetime
import sys
######################################
# model time range from param.nml
rundir = sys.argv[1] # location of hgrid.npz or hgrid.gr3 file
glofasdir = rundir # location of GLOFAS upriver area file
outdir = rundir # location to write preprocessed GLOFAS information to, this is the file used by the daily script
minarea = 0 # minimum cutoff for rivers to include based on watershed area in square km 
outsuffix = f'forecast-preproc'  #  appended to image filenames

# read hgrid
try:
  hgrid = pl.read(os.path.join(rundir,'hgrid.npz'))
except:
  hgrid = pl.read(os.path.join(rundir,'hgrid.gr3'))
hgrid.compute_bnd()

# generate shapely polygon of grid boundary
I = hgrid.bndinfo
prj='epsg:4326'
x = I.x
nbn = I.nbn
y = I.y
ring = np.array([x[0:nbn[0]], y[0:nbn[0]]]).T
holes = []
for i in range(1,I.nb):
  st = np.sum(nbn[:i])
  en = np.sum(nbn[:i+1])
  hole = np.array([x[st:en],y[st:en]]).T
  holes.append(hole)
poly = shl.Polygon(ring,holes=holes)


plotbounds = {'Columbia':[235.5,238,44.5,47.5],'VancIs':[231,237, 48.25, 52.1]}
# read in GLOFAS upriver basin area
uparea = ras.open(os.path.join(glofasdir,'uparea_glofas_v4_0.nc'))
band0 = uparea.read(1,masked=True)
# convert from -180 -180 to 0-360 longitude rage
band1 = np.c_[band0[:,3600:],band0[:,:3600]] # swap order of hemispheres
tr = np.array(uparea.transform)
if tr[2] == -180:
  tr[2] = 0  # shift starting point of transform from 
transf360 = ras.Affine(*tr) # create new Affine transform
processplot = True  # set to False to deactivate images
plotarea = 'VancIs'
if processplot:
  f,ax = plt.subplots(1,1)
  rasplot.show(band1,transform=transf360,vmin=0,vmax=10000e6,ax=ax)
  ax.axis(plotbounds[plotarea])
  ax.set_title(f'GLOFAS discharge around {plotarea}')
  f.savefig(f'glofas_{plotarea}_{outsuffix}_step0.png')
  plt.close(f)
# mask area covered by SCHISM grid

rasfeat.rasterize([(poly,-3.4028235e+38)],transform=transf360, all_touched=True, out=band1)
if processplot:
  f,ax = plt.subplots(1,1)
  rasplot.show(band1,transform=transf360,vmin=0,vmax=10000e6,ax=ax)
  ax.axis(plotbounds[plotarea])
  ax.set_title('burn in SCHISM mesh')
  f.savefig(f'glofas_{plotarea}_{outsuffix}_step0.png')
  plt.close(f)

# find river outlets by finding cells with local maximum uparea
nbh = 3
upmax = ndimage.maximum_filter(band1,nbh) # 3x3 windowed max onto uparea
maxima = upmax==band1 # find local maxima in uparea
riverendA = maxima * band1/1e6 #  convert from m^2 to km^2
ind0,ind1 = np.where(riverendA > minarea) # limit river ends to watersheds > [minarea] km ^2
lon, lat = ras.transform.xy(transf360,ind0,ind1) # get lat-lon of riverends
riverendarea = riverendA[ind0,ind1]
riverends = shl.points(np.array(list(zip(lon,lat)))) # contruct shapely set of points of river ends

outerbuf = shl.buffer(poly,1) # extend the schism grid polygon by 1 degree
D = outerbuf.distance(riverends) # find distance from riverends to expanded polygon
nearrivers =D ==0  # limit riverends to within 1 degree from schism grid

lonN = lon[nearrivers]
latN = lat[nearrivers]
riverCoords = np.c_[lonN,latN]
if processplot:
  f,ax = plt.subplots(1,1)
  rasplot.show(riverendA,transform=transf360,vmin=0,vmax=10000,ax=ax)
  plt.plot(lonN,latN,'m*')
  ax.axis(plotbounds[plotarea])
  ax.set_title(f'river ends {plotarea}')
  f.savefig(f'glofas_{plotarea}_{outsuffix}_step2.png')
  plt.close(f)
# find boundary element side-center nearest each river, maintaining connection to element id
hgrid.compute_ctr() # compute element centers, used for finding elements in region
hgrid.compute_side(2)   # compute table of side elements and side center coordinates
# 2) Generate element centroid KDTree
boundedges = np.where(hgrid.isdel[:,1] == -1)[0]  # isdel is the table of elements for each side, sides that only have 1 element have -1 as the second element, and are boundary sides
sidecenters = np.c_[hgrid.xcj[boundedges], hgrid.ycj[boundedges]]   # xcj,ycj is the side center point
# create a KDTree for boundary side centers
tree = cKDTree(sidecenters) 

# search for nearest boundary side center for all river
_, idxs = tree.query(np.vstack(riverCoords), workers=-1)
# idxs is the index of the side nearest to the river end points in the data subset
# plug idxs back into boundedges and isdel to get the element that contains the sides
Erivers = hgrid.isdel[boundedges[idxs],0]
######################################


# indices for final rivers into original full domain, including openbnd
ind0N = ind0[nearrivers]
ind1N = ind1[nearrivers]
lonN = lon[nearrivers]
latN = lat[nearrivers]
areaN = riverendarea[nearrivers]

# raster of final locations, including open boundaries, just used for process visualization
final =band1*0
final[ind0N,ind1N] = band1[ind0N,ind1N]
#
#
ncid = nc.Dataset(os.path.join(glofasdir,'uparea_glofas_v4_0.nc'))
lon = ncid.variables['longitude'][:]
g = lon<=0
lon[g] += 360
print(lon[~g])
lon = np.r_[lon[~g],lon[g]] # rearrange longitude to match data files
lat = ncid.variables['latitude'][:]

# save preprocessed results, including selected hgrid information
np.savez(os.path.join(glofasdir,'preproc-glofas.npz'),hgrid=hgrid,Erivers=Erivers,idxs=idxs,ind0N=ind0N,ind1N=ind1N, lonN=lonN,latN=latN,lon=lon,lat=lat,areaN=areaN,final=final,
  xctr = hgrid.xctr,
  yctr= hgrid.yctr,
  elnode = hgrid.elnode,
  iobn =  hgrid.iobn ,
  ne = hgrid.ne, x = hgrid.x, y = hgrid.y)
print(f'saved preproc-glofas.npz')
