import numpy as np
import json
import argparse
import geopandas
import shapely as shl
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pylib as pl
import os.path
import netCDF4 as nc

def getreaches(fname, R = [], E=[]):
   with open(fname) as f:
     print(fname)
     d = json.load(f)
     for k,v in d.items():
        for vv in v: R.append(vv)
        E.append(int(k)-1)
   return R,E,d 

def loadSourceSink(fn):
    d = np.loadtxt(fn,dtype = int)
    source = d[1:d[0]+1] - 1 # one based to zero based, first value is number of locations
    return source
'''
Matches NWM reaches to specified elements in source_sink.in
'''

regions = ['alaska', 'conus','hawaii']
nwmfs = {'conus':'NWMv3.0_conus_2017-11-01-2019-01-01_daily.nc', 
          'hawaii':'NWMv3.0_hawaii_2012-11-01-2014-11-01_daily.nc',
          'alaska':'NWMv3.0_alaska_2017-11-01-2019-01-01_daily.nc'}
distcut = {'hawaii':10000,'alaska':35000,'conus':10000}
basedir = '/sciclone/home/cseaton/data10/'
rundir = os.path.join(basedir,'pacific/GLOFAS/source_sink/')
if os.path.exists(os.path.join(rundir,'hgrid.npz')):
   hgrid = pl.read(os.path.join(rundir,'hgrid.npz'))
else:
   hgrid = pl.read(os.path.join(rundir,'hgrid.gr3'))
hgrid.compute_ctr()
hgrid.compute_bnd()
Gdir = os.path.join(basedir,'pacific/GLOFAS/source_sink/')
Gbp = pl.read(os.path.join(Gdir,'vsource.0.bp'))
for region in regions:
  if region == 'alaska' or region == 'hawaii':
    westcoast = np.logical_and(Gbp.x > 200,Gbp.x  < 230)
  if region == 'conus':
    westcoast = np.logical_and(Gbp.x > 220,Gbp.x  < 260)
  GE = loadSourceSink(os.path.join(Gdir,'source_sink.in.0'))
  ycG = hgrid.yctr[GE[westcoast]]
  xcG = hgrid.xctr[GE[westcoast]]
  GEwest = GE[westcoast]
  Gpts = np.c_[Gbp.x[westcoast],Gbp.y[westcoast]]
  G_shape = [shl.Point(s) for s in Gpts ]
  
  prjstr = {'conus':'EPSG:32610',
          'hawaii':'EPSG:32606',
   'alaska':'EPSG:32606'}[region]
  
  prj='epsg:4326'
  
  gdfGpts= geopandas.GeoDataFrame({ 'geometry':G_shape,'meanQ':Gbp.z[westcoast]}, crs="EPSG:4326")
  gdfGpts_paclcc = gdfGpts.to_crs(prjstr)
  Glcc = np.c_[gdfGpts_paclcc.geometry.x,gdfGpts_paclcc.geometry.y]
  
  Ndir = os.path.join(basedir,'NWM/')
  
  ncid = nc.Dataset(os.path.join(Ndir,nwmfs[region]))
  discharge = ncid.variables['streamflow'][:,:]
  meandis = np.mean(discharge,axis=0)
  maxdis = np.max(discharge,axis=0)
  maxdis[maxdis<1] = 1
  meandis[meandis<1] = 1
  
  dlon = ncid.variables['longitude'][:] + 360
  dlat = ncid.variables['latitude'][:]
  reaches = ncid.variables['feature_id'][:]
  #tree = KDTree(Nlcc)
  # search for nearest boundary side center for all river
  #_, idxs = tree.query(np.vstack(Glcc), workers=-1)
  #dist = np.sqrt(np.square(Glcc[:,0] - Nlcc[idxs,0]) + np.square(Glcc[:,1] - Nlcc[idxs,1]))
  #g = dist< 10000
  #fix = ~g
  #NMatch[g,:] = Npts[idxs[g],:]
  
  aNpts = np.c_[dlon,dlat]
  aN_shape = [shl.Point(s) for s in aNpts ]
  gdfaNpts= geopandas.GeoDataFrame({ 'geometry':aN_shape,'meanQ':meandis}, crs="EPSG:4326")
  gdfaNpts_paclcc = gdfaNpts.to_crs(prjstr)
  aNlcc = np.c_[gdfaNpts_paclcc.geometry.x,gdfaNpts_paclcc.geometry.y]
  fix = np.ones(Glcc.shape[0],dtype=bool)
  unusedN = np.ones(aNpts.shape[0],dtype=bool)
  NMatch = Gpts *np.nan
  GReach = Gpts[:,0] *np.nan
  Nind = np.array(Gpts[:,0]*0,dtype=int)
  cutoff = 1000
  maxcutoff = 10000
  f,ax = plt.subplots()
  xl = ax.get_xlim()
  yl = ax.get_ylim()
  ax.set_xlim(xl)
  ax.set_ylim(yl)
  ax.plot(Gpts[fix,0],Gpts[fix,1],'ok')
  ax.plot(aNpts[:,0],aNpts[:,1],'.c',zorder=1)
  for i in range(len(Gpts)): ax.annotate(f'{GEwest[i]+1}',Gpts[i,:],c='brown')
  ax.axis('equal')
  hgrid.plot_bnd(ax=ax)
  meanbounds = {'conus':1,'alaska':10,'hawaii':1}
  maxbounds = {'conus':1,'alaska':1,'hawaii':1}
  while np.any(fix) and cutoff >= meanbounds[region] and maxcutoff >= maxbounds[region]:
    gn = np.where(np.logical_and(np.logical_or(meandis >= cutoff, maxdis > maxcutoff),unusedN))[0]
    if len(gn)  > 0:
      ax.scatter(dlon[gn],dlat[gn],(np.log10(maxdis[gn])+1)*20,np.log10(meandis[gn]),cmap='cool',vmin=0,vmax=5)
      for i in range(len(gn)): ax.annotate(f'{reaches[gn[i]]}\nmean {meandis[gn[i]]:4.0f} max {maxdis[gn[i]]:4.0f}',(dlon[gn[i]],dlat[gn[i]]))
      unusedN[gn] = False
      ff = np.where(fix)[0]
      tree = KDTree(Glcc[ff,:])
      Ntree = KDTree(aNlcc[gn,:])
      _, idxs = tree.query(np.vstack(aNlcc[gn,:]), workers=-1) # for each NWM reach in this discharge range, find closest GloFAS based element
      gg = ff[idxs]
      dist = np.sqrt(np.square(Glcc[gg,0] - aNlcc[gn,0]) + np.square(Glcc[gg,1] - aNlcc[gn,1])) # find distance between NWM reaches and GloFAS element
      g = np.where(dist< distcut[region])[0]  # exclude distant reaches
      if len(g)==0:
          continue
      ind = np.unique(gg[g])  # find elements that will be assigned NWM reaches
      for i in ind: # for each element to be assigned find closest NWM reach of NWM reaches that are closer to the element than to other elements
          gi = g[np.where(gg[g]==i)[0]]
          shortest = gn[gi[np.argmin(dist[gi])]]
          NMatch[i,:] = aNpts[shortest,:]
          GReach[i] = reaches[shortest]
          Nind[i] = shortest
      ax.set_title('cutoff = %d' % cutoff)
      ax.plot(np.c_[Gpts[ind,0],NMatch[ind,0]].T,np.c_[Gpts[ind,1],NMatch[ind,1]].T,'m')
      fix[ind] = False
    cutoff *= 0.8
    maxcutoff *= 0.8
    print(cutoff)
  ax.plot(Gpts[fix,0],Gpts[fix,1],'py')
  
  f,ax = plt.subplots()
  ax.plot(Gpts[:,0],Gpts[:,1],'o',markersize=15)
  E2Reach = {}
  for i in range(len(GReach)):
    if not np.isnan(GReach[i]):
      E2Reach['%d' % (GEwest[i]+1)]= [GReach[i]]
      plt.annotate(f'{GReach[i]}',(xcG[i],ycG[i]))
  ax.plot(xcG,ycG,'+',markersize=15)
  with open(f'NWM_GLOFAS_{region}.json', 'w') as f:
      json.dump(E2Reach, f, indent=4)
