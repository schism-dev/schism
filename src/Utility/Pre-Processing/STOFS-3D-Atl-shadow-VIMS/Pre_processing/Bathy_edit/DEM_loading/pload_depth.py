#!/usr/bin/env python3
'''
load bathymetry for NWM model grid
'''
import os
import re
from pylib import *
import time
import errno
import json
from glob import glob

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
grd='hgrid.ll'  #grid name
grdout='hgrid.ll.dem_loaded.gr3' #grid name with depth loaded

#parameter
regions=[
    '/sciclone/schism10/Hgrid_projects/DEMs/regions/'+i for i in ("min_5m_ll_noPR.reg","SabinePass.reg","BergenPoint.reg","Washington_3.reg",
    "Elk_river.reg","Hudson_river.reg","James_river.reg","NorthEast_river.reg",
    "Rappahannock_river.reg","Susquehanna_river.reg","York_river.reg",
    "Androscoggin_Kennebec_rivers.reg","Merrimack_river.reg","Patuxent_river.reg",
    "Penobscot_river.reg","Saco_river.reg","StCroix_river.reg","Oyster_landing.reg",
    "st_lawrence1.reg", "st_lawrence2.reg", "st_lawrence3.reg",
)] #regions for modifying depth
rvalues=(
    5,7,5,15,
    2,16,14,5,
    6,10,10,
    3,3,5,
    5,3,5,1,
    10,10,10
) #minimum depth in regions

dem_info_dict = json.load(open('DEM_info.json'))  # existing json containing the DEM info
headers=tuple(dem_info_dict.keys()) #DEM headers used for glob'ing the DEM files; python3.6+ only, which preserves the order of the keys

sdir=r'/sciclone/schism10/Hgrid_projects/DEMs/npz2/'  #directory of DEM data
reverse_sign=1  #invert depth sign

#resource requst
walltime='01:30:00'
#qnode='x5672'; nnode=4; ppn=4       #hurricane, ppn=8
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
qnode='vortex'; nnode=6; ppn=12    #vortex, ppn=12
#qnode='femto'; nnode=1; ppn=32     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='frontera'; nnode=1; ppn=56  #frontera, ppn=56 (flex,normal)
#qnode='mistral'; nnode=1; ppn=36   #mistral, ppn=36
#qnode='stampede2'; nnode=1; ppn=48 #stampede2, ppn=48 (skx-normal,skx-dev,normal,etc)

#additional information:  frontera,mistral,stampede2
qname='flex'                        #partition name
account='TG-OCE140024'              #stampede2: NOAA_CSDL_NWI,TG-OCE140024; mistral: gg0028

jname='load_dem' #job name
ibatch=0; scrout='screen.out'; bdir=os.path.abspath(os.path.curdir)
#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun)
#-----------------------------------------------------------------------------
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt=fmt,qname=qname)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
bdir=os.getenv('bdir'); os.chdir(bdir) #enter working dir
comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
#check existence of region files
if myrank==0:
    for region in regions:
        if not os.path.exists(region):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), region)

#get all DEM files and distribute jobs
fnames0=array([i for i in os.listdir(sdir) if i.endswith('.npz')])
#check existence of DEM tiles
if myrank==0:
  for header in headers:
      fnames_sub=array([i for i in fnames0 if i.startswith(header)])
      if len(fnames_sub)==0:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f"{header}*")

#filter with headers, and sort by id numbers
fnames_sort=[]
for header in headers:
    fnames_sub=array([fname for fname in fnames0 if fname.startswith(header)])
    if len(fnames_sub)==1: fnames_sort.extend(fnames_sub); continue

    #get id number, which is the last number before the suffix ".*"
    fid = array([int(re.findall(r'(\d+)', i)[-1]) for i in fnames_sub])
    sind=argsort(fid); fnames_sub=fnames_sub[sind]; fnames_sort.extend(fnames_sub)
    pass
fnames_sort=array(fnames_sort)

#to exactly match the order to old method
#switch the order of (db_ll_1.npz and db_ll1.npz), ('cdem13_DE_navd88_2016' and 'cdem13_PA_navd88_2010')
#dn1='db_ll1.npz'; dn2='db_ll_1.npz'
#if (dn1 in fnames_sort)*(dn2 in fnames_sort):
#    fid1=nonzero(fnames_sort==dn1)[0][0]; fid2=nonzero(fnames_sort==dn2)[0][0]
#    fnames_sort[fid1]=dn2; fnames_sort[fid2]=dn1
#dn1='cdem13_PA_navd88_2010.npz'; dn2='cdem13_DE_navd88_2016.npz'
#if (dn1 in fnames_sort)*(dn2 in fnames_sort):
#    fid1=nonzero(fnames_sort==dn1)[0][0]; fid2=nonzero(fnames_sort==dn2)[0][0]
#    fnames_sort[fid1]=dn2; fnames_sort[fid2]=dn1

#distribute jobs
fnames=[]; inum=[]
for m,fname in enumerate(fnames_sort):
    if m%nproc==myrank: fnames.append(fname); inum.append(m)
fnames=array(fnames)

#read hgrid
if grd.endswith('npz'):
   gd=loadz(grd).hgrid
elif grd.endswith('gr3') or grd.endswith('ll'):
   gd=read_schism_hgrid(grd)
else:
   sys.exit('wrong format of grid: {}'.format(grd)); sys.stdout.flush()

#load bathymetry on each core
S=zdata(); S.dp=dict(); S.sind=dict()
for m,fname in enumerate(fnames):
   bname=fname.split('.')[0]

   #interpolate depth
   while(True):
       try:
          dpi,sindi=load_bathymetry(gd.x,gd.y,'{}/{}'.format(sdir,fname),fmt=1)
          break
       except:
           time.sleep(15)

   # reverse depth sign
   if reverse_sign==1:
      dpi=-dpi

   #save results
   S.dp[bname]=dpi; S.sind[bname]=sindi
   print('finished reading {},: {}, myrank={}'.format(fname,inum[m],myrank)); sys.stdout.flush()
#savez('S_{}'.format(myrank),S)

#gather results
comm.Barrier()
sdata=comm.gather(S,root=0)

if myrank==0:
   #combine
   S=zdata(); S.dp=dict(); S.sind=dict()
   for i in arange(nproc):
       Si=sdata[i]
       S.dp={**S.dp,**Si.dp}
       S.sind={**S.sind,**Si.sind}

   #load bathymetry
   did=ones(gd.np).astype('int'); dname=[]
   for i,fname in enumerate(fnames_sort):
       bname=fname.split('.')[0]
       sind=S.sind[bname]; dp=S.dp[bname]
       gd.dp[sind]=dp; did[sind]=i+1
       dnamei=[k for k in fnames0 if k.startswith(fname)][0]; dname.append(dnamei)

   #reverse depth sign
   if reverse_sign==1:
      gd.dp=-gd.dp

   #applying minimum depth
   if regions is not None:
      for i, region in enumerate(regions):
          if not os.path.exists(region): continue
          depth_min=rvalues[i]
          bp=read_schism_bpfile(region,fmt=1)
          sind=inside_polygon(c_[gd.x,gd.y], bp.x,bp.y)
          fp=(sind==1)*(gd.dp<depth_min); gd.dp[fp]=depth_min
          print('finished applying min depth={}: {}'.format(depth_min,region)); sys.stdout.flush()

   #save grid
   if grdout.endswith('npz'):
      S=zdata(); S.hgrid=gd; savez(grdout,S)
   else:
      gd.write_hgrid(grdout)
   fid=open('{}_dem_id'.format(grdout),'w+'); [fid.write('{}\n'.format(i)) for i in did]; fid.close()
   fid=open('{}_dem_name'.format(grdout),'w+'); [fid.write('{}: {}\n'.format(i+1,k)) for i,k in enumerate(dname)]; fid.close()

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora'] else os._exit(0)
