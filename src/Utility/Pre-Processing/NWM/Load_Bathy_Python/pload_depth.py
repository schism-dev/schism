#!/usr/bin/env python3
'''
load bathymetry for NWM model grid 
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
grd='hgrid.ll'  #grid name
grdout='hgrid.ll.new' #grid name with depth loaded

#parameter
regions=("min_5m_ll.reg","SabinePass.reg","BergenPoint.reg","Washington_3.reg") #regions for modifying depth
rvalues=(5,7,5,15) #minimum depth in regions
headers=("etopo1","crm_3arcs","cdem13_","dem_continetalus_southcarolina","North_Carolina_USGS_3m",
         "al_ll","nc_ll","fl_ll","gulf_1_dem_usgs","gulf_3_demcombined_ll","ge_ll","sc_ll",
         "cb_ll","db_ll","new_england_topobathy_dem_3m_dd","Tile3_R3_DEMv2")
sdir=r'/sciclone/data10/wangzg/DEM/npz'  #directory of DEM data
reverse_sign=1  #invert depth sign

#resource requst 
walltime='00:10:00'
qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=2; ppn=12    #vortex, ppn=12
#qnode='femto'; nnode=2; ppn=12     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='skylake'; nnode=2; ppn=36   #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2    #viz3,haswell, ppn=24,or 28

jname='load_dem' #job name
ibatch=1; scrout='screen.out'; bdir=os.path.abspath(os.path.curdir)
#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt=fmt)
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
#get all DEM files and distribute jobs
fnames0=array([i for i in os.listdir(sdir) if i.endswith('.npz')])

#filter with headers, and sort by id numbers
fnames_sort=[]
for header in headers:
    fnames_sub=array([i for i in fnames0 if i.startswith(header)])
    if len(fnames_sub)==1: fnames_sort.extend(fnames_sub); continue

    #get id number 
    fid=array([i.replace('tif.','').replace('.','_')[len(header):].split('_')[-2] for i in fnames_sub]).astype('int')
    sind=argsort(fid); fnames_sub=fnames_sub[sind]; fnames_sort.extend(fnames_sub)
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
S=npz_data(); S.dp=dict(); S.sind=dict()
for m,fname in enumerate(fnames):
    bname=fname.split('.')[0]

    #interpolate depth
    while(True):
        try:
           dpi,sindi=load_bathymetry(gd.x,gd.y,'{}/{}'.format(sdir,fname),fmt=1)
           break
        except:
            time.sleep(15)

    #save results
    S.dp[bname]=dpi; S.sind[bname]=sindi
    print('finished reading {},: {}, myrank={}'.format(fname,inum[m],myrank)); sys.stdout.flush()
#save_npz('S_{}'.format(myrank),S)

#gather results
comm.Barrier()
sdata=comm.gather(S,root=0)

if myrank==0:
   #combine
   S=npz_data(); S.dp=dict(); S.sind=dict()
   for i in arange(nproc):
       Si=sdata[i]
       S.dp={**S.dp,**Si.dp}
       S.sind={**S.sind,**Si.sind}

   #load bathymetry
   for i,fname in enumerate(fnames_sort):
       bname=fname.split('.')[0]
       sind=S.sind[bname]; dp=S.dp[bname]
       gd.dp[sind]=dp

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
      S=npz_data(); S.hgrid=gd; save_npz(grdout,S)
   else:
      gd.write_hgrid(grdout)
   
#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora'] else os._exit(0)
