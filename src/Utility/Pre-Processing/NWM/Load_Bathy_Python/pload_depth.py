#!/usr/bin/env python3
'''
load bathymetry for NWM model grid 
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
#grd='/sciclone/data10/wangzg/NWM/grid/v7_3D/NWM_3D_v7.ll'  #grid name
#grdout='hgrid.ll' #grid name with depth loaded
grd='/sciclone/data10/wangzg/NWM/grid/v17/NWM_2D_v17.ll'  #grid name
grdout='/sciclone/data10/wangzg/NWM/grid/v17/hgrid.ll' #grid name with depth loaded

#parameter
regions=("min_5m_ll.reg","SabinePass.reg","BergenPoint.reg","Washington_3.reg") #regions for modifying depth
rvalues=(5,7,5,15) #minimum depth in regions
headers=("etopo1","crm_3arcs","cdem13_","dem_continetalus_southcarolina","North_Carolina_USGS_3m",
         "al_ll","nc_ll","fl_ll","gulf_1_dem_usgs","gulf_3_demcombined_ll","ge_ll","sc_ll",
         "cb_ll","db_ll","new_england_topobathy_dem_3m_dd","Tile3_R3_DEMv2")
sdir=r'/sciclone/data10/wangzg/DEM/npz'  #directory of DEM data
reverse_sign=1  #invert depth sign

#job name and time
jname='load_dem' #job name
walltime='01:00:00' 

#resource requst 
#qnode='bora'; nnode=2; ppn=5      #bora, ppn=20
#qnode='vortex'; nnode=2; ppn=4   #vortex, ppn=12
qnode='x5672'; nnode=20; ppn=4     #hurricane, ppn=8
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='femto'; nnode=1; ppn=2      #femto,ppn=32, not working yet

#-----------------------------------------------------------------------------
#pre-processing
#-----------------------------------------------------------------------------
nproc=nnode*ppn
bdir=os.path.abspath(os.path.curdir)

#-----------------------------------------------------------------------------
#on front node; submit jobs in this section
#-----------------------------------------------------------------------------
if os.getenv('param')==None and os.getenv('job_on_node')==None:
    args=sys.argv
    param=[bdir,args[0]]
    
    #submit job on node
    if qnode=='femto': 
        scode='sbatch --export=param="{} {}" -J {} -N {} -n {} -t {} {}'.format(*param,jname,nnode,nproc,walltime,args[0])
    else:
        scode='qsub {} -v param="{} {}", -N {} -j oe -l nodes={}:{}:ppn={} -l walltime={}'.format(args[0],*param,jname,nnode,qnode,ppn,walltime)
    print(scode); os.system(scode)
    os._exit(0)

#-----------------------------------------------------------------------------
#still on front node, but in batch mode; running jobs in this section
#-----------------------------------------------------------------------------
if os.getenv('param')!=None and os.getenv('job_on_node')==None:
    param=os.getenv('param').split();
    param=[int(i) if i.isdigit() else i for i in param]
    bdir=param[0]; bcode=param[1]
    os.chdir(bdir)

    if qnode=='bora':
       rcode="mpiexec -x job_on_node=1 -x bdir='{}' -n {} {} >& screen.out".format(bdir,nproc,bcode)
    elif qnode=='femto':
       pypath='/sciclone/home10/wangzg/bin/pylibs/Scripts/:/sciclone/home10/wangzg/bin/pylibs/Utility/'
       rcode="srun --export=job_on_node=1,bdir='{}',PYTHONPATH='{}' {} >& screen.out".format(bdir,pypath,bcode)
    elif qnode=='x5672' or qnode=='vortex' or qnode=='potomac' or qnode=='james':
       rcode="mvp2run -v -e job_on_node=1 -e bdir='{}' {} >& screen.out".format(bdir,bcode)
    print(rcode); os.system(rcode); sys.stdout.flush()
    os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
#enter working dir
bdir=os.getenv('bdir'); os.chdir(bdir)

#get nproc and myrank
comm=MPI.COMM_WORLD
nproc=comm.Get_size()
myrank=comm.Get_rank()
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
    fid=array([i.replace('.','_')[len(header):].split('_')[-2] for i in fnames_sub]).astype('int')
    sind=argsort(fid); fnames_sub=fnames_sub[sind]; fnames_sort.extend(fnames_sub)
fnames_sort=array(fnames_sort)

#to exactly match the order to old method
#switch the order of (db_ll_1.npz and db_ll1.npz), ('cdem13_DE_navd88_2016' and 'cdem13_PA_navd88_2010')
dn1='db_ll1.npz'; dn2='db_ll_1.npz'
if (dn1 in fnames_sort)*(dn2 in fnames_sort):
    fid1=nonzero(fnames_sort==dn1)[0][0]; fid2=nonzero(fnames_sort==dn2)[0][0]
    fnames_sort[fid1]=dn2; fnames_sort[fid2]=dn1
dn1='cdem13_PA_navd88_2010.npz'; dn2='cdem13_DE_navd88_2016.npz'
if (dn1 in fnames_sort)*(dn2 in fnames_sort):
    fid1=nonzero(fnames_sort==dn1)[0][0]; fid2=nonzero(fnames_sort==dn2)[0][0]
    fnames_sort[fid1]=dn2; fnames_sort[fid2]=dn1

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
if qnode=='x5672' or qnode=='james':
   os._exit(0)
else:
   sys.exit(0)
