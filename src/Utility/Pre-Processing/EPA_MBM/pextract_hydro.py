#!/usr/bin/env python3
'''
  extract time series of points@xyz or transects@xy from SCHISM outputs
'''
from pylib import *
from mpi4py import MPI

#-----------------------------------------------------------------------------
#Input
#hpc: kuro, femto, bora, potomac, james, frontera, levante, stampede2
#ppn:  64,   32,    20,    12,     20,     56,      128,      48
#-----------------------------------------------------------------------------
run='/sciclone/data10/wangzg/CBP/RUN09c'

#input for different modules: flag,svars, bpfile,fname,itype
dir0='/sciclone/data10/wangzg/CBP/Results/Station/'
inputs=(
  (1,['salt','temp','salt_elem','temp_elem'],dir0+'station_salt.bp','salt_temp',0), #salt@xyz
  (1,['salt','salt_elem','elev'],dir0+'station_channel.bp','salt_profile',1), #salt profile@xy
  (1,['elev'],dir0+'station_elev.bp','elev',0), #elev @xy
  (1,['sed_tconc','sed_dp','sed_str','sed_rough','sed_eflux', 'sed_dflux'],dir0+'station_sed.bp','sed_profile',1), #sediment profile @xy
  (1,['wwm_WVHT','wwm_DPD','wwm_Tp'],dir0+'station_wave.bp','wave',0), #elev @xy
  (0,['sed_conc1','sed_conc2','sed_conc3','sed_conc4','sed_tconc','sed_dp','sed_str','sed_rough','sed_eflux',
      'sed_dflux','sed_frac1','sed_frac2','sed_frac3','sed_frac4',],dir0+'station_sed.bp','sed_profile_ext',1), #detailed sediment profile @xy
)

#optional
#itype=1         #0: time series of points @xyz;  1: time series of trasects @xy
#ifs=0           #0: refer to free surface; 1: fixed depth
#stacks=[1,3]    #outputs stacks to be extracted
#nspool=12       #sub-sampling frequency within each stack (1 means all)
#mdt=1           #time window (day) for averaging output
#rvars=['elev','salt','hvel','NO3'] #rname the varibles 
#prj=['epsg:26918','epsg:4326']  #projections used to transform coord. in station.bp

#hpc resource requst
walltime='24:00:00'; nnode=1;  ppn=32

#optional: (frontera,levante,stampede2,etc.)
ibatch     =1              #0: serial mode;  1: parallel mode
qnode      =None           #specify node name, or default qnode based on HOST will be used
qname      =None           #partition name
account    =None           #account name
reservation=None           #reservation information

#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
brun=os.path.basename(run); jname='Rd_'+brun; scrout='screen_{}.out'.format(brun); bdir=os.path.abspath(os.path.curdir)
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]; os.environ['qnode']=get_qnode(qnode)
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt,'param',qname,account,reservation)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
bdir=os.getenv('bdir'); os.chdir(bdir) #enter working dir
odir=os.path.abspath('.')
if ibatch==0: nproc=1; myrank=0
if ibatch==1: comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()
if myrank==0 and (not fexist(odir)): os.mkdir(odir)

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
sdir=run+'/outputs'                                              #output directory
if 'ifs' not in locals(): ifs=0                                  #refer to free surface
if 'nspool' not in locals(): nspool=1                            #subsample
if 'prj' not in locals(): prj=None                               #projections
if 'mdt' not in locals(): mdt=None                               #average
modules, outfmt, dstacks, dvars, dvars_2d = get_schism_output_info(sdir,1)     #schism outputs information
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks #check stacks

#read model grid
fgz=run+'/grid.npz'; fgd=run+'/hgrid.gr3'; fvd=run+'/vgrid.in'
gd=loadz(fgz,'hgrid') if fexist(fgz) else read_schism_hgrid(fgd); gd.compute_bnd()
vd=loadz(fgz,'vgrid') if fexist(fgz) else read_schism_vgrid(fvd); sys.stdout.flush()

#extract results for each module
for nn, [iflag,svars,bpfile,sname,itype] in enumerate(inputs):
    if iflag==0: continue
    sname=odir+'/'+sname; rvars=svars

    #extract results for different bpfiles
    irec=0; oname=odir+'/.schout'
    for svar in svars: 
       ovars=get_schism_var_info(svar,modules,fmt=outfmt)
       if ovars[0][1] not in dvars: continue 
       for istack in stacks:
           fname='{}_{}_{}'.format(oname,svar,istack); irec=irec+1; t00=time.time()
           if irec%nproc==myrank: 
              try:
                 read_schism_output(run,svar,bpfile,istack,ifs,nspool,fname=fname,hgrid=gd,vgrid=vd,fmt=itype,prj=prj,mdt=mdt)
                 dt=time.time()-t00; print('finishing reading {}_{}.nc on myrank={}: {:.2f}s'.format(svar,istack,myrank,dt)); sys.stdout.flush()
              except:
                 pass
    
    #combine results
    if ibatch==1: comm.Barrier()
    if myrank==0:
       S=zdata(); S.time=[]; fnames=[]
       for i,[k,m] in enumerate(zip(svars,rvars)):
           data=[]; mtime=[]
           for istack in stacks:
               fname='{}_{}_{}.npz'.format(oname,k,istack)
               if not fexist(fname): continue
               C=loadz(fname); datai=C.__dict__[k]; fnames.append(fname)
               data.extend(datai.transpose([1,0,*arange(2,datai.ndim)])); mtime.extend(C.time)
           if len(data)>0: S.__dict__[m]=array(data).transpose([1,0,*arange(2,array(data).ndim)])
           if len(mtime)>len(S.time): S.time=array(mtime)
       S.bp=read_schism_bpfile(bpfile)
       for pn in ['param','icm','sediment','cosine','wwminput']:
           if fexist('{}/{}.nml'.format(run,pn)): S.__dict__[pn]=read_schism_param('{}/{}.nml'.format(run,pn),3)
       savez(sname,S)
       for i in fnames: os.remove(i)
    if ibatch==1: comm.Barrier()

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
