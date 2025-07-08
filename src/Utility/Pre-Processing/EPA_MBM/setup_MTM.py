#!/usr/bin/env python3
#written by Zhengui Wang on 2/1/2024
#Note: This script and associated databases are for internal use within VIMS MTM team (Rappahannock River) 
#      and not to be shared without the constent from all authors
from pylib import *
import subprocess as sbp
from shutil import copyfile
from glob import glob
p=zdata(); p.flag={}

#----------------------------------------------------------------------
#Inputs:
#  base:  reference run for current setup
#  flag=0:
#       if base is not None:  create link to the input file in base run
#       if base is None: skip this file
#  flag=1: re-generate input file
#----------------------------------------------------------------------
p.StartT=datenum(1991,1,1);  p.EndT=datenum(2000,12,31) #simulation time

p.base=None
p.grid_dir='/sciclone/data10/wangzg/MTM/grid/v0' #directory of hgrid & vgrid; p.grid_dir=p.base if it is None

p.flag['elev2D.th.nc']         = 0 #hydro
p.flag['TEM_3D.th.nc']         = 0 #hydro
p.flag['TEM_nu.nc']            = 0 #hydro
p.flag['SAL_3D.th.nc']         = 0 #hydro
p.flag['SAL_nu.nc']            = 0 #hydro
p.flag['uv3D.th.nc']           = 0 #hydro
p.flag['sflux']                = 0 #hydro
p.flag['albedo.gr3']           = 0 #hydro
p.flag['drag.gr3']             = 0 #hydro
p.flag['watertype.gr3']        = 0 #hydro
p.flag['diffmin.gr3']          = 0 #hydro
p.flag['diffmax.gr3']          = 0 #hydro
p.flag['shapiro.gr3']          = 0 #hydro
p.flag['TEM_nudge.gr3']        = 0 #hydro
p.flag['SAL_nudge.gr3']        = 0 #hydro
p.flag['tvd.prop']             = 0 #hydro
p.flag['windrot_geo2proj.gr3'] = 0 #hydro
p.flag['bctides.in']           = 0 #hydro
p.flag['hotstart.nc']          = 0 #hydro
p.flag['source.nc']            = 0 #hydro

#databases
p.bdir='/sciclone/data10/wangzg/MTM/database/'         #dir of setup files
p.MBM       = p.bdir+'hydro_MBM/RUN07b/'               #Main Bay Model hydro outputs
p.sflux     = p.bdir+'sflux_narr_subdomain'            #sflux database
p.source    = p.bdir+'load_p6_v3.npz'                  #CBP watershed sources
p.region    = p.bdir+'region/'                         #regions files
p.outdir    = '/sciclone/pscr/{}/MTM'.format(os.environ['USER']) #parental direcotry of outputs
p.bndfile   = None                                     #boundary conditions from MBM  

#----------------------------------------------------------------------
#parallel processing to speed up
#----------------------------------------------------------------------
try:
   comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
except:
   nproc=1; myrank=0

#----------------------------------------------------------------------
#model configuration (don't change this line)
#----------------------------------------------------------------------
if myrank==0:
   #check reference run
   if (p.base is not None) and (not fexist(p.base)):
      sys.exit('reference run not exist: {}'.format(p.base))
   
   #reset values
   if p.grid_dir is None: p.grid_dir=p.base
   
   #define function to make symbolic links to the files from base, and  check each input file, and make links
   def linkfile(fname):
       sname='{}/{}'.format(p.base,fname)
       if fexist(sname) and (not fexist(fname)):
          os.symlink(os.path.relpath(os.path.realpath(sname)),fname)
   
   for fname in p.flag:
       if (p.base is not None) and p.flag[fname]==0 : linkfile(fname)
       if p.flag[fname]==1 and fexist(fname): 
          if os.path.isfile(fname):
             os.remove(fname)
          else:
             os.system('rm -rf {}'.format(fname))
   
   #create outputs directory
   if not fexist('outputs'):
      outdir='{}/{}/outputs'.format(p.outdir,os.path.basename(os.path.abspath('.')))
      try:
         os.system('mkdir -p {}'.format(outdir))
      except:
         pass
      if not fexist(outdir): sys.exit('check p.outdir path: {}'.format(outdir))
      os.symlink(outdir,'outputs')
   
#----------------------------------------------------------------------
#hgrid.gr3, hgrid.ll, vgrid.in, grid.npz
#----------------------------------------------------------------------
for fname in ['hgrid.gr3','hgrid.ll','vgrid.in','grid.npz']:
    if myrank==0:
       print('writing '+fname)
       sname=os.path.relpath(os.path.realpath('{}/{}'.format(p.grid_dir,fname)))
       if not fexist(sname): sys.exit('{}/{} not exist'.format(p.grid_dir,fname))
       if os.path.islink(fname) or (fexist(fname) and sname!=fname): os.remove(fname)
       if sname!=fname: os.symlink(sname,fname)
if nproc>1: comm.Barrier()
p.grd='grid.npz'; gd=read(p.grd,'hgrid')

#----------------------------------------------------------------------
#parameter files
#----------------------------------------------------------------------
fname='param.nml'; sname='{}/{}'.format(p.base,fname)
if not fexist(fname) and myrank==0:
   if fexist(sname):
      copyfile(sname,fname); print('copy {} from {}'.format(fname,p.base))
   else:
      copyfile(p.bdir+'param/'+fname,fname); print('writing {}'.format(fname))

#change parameters
fname='param.nml'; sname='{}/{}'.format(p.base,fname)
if not fexist(sname) and myrank==0:
   chparam(fname,'start_year',num2date(p.StartT).year)
   chparam(fname,'start_month',num2date(p.StartT).month)
   chparam(fname,'start_day',num2date(p.StartT).day)
   chparam(fname,'rnday',int(p.EndT-p.StartT))

#----------------------------------------------------------------------
#bctides.in
#----------------------------------------------------------------------
fname='bctides.in'
if p.flag[fname]==1  and myrank==0:
   print('writing '+fname)
   fid=open('./bctides.in','w+')
   fid.write('{} GMT\n'.format(num2date(p.StartT).strftime('%m/%d/%Y %H:%M:%S')))
   lstr='0 40. ntip\n0 nbfr\n1 nope\n{:d} 4 0 4 4'.format(gd.nobn[0])
   bstr='\n1.0 !TEM nudge\n1.0 !SAL nudge'
   fid.write(lstr+bstr); fid.close()

#----------------------------------------------------------------------
#gr3,prop files
#----------------------------------------------------------------------
#constant *.gr3 file
fnames=('albedo.gr3','drag.gr3','watertype.gr3','windrot_geo2proj.gr3','diffmin.gr3', 'diffmax.gr3',)
fvalues=(0.1,   2.5e-3,   4,     0.0,    1e-6, 0.1,)
for fname,fvalue in zip(fnames,fvalues):
    if p.flag[fname]==1 and myrank==0:
       print('writing '+fname)
       gd.save(fname,value=fvalue)

#shapiro.gr3
fname='shapiro.gr3'
if p.flag[fname]==1 and myrank==0:
   print('writing '+fname)
   shapiro_max,threshold_slope=0.5,0.5
   slope=gd.compute_gradient(fmt=2)[2]; fvi=shapiro_max*tanh(2*slope/threshold_slope)
   gd.save('shapiro.gr3',value=fvi)

#TEM_nudge.gr3, SAL_nudge.gr3
for fname,sname in zip(['TEM_nudge.gr3','SAL_nudge.gr3',],['TEM_3D.th.nc','SAL_3D.th.nc']):
    if (p.flag[fname]==1 or (p.flag[sname]==1 and (not fexist(fname)))) and myrank==0:
       print('writing '+fname)
       dist=abs((gd.x+1j*gd.y)[:,None]-(gd.x[gd.iobn[0]]+1j*gd.y[gd.iobn[0]])[None,:]).min(axis=1)
       x1,x2= [2e3,5e3]
       vi=1.15740741e-05; fvi=vi*(dist-x2)/(x1-x2); fvi[fvi<0]=0; fvi[fvi>vi]=vi
       gd.save(fname,value=fvi,outfmt='{:16.8e}')

#tvd.prop
fname='tvd.prop'
if p.flag[fname]==1 and myrank==0:
   print('writing '+fname)
   gd.compute_ctr(); sindp=read(p.region+'stream_head.reg').inside(c_[gd.xctr,gd.yctr])
   fvi=ones(gd.ne); fvi[sindp]=0; gd.write_prop('tvd.prop',value=fvi,fmt='{:d}')

#----------------------------------------------------------------------
#parse MBM setup
#----------------------------------------------------------------------
#p.MBM='/sciclone/data10/wangzg/CBP/RUN07a/'

#read grid info
grd0=p.MBM+'grid.npz';  gd0=read(grd0,'hgrid'); vd0=read(grd0,'vgrid')
gd=read(p.grd,'hgrid'); vd=read(p.grd,'vgrid')

MBM=get_schism_output_info(p.MBM+'outputs',3); [MBM.modules,MBM.outfmt,stacks,MBM.svars,MBM.svars_2d]=get_schism_output_info(p.MBM+'outputs')
if p.StartT<(MBM.StartT-MBM.dt-1.0) or p.EndT>MBM.EndT:
   sys.exit('time range ({},{}) beyond MBM outputs ({},{})'.format(*num2date([p.StartT,p.EndT,MBM.StartT,MBM.EndT])))

#----------------------------------------------------------------------
#hotstart.nc
#----------------------------------------------------------------------
fname='hotstart.nc'
if p.flag[fname]==1 and myrank==0:
   print('writing '+fname)

   #grid info
   np,ne,ns,nvrt=gd.np,gd.ne,gd.ns,vd.nvrt
   pxy=c_[gd.x,gd.y]; pie,pip,pacor=gd0.compute_acor(pxy) #interp indices
   fp=pie==-1; pie[fp]=near_pts(pxy[fp],c_[gd0.xctr,gd0.yctr]); pip[fp,0]=near_pts(pxy[fp],c_[gd0.x,gd0.y])

   #write hotstart.nc template
   from netCDF4 import Dataset
   fid=Dataset(fname,'w',format='NETCDF4'); fid.setncattr('file_format','NETCDF4')
   dnames=['node','elem','side','nVert','ntracers','one']; dims=[np,ne,ns,nvrt,2,1]
   vnames=(*['time','iths','ifile','nsteps_from_cold'], *['idry','eta2','cumsum_eta'], *['q2','xl','dfv','dfh','dfq1','dfq2'],
           *['su2','sv2'], *['tr_nd','tr_nd0'], 'idry_e', 'idry_s', 'we', 'tr_el')
   vdims =(*['one']*4, *['node']*3, *[['node','nVert']]*6, 
           *[['side','nVert']]*2, *[['node','nVert','ntracers']]*2 , 'elem', 'side', ['elem','nVert'], ['elem','nVert','ntracers'])
   vtypes=(*[2,1,1,1], *[1,2,2],  *[2,2,2,2,2,2], *[2,2],*[2,2], 1,1,2,2); vtypes=['int32' if i==1 else 'float' for i in vtypes]
   for dn,ds in zip(dnames,dims): fid.createDimension(dn,ds)  #define dimensions
   for vn,vdm,vt in zip(vnames,vdims,vtypes): #define variables
       vds=[vdm] if isinstance(vdm,str) else vdm;  vdi=[dims[dnames.index(i)] for i in vds] #get dimension
       fid.createVariable(vn,dtype(vt),vdm,zlib=True); fid.variables[vn][:]=(zeros(vdi) if vn!='ifile' else ones(vdi)).astype(vt) 
   fid.close()
   
   #find the outputs stack and record     
   sid=argmin(abs(p.StartT-MBM.mts)); istack,irec=MBM.istack[sid],MBM.irec[sid]; fid=read(fname,1,mode='r+') 
  
   #elevation
   if 'elevation' in MBM.svars:
      C=read(MBM.outdir+'/out2d_{}.nc'.format(istack),1)
      fid.variables['eta2'][:]=(array(C.variables['elevation'][irec])[pip]*pacor).sum(axis=1); C.close()

   #interpolation for 3D variables (ST)
   zs=(vd0.compute_zcor(gd0.dp)[pip]*pacor[...,None]).sum(axis=1).T; zcor=vd.compute_zcor(gd.dp).T; ds=zcor.shape
   svars=['temperature','temperatureAtElement','salinity','salinityAtElement']
   for m,svar in enumerate(svars):
      if (svar not in MBM.svars) or (m==1 and (svars[0] in MBM.svars)) or (m==3 and (svars[2] in MBM.svars)): continue
      C=read(MBM.outdir+'/{}_{}.nc'.format(svar,istack),1) #extract values @node or element
      vs=(array(C.variables[svar][irec])[pip]*pacor[...,None]).sum(axis=1).T if (m in [0,2]) else array(C.variables[svar][irec])[pie].T; C.close() 
      for k in arange(vd0.nvrt-1)[::-1]: fpn=abs(vs[k-1])>1e6; vs[k-1,fpn]=vs[k,fpn]  #extend bottom values
      data=interp_vertical(vs,zs,zcor).T #interpolation in the vertical

      #write to hotstart
      idm=0 if m in [0,1] else 1
      fid.variables['tr_nd'][...,idm]=data; fid.variables['tr_nd0'][...,0]=data
      fid.variables['tr_el'][...,idm]=gd.interp_node_to_elem(value=data)
   fid.close()

#----------------------------------------------------------------------
#boundary and nudge conditions: 'elev2D.th.nc','TEM_3D.th.nc','TEM_nu.nc','SAL_3D.th.nc','SAL_nu.nc'
#----------------------------------------------------------------------
fnames=['elev2D.th.nc','uv3D.th.nc','TEM_3D.th.nc','SAL_3D.th.nc','TEM_nu.nc','SAL_nu.nc']; dts=[1/24,1/24,1,1,1,1]; bname='bndfile.npz'

#create database for boundary conditions
if (p.bndfile is None) and (not fexist(bname)): 
   C=zdata(); C.time={}; C.data={}; C.xy={}; C.vars=[]; C.MBM=p.MBM; C.StartT=MBM.StartT
   sgrid=gd0.scatter_to_grid(fmt=2)
   for m,[fname,dt] in enumerate(zip(fnames,dts)):
       if p.flag[fname]==0: continue
   
       #check variables to be extracted
       if fname=='elev2D.th.nc':    mvar='elevation'
       if fname.startswith('TEM'):  mvar='temp' if 'temperature' in MBM.svars else 'temp_elem'
       if fname.startswith('SAL'):  mvar='salt' if 'salinity' in MBM.svars else 'salt_elem'
       if fname.startswith('uv3D'): mvar='hvel' if 'horizontalVelY' in MBM.svars else 'hvel_side'
   
      #pts
       if fname.endswith('nu.nc'):
          gdn=read(fname[:3]+'_nudge.gr3'); sindp=nonzero(gdn.dp!=0)[0]; bxy=c_[gdn.x,gdn.y][gdn.dp!=0]
       else:
          sindp=gd.iobn[0]; bxy=c_[gd.x,gd.y][sindp]
       vname=fname.split('.')[0];  C.xy[vname]=bxy; C.vars.append(vname)
   
       #extract results 
       fp=(MBM.mts>=(p.StartT-MBM.dt))*(MBM.mts<=(p.EndT+MBM.dt)); fs=MBM.istack[fp]; istacks=arange(fs.min(),fs.max()+1)
       for istack in istacks:
           rname='.schout_{}_{}.npz'.format(vname,istack)
           if not fexist(rname) and istack%nproc==myrank:
              print('reading {}, stack={}'.format(rname,istack)); sys.stdout.flush()
              read_schism_output(p.MBM,mvar,bxy,istack,fmt=1,sname='data',fname=rname,nspool=max(1,round(dt/MBM.dt)),sgrid=sgrid)
       
       #collect results
   if nproc>1: comm.Barrier()
   if myrank==0:
      for m,[fname,dt] in enumerate(zip(fnames,dts)):
          mtime,mdata=[],[]; vname=fname.split('.')[0]
          for dn in glob('.schout_{}_*.npz'.format(vname)):
              s=read(dn); mtime.extend(s.time); mdata.extend(s.data.transpose([1,0,*arange(2,s.data.ndim)])); os.remove(dn)
          sind=argsort(mtime); C.time[vname]=array(mtime)[sind]; C.data[vname]=array(mdata)[sind]
      p.bndfile=bname; C.save(bname)

#write boundary conditions
if myrank==0:
   C=read(p.bndfile if (p.bndfile is not None) else bname) 
   for m,[fname,dt] in enumerate(zip(fnames,dts)):
       if p.flag[fname]==0: continue
       print('writing '+fname)

       #check database
       vname=fname.split('.')[0]
       if vname not in C.time: sys.exit('{} database not exist'.format(fname)) 
       mtime=C.StartT+C.time[vname]; mdata=C.data[vname]; bxy0=C.xy[vname]
       if p.StartT<(mtime[0]-3) or p.EndT>(mtime[-1]+3): 
          sys.exit('check time range of database {} [{}, {}]: [{}, {}]'.format(fname, num2date(mtime[0]),num2date(mtime[-1]),num2date(p.StartT),num2date(p.EndT)))
       fpt=(mtime>=(p.StartT-1))*(mtime<=(p.EndT+1)); mtime=mtime[fpt]; mdata=mdata[fpt]

       #bnd pts
       if fname.endswith('nu.nc'):
          gdn=read(fname[:3]+'_nudge.gr3'); sindp=nonzero(gdn.dp!=0)[0]; bxy=c_[gdn.x,gdn.y][gdn.dp!=0]
       else:
          sindp=gd.iobn[0]; bxy=c_[gd.x,gd.y][sindp]
       if len(bxy0)!=len(bxy): sys.exit('number to boundary points not consistent in database {} {}'.format(len(bxy0),len(bxy)))
       dx=abs((bxy0[:,0]-bxy[:,0])+1j*(bxy0[:,1]-bxy[:,1])); nobn=len(sindp)
       if dx.max()>50: sys.exit('the location of {} points not match database'.format(fname))
      
       #interp in the vertical and in time 
       t0=time.time()
       if fname!='elev2D.th.nc':
          ndim=mdata.ndim; mdata=mdata.transpose([2,1,0,*arange(3,ndim)]); zcor=vd.compute_zcor(gd.dp)[sindp].T; nvrt=vd.nvrt
          pie,pip,pacor=gd0.compute_acor(bxy); zs=(vd0.compute_zcor(gd0.dp)[pip]*pacor[...,None]).sum(axis=1).T 
          mdata=interp_vertical(mdata,zs,zcor).transpose([2,1,0,*arange(3,ndim)])
       btime=arange(p.StartT,p.EndT+dt,dt); bdata=interpolate.interp1d(mtime,mdata,fill_value="extrapolate",axis=0)(btime); nt=len(btime)
       dt0=time.time()-t0; print(fname,dt0)

       #create NETCDF
       exec('bdata=bdata[...{}]'.format(',None'*(4-bdata.ndim))) #expand dimension
       if fname.endswith('th.nc'): dimname=['time','nOpenBndNodes', 'nLevels', 'nComponents', 'one']; dims=[*bdata.shape,1]; nvar='time_series' 
       if fname.endswith('nu.nc'): dimname=['time','node','nLevels','one']; dims=[*bdata.shape]; nvar='tracer_concentration'
       nd=zdata(); nd.dimname=dimname; nd.dims=dims; ndict=nd.__dict__
       z=zdata(); z.dimname=('one',);    z.val=array(dt*86400);         nd.time_step=z
       z=zdata(); z.dimname=('time',);   z.val=(btime-btime[0])*86400;  nd.time=z
       z=zdata(); z.dimname=dimname[:4]; z.val=bdata.astype('float32'); ndict[nvar]=z
       if fname.endswith('nu.nc'): z=zdata(); z.dimname=('node',); z.val=sindp+1; nd.map_to_global_node=z
       WriteNC(fname,nd,zlib=True)

#----------------------------------------------------------------------
#source.nc
#Note: same as sbp.run('./write_source.py'); create source.nc
#----------------------------------------------------------------------
fname='source.nc'
if p.flag[fname]==1 and myrank==0:
    print('writing '+fname)

    #--------------------------------------------------------------------------
    #input for NPS and PS loading
    #--------------------------------------------------------------------------
    svars=['temp','salt','clay1','clay2','silt','sand','PB1','PB2','PB3','RPOC','LPOC','DOC',
           'RPON','LPON','DON','NH4','NO23','RPOP','LPOP','DOP','PO4','COD','DO','SRPOC','SRPON','SRPOP','PIP']

    #read data
    C=loadz(p.source); nps,ps,sho=C.nps,C.ps,C.sho
    gd=loadz(p.grd).hgrid; gd.x,gd.y=gd.lon,gd.lat; gd.compute_ctr()

    #--------------------------------------------------
    #remove high nutrient concentration in NPS and PS
    for S in [nps,ps]:
        fpf=S.flow==0; fpc=(S.TN!=0)|(S.TP!=0); S.flow[fpf]=1e-12; S.flow[fpf*fpc]=1e-3 #add some base flow

        #smooth the flow and loading
        for vn, vm in zip(['TN','TP','ORGN','ORGP','NH4','NO3','PO4','PIP'], [100,50,100,100,100,100,50,50]):
            if S==ps and vn=='PIP': vn='BOD'
            for m in arange(len(S.sname)):
                flow=S.flow[m].copy(); mass=S.__dict__[vn][m].copy(); mts=nonzero(mass/(86.4*flow)>vm)[0]
                if len(mts)==0: continue

                #fix high concentration at each point
                fmin=0.1
                for n in mts:
                    fm=mass[n]/(86.4*vm); flow[n]=max(flow[n],min(fm,fmin)) #total flow needed, and add extra flow <0.1
                    if fm<fmin: continue #if the needed flow is <fmin

                    #if the needed flow is > fmin. Then distribute the total loading
                    for k in arange(1,15): #search flow around anomly
                        t1=(n-k); t2=(n+k+1)
                        if sum(flow[t1:t2])>fm: break
                    ftotal=sum(flow[t1:t2])  #total flow available in the time range

                    #distribute loading
                    flow[t1:t2]=flow[t1:t2]+max((fm-ftotal),0)/(t2-t1) #add extra flow if ftotal<fm
                    dm=mass[n]; mass[n]=0; mass[t1:t2]=mass[t1:t2]+dm*flow[t1:t2]/flow[t1:t2].sum() #distribute the loading proportionally to flow

                S.flow[m]=flow; S.__dict__[vn][m]=mass #update flow and mass
    #--------------------------------------------------
    #filter source based on CBSEG information
    #CBSEG=['WSTMH','SOUMH','SEVMH','SASOH','RHDMH','PATMH','NORTF','MIDOH','MAGMH','LCHMH','GUNOH','ELKOH','EASMH','BACOH','BOHOH', #upper bay
    #       'CHSTF','CHSOH','CHSMH','CHOTF','CHOOH','CHOMH2','CB4MH','CB3MH','CB2OH','CB1TF','C&DOH_MD','C&DOH_DE','BSHOH','CHOMH1']
    CBSEG=['RPPTF','RPPMH','RPPOH','CB6PH','CB5MH_VA'] #Rappahannock
    nRiverSeg =[i for i,k in zip(C.RiverSeg,C.CBSEG_92)  if (k not in CBSEG)]

    #pre-proc
    mtime=arange(p.StartT,p.EndT); nt,ntr,nsource=len(mtime),len(svars),len(C.sinde)
    fpn=sort(nonzero((nps.time>=mtime[0])*(nps.time<=mtime[-1]))[0])  #time flag for NPS
    fpp=sort(nonzero((ps.time>=mtime[0])*(ps.time<=mtime[-1]))[0])    #time flag for PS
    fps=sort(nonzero((sho.time>=mtime[0])*(sho.time<=mtime[-1]))[0])  #time flag for SHO
    sinde=near_pts(C.xy,c_[gd.xctr,gd.yctr]) #grid elements

    flow=zeros([nt,nsource]); trs=zeros([nt,ntr,nsource])
    #--------------------------------------------------------------------------
    #NPS input: partitioning coefficients for each river basin
    #--------------------------------------------------------------------------
    c2chl=array([0.075,0.06,0.06]); n2c=array([0.167,0.167,0.167]); p2c=array([0.0125,0.022,0.0125]) #convert ratios
    rivers= ['general','susquehanna','potomac','choptank','james','york','rappahannock','patuxent','chester']; nr=len(rivers)
    vsstc = [2.5,       2.5,          2.5,      2.5,       2,5,    2.5,   2.5,           2.5,       2.5]
    c2n   = [8,         8,            8,        8,         8,      8,     8,             6.0,       6.0]
    fDOC  = [0.7,       0.6,          0.74,     0.67,      0.79,   0.79,  0.64,          0.74,      0.54]
    fLPOC = [0.15,      0.15,         0.15,     0.15,      0.15,   0.15,  0.15,          0.15,      0.15]
    fRPOC = [0.35,      0.35,         0.35,     0.35,      0.35,   0.35,  0.35,          0.35,      0.35]
    fDON  = [0.7,       0.6,          0.74,     0.67,      0.79,   0.79,  0.64,          0.74,      0.54]
    fLPON = [0.15,      0.15,         0.15,     0.15,      0.15,   0.15,  0.15,          0.15,      0.15]
    fRPON = [0.45,      0.45,         0.45,     0.45,      0.45,   0.45,  0.45,          0.45,      0.45]
    fDOP  = [0.35,      0.35,         0.35,     0.355,     0.39,   0.484, 0.228,         0.308,     0.261]
    fLPOP = [0.30,      0.3,          0.30,     0.30,      0.30,   0.3,   0.30,          0.30,      0.30]
    fRPOP = [0.4,       0.4,          0.4,      0.4,       0.4,    0.4,   0.40,          0.4,       0.40]
    fPIP  = [0.6,       0.58,         0.47,     0.7,       0.6,    0.6,   0.6,           0.6,       0.6]
    fPB1  = [0.1,       0.1,          0.1,      0.1,       0.1,    0.1,   0.1,           0.1,       0.1]
    fPB2  = [0.9,       0.9,          0.1,      0.9,       0.9,    0.9,   0.9,           0.9,       0.1]
    fPB3  = [0.01,      0.01,         0.9,      0.01,      0.01,   0.01,  0.01,          0.01,      0.9]

    #collect sources from NPS
    nps.rivers=array([C.rivers[list(C.RiverSeg).index(i)] for i in nps.sname])
    sflow=nps.flow[nps.rivers=='susquehanna'].sum(axis=0) #susquehanna river flow

    for k,[rats,ids] in enumerate(zip(C.nrat,C.nid)):   #each grid elements
        if len(ids)==0: continue                        #no NPS into this elements
        for n,[rat,id] in enumerate(zip(rats,ids)):     #each nps
            if nps.sname[id] in nRiverSeg: continue     #skip if not inside CBSEG
            river=C.rivers[id]; rid=rivers.index(river) #river basin

            #add flow,temp,salt
            flowi=nps.flow[id,fpn]*rat
            temp=nps.temp[id,fpn]*flowi*86.4;  DOX=nps.DO[id,fpn]*flowi*86.4
            COD=zeros(nt); salt=zeros(nt);
            NH4=nps.NH4[id,fpn]*rat; NO3=nps.NO3[id,fpn]*rat; PO4=nps.PO4[id,fpn]*rat
            NH4[NH4<0]=0; NO3[NO3<0]=0; PO4[PO4<0]=0; temp[temp<0]=0; DOX[DOX<0]=0

            #compute algal carbon (mg[C]/L * flow)
            PB1=c2chl[0]*(fPB1[rid]*nps.chla[id,fpn]*flowi*86.4); PB1[PB1<0]=0
            PB2=c2chl[1]*(fPB2[rid]*nps.chla[id,fpn]*flowi*86.4); PB2[PB2<0]=0
            PB3=c2chl[2]*(fPB3[rid]*nps.chla[id,fpn]*flowi*86.4); PB3[PB3<0]=0

            #compute nitrogen
            ORGN=nps.ORGN[id,fpn]*rat-n2c[0]*PB1-n2c[1]*PB2-n2c[2]*PB3; ORGN[ORGN<0]=0
            DON=fDON[rid]*ORGN; PON=ORGN-DON
            RPON=fRPON[rid]*PON;  LPON=fLPON[rid]*PON;  SRPON=PON-RPON-LPON

            #compute carbon
            ORGC=c2n[rid]*(DON+PON); DOC=fDOC[rid]*ORGC; POC=ORGC-DOC
            RPOC=fRPOC[rid]*POC; LPOC=fLPOC[rid]*POC; SRPOC=POC-RPOC-LPOC

            #compute phosphorus
            ORGP=nps.ORGP[id,fpn]*rat+nps.PIP[id,fpn]*rat-p2c[0]*PB1-p2c[1]*PB2-p2c[2]*PB3; ORGP[ORGP<0]=0
            DOP=fDOP[rid]*ORGP; POP=ORGP-DOP; PIP=fPIP[rid]*POP
            RPOP=(1-fPIP[rid])*fRPOP[rid]*POP; LPOP=(1-fPIP[rid])*fLPOP[rid]*POP; SRPOP=(1-fPIP[rid])*POP-RPOP-LPOP

            #special treatment for susquehanna flow
            if river=='susquehanna':
                fpr=sflow>6500

                sfRPON=fRPON[rid]-0.001638*(sflow[fpr]-6500)/100; sfRPON[sfRPON<=0]=0
                sfLPON=fLPON[rid]-0.000749*(sflow[fpr]-6500)/100; sfLPON[sfLPON<=0]=0
                RPON[fpr]=sfRPON*PON[fpr]; LPON[fpr]=sfLPON*PON[fpr]; SRPON[fpr]=(PON-RPON-LPON)[fpr]

                sfRPOC=fRPOC[rid]-0.001330*(sflow[fpr]-6500)/100; sfRPOC[sfRPOC<=0]=0
                sfLPOC=fLPOC[rid]-0.000764*(sflow[fpr]-6500)/100; sfLPOC[sfLPOC<=0]=0
                RPOC[fpr]=sfRPOC*POC[fpr]; LPOC[fpr]=sfLPOC*POC[fpr]; SRPOC[fpr]=(POC-RPOC-LPOC)[fpr]

                sfRPOP=fRPOP[rid]-0.000949*(sflow[fpr]-6500)/100; sfRPOP[sfRPOP<=0]=0
                sfLPOP=fLPOP[rid]-0.001091*(sflow[fpr]-6500)/100; sfLPOP[sfLPOP<=0]=0
                RPOP[fpr]=sfRPOP*POP[fpr]; LPOP[fpr]=sfLPOP*POP[fpr]; SRPOP[fpr]=((1-fPIP[rid])*POP-RPOP-LPOP)[fpr]

            #compute SED
            sand=nps.sand[id,fpn]*rat; silt=nps.silt[id,fpn]*rat; clay1=0.0*nps.clay[id,fpn]*rat; clay2=1.0*nps.clay[id,fpn]*rat

            #add nps load
            flow[:,k]=flow[:,k]+flowi
            trs[:,:,k]=trs[:,:,k]+c_[temp,salt,clay1,clay2,silt,sand,PB1,PB2,PB3,RPOC,LPOC,DOC,
                                     RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4,COD,DOX,SRPOC,SRPON,SRPOP,PIP]

    #-------------------------------------------------------------------------
    # PS input: for adding PS load
    #-------------------------------------------------------------------------
    c2b=10.0; c2n=18.9    #conversion from BOD/Nitrogen to Carbon
    fRPOC=0.04; fLPOC=0.15; fDOC=0.8
    fRPON=0.28; fLPON=0.15; fDON=0.5
    fRPOP=0.42; fLPOP=0.07; fDOP=0.4

    for k,[rats,ids] in enumerate(zip(C.prat,C.pid)):   #each grid elements
        if len(ids)==0: continue                        #no PS into this elements
        for n,[rat,id] in enumerate(zip(rats,ids)):     #each nps
            if C.ps.sname[id][-13:] in nRiverSeg: continue #skip if not inside CBSEG

            #add flow,temp,salt
            flowi=ps.flow[id,fpp]*rat; zn=zeros(nt)
            temp=zn; salt=zn; PB1=zn; PB2=zn; PB3=zn; PIP=zn
            COD=zn; DOX=ps.DO[id,fpp]*rat; TSS=ps.TSS[id,fpp]*rat
            clay1,clay2,silt,sand=0.05*TSS,0.05*TSS,0.45*TSS,0.45*TSS

            #compute carbon, nitrogen, phosphorus
            NH4=ps.NH4[id,fpp]*rat; NO3=ps.NO3[id,fpp]*rat;   PO4=ps.PO4[id,fpp]*rat
            BOD=ps.BOD[id,fpp]*rat; ORGN=ps.ORGN[id,fpp]*rat; ORGP=ps.ORGP[id,fpp]*rat
            NH4[NH4<0]=0; NO3[NO3<0]=0; PO4[PO4<0]=0; BOD[BOD<0]=0; ORGN[ORGN<0]=0; ORGP[ORGP<0]=0

            fp=(BOD==0)*(ORGN>0); ORGC=c2b*BOD;  ORGC[fp]=c2n*ORGN[fp]
            RPOC=fRPOC*ORGC; LPOC=fLPOC*ORGC; DOC=fDOC*ORGC; SRPOC=ORGC-RPOC-LPOC-DOC
            RPON=fRPON*ORGN; LPON=fLPON*ORGN; DON=fDON*ORGN; SRPON=ORGN-RPON-LPON-DON
            RPOP=fRPOP*ORGP; LPOP=fLPOP*ORGP; DOP=fDOP*ORGP; SRPOP=ORGP-RPOP-LPOP-DOP

            #add ps load
            flow[:,k]=flow[:,k]+flowi
            trs[:,:,k]=trs[:,:,k]+c_[temp,salt,clay1,clay2,silt,sand,PB1,PB2,PB3,RPOC,LPOC,DOC,
                                     RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4,COD,DOX,SRPOC,SRPON,SRPOP,PIP]

    #-------------------------------------------------------------------------
    # adding SHO load
    #-------------------------------------------------------------------------
    fRPOC=0.0;  fLPOC=0.0;  fDOC=0.0  #partition of TC
    fRPON=0.0;  fLPON=0.0;  fDON=0.0  #partition of ORGN (=TN)
    fRPOP=0.0;  fLPOP=0.0;  fDOP=0.0  #partition of ORGP
    for k,[rats,ids] in enumerate(zip(C.srat,C.sid)):   #each grid elements
        if len(ids)==0: continue                        #no SHO into this elements
        for n,[rat,id] in enumerate(zip(rats,ids)):     #each nps
            if C.sho.sname[id][-13:] in nRiverSeg: continue #skip if not inside CBSEG

            #add temp,salt, phyto, dissolved nutrients, DO
            zn=zeros(nt); temp=zn; salt=zn; PB1=zn; PB2=zn; PB3=zn; NO3=zn; NH4=zn; PO4=zn; COD=zn; DOX=zn

            #for sediment
            clay=sho.clay[id,fps]*rat; clay1,clay2=0.5*clay,0.5*clay; silt=sho.silt[id,fps]*rat; sand=sho.sand[id,fps]*rat

            #carbon
            ORGC=sho.TC[id,fps]*rat;   RPOC=fRPOC*ORGC; LPOC=fLPOC*ORGC; DOC=fDOC*ORGC; SRPOC=ORGC-RPOC-LPOC-DOC
            ORGN=sho.ORGN[id,fps]*rat; RPON=fRPON*ORGN; LPON=fLPON*ORGN; DON=fDON*ORGN; SRPON=ORGN-RPON-LPON-DON
            ORGP=sho.ORGP[id,fps]*rat; RPOP=fRPOP*ORGP; LPOP=fLPOP*ORGP; DOP=fDOP*ORGP; SRPOP=ORGP-RPOP-LPOP-DOP; PIP=sho.PIP[id,fps]*rat

            #add SHO load
            trs[:,:,k]=trs[:,:,k]+c_[temp,salt,clay1,clay2,silt,sand,PB1,PB2,PB3,RPOC,LPOC,DOC,
                                     RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4,COD,DOX,SRPOC,SRPON,SRPOP,PIP]

    #-------------------------------------------------------------------------
    #convert mass load to concentration
    #-------------------------------------------------------------------------
    #combine loading at each element
    csinde=unique(sinde); nsource=len(csinde); cflow=zeros([nt,nsource]).astype('float32'); ctrs=zeros([nt,ntr,nsource]).astype('float32')
    for k,ie in enumerate(csinde): fpe=sinde==ie; cflow[:,k]=flow[:,fpe].sum(axis=1); ctrs[...,k]=trs[...,fpe].sum(axis=2)
    sinde,flow,trs=csinde,cflow,ctrs; csinde,cflow,ctrs=0,0,0

    #compute conc.
    flow[flow<1.e-8]=1.e-8; trs=trs/(86.4*flow[:,None,:])
    tri=trs[:,-5]; tri[tri>20]=20; trs[:,-5]=tri  #limit DO value
    trs[:,2:6]=trs[:,2:6]/1e3 #convert sediment concentration to g/L

    #check modules
    trs=trs[:,:2]; ntr=2
    #if p.flag['SED']==0: ntr=ntr-4; tid=r_[arange(2),arange(6,27)]; trs=trs[:,tid,:] #remove SED setup
    #if p.flag['ICM']==0: ntr=ntr-21; trs=trs[:,:-21,:] #remove ICM setup
    #if p.flag['ICM'] in [2,20]: ntr=ntr-4; trs=trs[:,:-4,:] #remove ICM's CBPO setup

    #write source.nc
    nd=zdata()
    nd.dimname=['nsources','nsinks','ntracers','time_msource','time_vsource','time_vsink','one']
    nd.dims=[nsource,0,ntr,nt,nt,0,1]

    vi=zdata(); vi.dimname=('nsources',); vi.val=sinde+1; nd.source_elem=vi
    vi=zdata(); vi.dimname=('time_vsource','nsources'); vi.val=flow.astype('float32'); nd.vsource=vi
    vi=zdata(); vi.dimname=('time_msource','ntracers','nsources'); vi.val=trs.astype('float32'); nd.msource=vi
    vi=zdata(); vi.dimname=('one',); vi.val=array(86400.0); nd.time_step_msource=vi
    vi=zdata(); vi.dimname=('one',); vi.val=array(86400.0); nd.time_step_vsource=vi
    vi=zdata(); vi.dimname=('one',); vi.val=array(1e10);    nd.time_step_vsink=vi
    WriteNC('source.nc',nd)

#----------------------------------------------------------------------
#sflux
#----------------------------------------------------------------------
fname='sflux'
if p.flag[fname]==1 and myrank==0:
   print('writing '+fname)

   #parameters
   tdir='./sflux'      #target dir
   itag=1              #itag=[1 or 2],for sflux_air_itag.0691.nc
   tm=[datenum(1990,1,1),datenum(2022,10,31)]

   #make links
   bsdir=os.path.abspath(os.path.curdir); tdir=os.path.abspath(tdir)
   if fexist(tdir): os.rmdir(tdir)
   os.mkdir(tdir); os.chdir(tdir)
   mtime=arange(max(tm[0],p.StartT-2),min(tm[1],p.EndT+2)); svars=['air','rad','prc']
   for irec,ti in enumerate(mtime):
       #link each file
       year=num2date(ti).year; month=num2date(ti).month; day=num2date(ti).day
       for m,svar in enumerate(svars):
           fname='{}/sflux_{}.{:04d}_{:02d}_{:02d}.nc'.format(p.sflux,svar,year,month,day)
           os.symlink(os.path.relpath(fname),'sflux_{}_{}.{:04d}.nc'.format(svar,itag,irec+1))
           if m==1 and month==1 and day==1: print('    sflux: {:04d}-{:02d}-{:02d}'.format(year,month,day))
   #write sflux_inputs.txt
   fid=open('{}/sflux_inputs.txt'.format(tdir),'w+'); fid.write('&sflux_inputs\n   \n/'); fid.close()
   os.chdir(bsdir)
