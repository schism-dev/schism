#!/usr/bin/env python3
#written by Zhengui Wang on 9/20/2022
#Note: This script is still under development/testing, and is for VIMS SCHISM use only.
#      When the development is done, we will make it public at approperiate time.
from pylib import *
import subprocess as sbp
from shutil import copyfile
p=zdata(); p.flag={}

#----------------------------------------------------------------------
#Inputs:
#  base:  reference run for current setup
#  flag=0:
#       if base is not None:  create link to the input file in base run
#       if base is None: skip this file
#  flag=1: re-generate input file
#----------------------------------------------------------------------
p.StartT=datenum(1991,1,1);  p.EndT=datenum(1995,12,31) #simulation time

p.base= 'None' # reference run
p.grid_dir='/sciclone/data10/wangzg/CBP/grid/v9b' #directory of hgrid & vgrid; p.grid_dir=p.base if it is None

#models
p.flag['ICM']= 0  #ICM   model (1: 21 variables; 10: 21-variable offline mode; 2: 17 variables; 20: 17-variable offline mode)
p.flag['SED']= 0  #SED3D model
p.flag['WWM']= 0  #Wave  model

#sub-modules
p.flag['VEG'] = 0 #Hydro Vegetation module
p.flag['CLAM']= 0 #ICM: oyster/clam model
p.flag['SAV'] = 0 #ICM: SAV model
p.flag['WET'] = 0 #ICM: tidal-wetland model

#model inputs
p.flag['elev2D.th.nc']         = 0 #hydro
p.flag['TEM_3D.th.nc']         = 0 #hydro
p.flag['TEM_nu.nc']            = 0 #hydro
p.flag['SAL_3D.th.nc']         = 0 #hydro
p.flag['SAL_nu.nc']            = 0 #hydro
p.flag['sflux']                = 0 #hydro
p.flag['albedo.gr3']           = 0 #hydro
p.flag['watertype.gr3']        = 0 #hydro (1: max=7; 2: max=8)
p.flag['diffmin.gr3']          = 0 #hydro
p.flag['diffmax.gr3']          = 0 #hydro
p.flag['shapiro.gr3']          = 0 #hydro
p.flag['tvd.prop']             = 0 #hydro
p.flag['windrot_geo2proj.gr3'] = 0 #hydro
p.flag['veg_*.gr3']            = 0 #hydro: veg module (SAV drag effect)
p.flag['rough.gr3']            = 0 #hydro,SED
p.flag['bctides.in']           = 0 #hydro,SED,ICM
p.flag['hotstart.nc']          = 0 #hydro,SED,ICM
p.flag['source.nc']            = 0 #hydro,SED,ICM
p.flag['ICM_3D.th.nc']         = 0 #ICM
p.flag['ICM_nu.nc']            = 0 #ICM
p.flag['ICM_sflux.th.nc']      = 0 #ICM
p.flag['ICM_param.nc']         = 0 #ICM
p.flag['SED_hvar_*.ic']        = 0 #SED
p.flag['bed_frac_*.ic']        = 0 #SED
p.flag['bedthick.ic']          = 0 #SED
p.flag['hgrid_WWM.gr3']        = 0 #WWM
p.flag['wwmbnd.gr3']           = 0 #WWM

#databases
p.bdir='/sciclone/data10/wangzg/CBP/database/'         #MBM database dir
p.MBM_init  = p.bdir+'MBM_init.npz'                    #MBM database: multiple datasets
p.source    = p.bdir+'load_p7_v4.npz'                  #CBP watershed sources
p.atmdep    = p.bdir+'atm_load.npz'                    #CBP atmospheric deposition
p.hycom     = p.bdir+'hycom.nc'                        #HYCOM database
p.sflux     = p.bdir+'sflux_narr_subdomain'            #sflux database
p.WW3       = p.bdir+'WW3'                             #WW3 wave forcing
p.hydro_out = p.bdir+'hydro/RUN11/outputs'             #hydro_out for offline ICM model
p.region    = p.bdir+'region/'                         #region files
p.outdir    = '/sciclone/scr-lst/{}/CBP'.format(os.environ['USER']) #parental direcotry of outputs
p.dt_offline= 1800 #sec: time step for offline ICM mode

#----------------------------------------------------------------------
#model configuration (***no change beyond here***)
#----------------------------------------------------------------------
#pre-proc: read grid, and define alias functions
if p.base=='None': p.base=None
if p.grid_dir=='None': p.grid_dir=None
M=read(p.MBM_init,1) #MBM database
M.f2=M.gd.interp #2D interpolator
M.f3=lambda x: interp_schism_3d(M.gd,M.vd,x.T,(gd,vd))  #3D interpolator
def add_var(Z,svar,dim,data): v=zdata(); v.dimname=dim; v.val=data; Z.attr(svar,v)  #define fun to add variables
relpath=lambda x: os.path.relpath(os.path.realpath(x))
def pattr(x,a): p.flag[x]=a

#check reference run, and reset grid_dir
if (p.base is not None) and (not fexist(p.base)): sys.exit('reference run not exist: {}'.format(p.base))
if p.grid_dir is None: p.grid_dir=p.base
ficm=['ICM_3D.th.nc','ICM_nu.nc','ICM_sflux.th.nc','ICM_param.nc','CLAM','SAV','WET']
fsed=['SED_hvar_*.ic','bed_frac_*.ic','bedthick.ic','SED_nudge.gr3']; fwwm=['hgrid_WWM.gr3','wwmbnd.gr3']; p.flag['SED_nudge.gr3']=0
[pattr(i,0) for i in ficm if p.flag['ICM']==0]; [pattr(i,0) for i in fsed if p.flag['SED']==0]; [pattr(i,0) for i in fwwm if p.flag['WWM']==0]

#link or remove files
for fname in p.flag:
    if ((fname in ficm) and p.flag['ICM']==0) or ((fname in fsed) and p.flag['SED']==0) or ((fname in fwwm) and p.flag['WWM']==0): continue
    sn='{}/{}'.format(p.base,fname); sns=glob(sn) if ('*' in fname) else [sn,]; fns=[os.path.basename(i) for i in sns]
    if fname.endswith('_nu.nc'): sns.append(sn[:-6]+'_nudge.gr3'); fns.append(fname[:3]+'_nudge.gr3') #nudge.gr3
    [os.symlink(relpath(s),f) for s, f in zip(sns,fns) if (p.base is not None) and p.flag[fname]==0 and fexist(s) and not fexist(f)]
    if p.flag[fname]==1 and fexist(fname): os.remove(fname) if os.path.isfile(fname) else os.system('rm -rf {}'.format(fname)) #remove file

#add hydro output for ICM offline mode; add istation.in
if p.flag['ICM'] in [10,20]: fn='hydro_out'; [os.remove(fn) if fexist(fn) else '']; os.symlink(p.hydro_out,fn)
if p.flag['ICM']!=0: schism_bpfile(*M.station_xyz.T,M.station_name).save('istation.in')

#create outputs directory
if not fexist('outputs'):
   outdir='{}/{}/outputs'.format(p.outdir,os.path.basename(os.path.abspath('.'))); os.symlink(outdir,'outputs')
   if not fexist(outdir): os.system('mkdir -p {}'.format(outdir))
   if not fexist('outputs'): sys.exit('check p.outdir path: {}'.format(outdir))

#hgrid.gr3, hgrid.ll, vgrid.in, grid.npz
p.grd='grid.npz'
for fname in ['hgrid.gr3','hgrid.ll','vgrid.in','grid.npz']:
    print('writing '+fname)
    sname=relpath('{}/{}'.format(p.grid_dir,fname))
    if not fexist(sname): sys.exit('{}/{} not exist'.format(p.grid_dir,fname))
    if os.path.islink(fname) or (fexist(fname) and sname!=fname): os.remove(fname)
    if sname!=fname: os.symlink(sname,fname)
gd,vd=grd(p.grd,fmt=2) #read grid

#----------------------------------------------------------------------
#parameter files
#----------------------------------------------------------------------
#copy parameter files
fnames=['param.nml','sediment.nml','wwminput.nml','icm.nml']
modules=['hydro','SED','WWM','ICM']; p.flag['hydro']=1
for fname,module in zip(fnames,modules):
    if fexist(fname) or p.flag[module]==0: continue
    sname='{}/{}'.format(p.base,fname); pname=p.bdir+'param/'+fname
    if fexist(sname):
       copyfile(sname,fname); print('copy {} from {}'.format(fname,p.base))
    else:
       copyfile(pname,fname); print('writing {}'.format(fname))

#change parameter values in [param.nml, wwminput.nml,icm.nml]
fname='param.nml'; sname='{}/{}'.format(p.base,fname); pm={}; pm['ihot']=1; pm['iveg']=p.flag['VEG'] #reset ihot=1
if not fexist(sname):
   if p.flag['SED']==1: pm['inu_tr(5)']=1; pm['nspool']=12; pm['iof_sed(27)']=1
   if p.flag['WWM']==1: pm['ics']=2; pm['nspool']=12; pm['icou_elfe_wwm']=1; pm['iwbl']=2
   if p.flag['ICM_nu.nc']==1: pm['inu_tr(7)']=2; pm['vnf1']=0.01
   if p.flag['SED']==1 and p.flag['WWM']==1: #for offline transport mode
      for i in [13,21,27,29,30]:  pm['iof_hydro({})'.format(i)]=1
   s=num2date(p.StartT); pm['start_year']=s.year; pm['start_month']=s.month; pm['start_day']=s.day; pm['rnday']=int(p.EndT-p.StartT)
if p.flag['ICM'] in [10,20]: #offline ICM mode
   pm['nspool']=int(6*3600/p.dt_offline); pm['ihfskip']=int(720*3600/p.dt_offline)
   pm['nhot_write']=pm['ihfskip']; pm['dt']=p.dt_offline; pm['wtiminc']=pm['dt']
   pm['itransport_only']=2; pm['nadv']=1; pm['itr_met']=3; pm['ielm_transport']=1; pm['max_subcyc']=30
for i in pm: chparam(fname,i,pm[i]) #update param.nml

fname='wwminput.nml'; sname='{}/{}'.format(p.base,fname); pm={}
if p.flag['WWM']==1 and (not fexist(sname)):
   pm['BEGTC_OUT']=num2date(p.StartT).strftime('%Y%m%d.000000'); pm['BEGTC']="'{}'".format(pm['BEGTC_OUT'])
   pm['ENDTC_OUT']=num2date(p.EndT).strftime('%Y%m%d.000000');   pm['ENDTC']="'{}'".format(pm['ENDTC_OUT'])
   for i in pm: chparam(fname,i,pm[i]) #update wwminput.nml

fname='icm.nml'; sname='{}/{}'.format(p.base,fname); pm={}; pm0=read(fname,1)
pm['iout_icm']=2; pm['iClam']=0; pm['isav_icm']=0; pm['imarsh_icm']=0
if p.flag['ICM']!=0: 
   if p.flag['ICM'] in [1,10] and (not fexist(sname)): pm['iSRM']=1
   if p.flag['SED']==1: pm['iKe']=1
   if p.flag['ICM_sflux.th.nc']==1: pm['isflux']=1
   if (p.flag['ICM'] in [1,2]): pm['nspool_icm']=24
   if p.flag['ICM'] in [10,20]: pm['nsub']=int(p.dt_offline/150); pm['iKe']=1; pm['nspool_icm']=1
   if p.flag['CLAM']==1: pm['iClam']=1; pm['cFc']='-999   '*5; pm['cMTB']='-999   '*5   
   if p.flag['SAV']==1: pm['isav_icm']=2; pm['sFc']=-999
   if p.flag['WET']==1: pm['imarsh_icm']=2; pm['vAw']=-999
   if not fexist(sname): #change WSP, WSPn, KC0/KN0/KP0
      for i in [3,4,6,7,11,12]: pm0['WSP'][i]=-999; pm0['WSPn'][i]=-999
      pm0['KC0'][2]=-999; pm0['KN0'][2]=-999; pm0['KP0'][2]=-999;
      for i in ['WSP','WSPn','KC0','KN0','KP0']: pm[i]='  '.join([str(k) for k in pm0[i]])
   for i in pm: chparam(fname,i,pm[i]) #update icm.nml

#----------------------------------------------------------------------
#bctides.in
#----------------------------------------------------------------------
fname='bctides.in'
if p.flag[fname]==1:
   print('writing '+fname)
   fid=open('./bctides.in','w+')
   fid.write('{} GMT\n'.format(num2date(p.StartT).strftime('%m/%d/%Y %H:%M:%S')))
   lstr='0 40. ntip\n0 nbfr\n1 nope\n{:d} 4 0 4 4'.format(gd.nobn[0])
   bstr='\n1.0 !TEM nudge\n1.0 !SAL nudge'
   if p.flag['SED']==1: lstr=lstr+' 3'; bstr=bstr+'\n0.1 !SED nudge'
   if p.flag['ICM']!=0: lstr=lstr+' 4'; bstr=bstr+'\n1.0 !ICM nudge'
   fid.write(lstr+bstr); fid.close()

#----------------------------------------------------------------------
#gr3,prop files
#----------------------------------------------------------------------
fnames=('albedo.gr3','windrot_geo2proj.gr3','bedthick.ic', 'diffmax.gr3',)
fvalues=(0.1,          0.0,    5.0,   0.1,)
for fname,fvalue in zip(fnames,fvalues):
    if p.flag[fname]==1:
       print('writing '+fname)
       gd.save(fname,value=fvalue)

fname='rough.gr3'
if p.flag[fname]!=0:
   print('writing '+fname)
   gd.save(fname,value=M.f2(gd.xy,value=M.rough),outfmt='{:16.8f}')

fname='watertype.gr3'
if p.flag[fname]!=0:
   print('writing '+fname)
   vi=M.f2(gd.xy,value=M.watertype.clip(0,7) if p.flag[fname]==1 else M.watertype)
   gd.save(fname,value=vi.astype('int'),outfmt='{:d}')

fname='diffmin.gr3'
if p.flag[fname]==1:
   print('writing '+fname)
   vi=ones(gd.np)*1e-6; vi[read(p.region+'lowerbay.reg').inside(gd.xy)]=8e-5
   gd.save(fname,value=vi)

fname='shapiro.gr3'
if p.flag[fname]==1:
   print('writing '+fname)
   shapiro_max,threshold_slope=0.5,0.5
   slope=gd.compute_gradient(fmt=2)[2]; vi=shapiro_max*tanh(2*slope/threshold_slope)
   gd.save('shapiro.gr3',value=vi)

fname='veg_*.gr3'
if p.flag[fname]==1:
   print('writing '+fname)
   for i in ['h','D','cd','N']: gd.save('veg_{}.gr3'.format(i),value=M.f2(gd.xy,M.attr('veg_'+i))) 

fname='SED_nudge.gr3'
if p.flag['SED']==1 and not fexist(fname): 
   dist=abs(gd.cxy[:,None]-gd.cxy[gd.iobn[0]][None,:]).min(axis=1); x1,x2=[1e3,3e4]
   vi=1.15740741e-05*((dist-x2)/(x1-x2)).clip(0,1); gd.save(fname,value=vi,outfmt='{:16.8e}') 

fname='tvd.prop'
if p.flag[fname]==1:
   print('writing '+fname)
   vi=ones(gd.ne); vi[read(p.region+'stream_head.reg').inside(gd.exy)]=0
   gd.save('tvd.prop',value=vi,fmt='{:d}')

fname='SED_hvar_*.ic'
if p.flag[fname]==1:
   print('writing '+fname)
   for i in arange(4): gd.save('SED_hvar_{}.ic'.format(i+1),value=0)

fname='bed_frac_*.ic'
if p.flag[fname]==1:
   print('writing '+fname)
   for i in arange(4): gd.save('bed_frac_{}.ic'.format(i+1),value=M.f2(gd.xy,value=M.SED3D_bedfrac[i]))

fname='hgrid_WWM.gr3'
if p.flag[fname]==1:
   print('writing '+fname)
   gdt=read(p.grd,'hgrid'); gdt.x,gdt.y=gdt.lon,gdt.lat; gdt.split_quads_wwm(fname)

fname='wwmbnd.gr3'
if p.flag[fname]==1:
   print('writing '+fname)
   gdt=read(p.grd,'hgrid'); gdt.dp[:]=0; gdt.dp[gd.iobn[0]]=2; gdt.save(fname)

#----------------------------------------------------------------------
#hotstart.nc
#----------------------------------------------------------------------
fname='hotstart.nc'
if p.flag[fname]==1:
   print('writing '+fname)

   #init. condition for salinity and temp
   tid=abs(M.temp_time-(p.StartT-datenum(num2date(p.StartT).year,1,1))).argmin()
   if (p.StartT>=M.salt_time.min())*(p.StartT<=M.salt_time.max()):
      sid=abs(M.salt_time-p.StartT).argmin()
   else:
      mm=array([num2date(i).month for i in M.salt_time]); it=pindex(mm==num2date(p.StartT).month); sid=it[abs(M.salt_time[it]-p.StartT).argmin()]
   trs=c_[M.f3(M.temp_data[tid])[...,None],M.f3(M.salt_data[sid])[...,None]].astype('float32')

   #write hotstart: init. hydro
   H=zdata(); H.file_format='NETCDF4'; nn,ne,ns,nvrt,ntr=gd.np,gd.ne,gd.ns,vd.nvrt,2
   H.dimname=['node','elem','side','nVert','ntracers','one','two','three','four','five']; H.dims=[nn,ne,ns,nvrt,ntr,1,2,3,4,5]
   for i,k in zip(['time','iths','ifile','nsteps_from_cold'],[0.0,0,1,0]): add_var(H,i,('one',),array(k)) #single 
   for i,k,m in zip(['idry_e','idry_s','idry'],['elem','side','node'],[ne,ns,nn]): add_var(H,i,(k,),zeros(m,'int32')) #idry_e,idry_s,idry
   for i in ['eta2','cumsum_eta']: add_var(H,i,('node',),zeros(nn,'float32')) #eta2, and cumsum_eta
   for i in ['we','su2','sv2','q2','xl','dfv','dfh','dfq1','dfq2']: #2d variables
       dm,npt=['elem',ne] if i=='we' else ['side',ns] if (i in ['su2','sv2']) else ['node',nn]; add_var(H,i,(dm,'nVert'),zeros([npt,nvrt],'float32'))

   #SED3D
   if p.flag['SED']==1:
      print('  add SED3D model')
      H.dimname.extend(['sed_class','nbed', 'MBEDP']); H.dims.extend([4,1,3]); ntr=ntr+4; H.dims[H.dimname.index('ntracers')]=ntr #update dims
      for i in ['SED3D_dp','SED3D_rough']: add_var(H,i,('node',),M.f2(gd.xy,value=M.attr(i)).astype('float32')) #dp, rough
      add_var(H,'SED3D_bed',('elem','nbed','MBEDP'), array([M.f2(gd.exy,value=i) for i in M.SED3D_bed]).T[:,None,:].astype('float32')) #bed
      add_var(H,'SED3D_bedfrac',('elem','nbed','sed_class'), array([M.f2(gd.exy,value=i) for i in M.SED3D_bedfrac]).T[:,None,:].astype('float32')) #bedfrac
      trs=resize(trs,[nn,nvrt,ntr]) #SED state variables

   #ICM
   if p.flag['ICM']!=0:
      print('  add ICM model')
      ntr_icm=21 if (p.flag['ICM'] in [1,10]) else 17
      svars=['PB1','PB2','PB3','RPOC','LPOC','DOC','RPON','LPON','DON','NH4','NO3','RPOP','LPOP','DOP','PO4','COD','DOX','SRPOC','SRPON','SRPOP','PIP']
      bvars=['bPOC','bPON','bPOP','btemp','bstc','bSTR','bThp','bTox','bNH4','bNH4s','bNO3','bPO4','bH2S','bCH4']
      ntr=ntr+ntr_icm; H.dims[H.dimname.index('ntracers')]=ntr #update dims 
      for i in bvars[:3]: add_var(H,i,('elem','three'),array([M.f2(gd.exy,k) for k in M.attr(i)]).T) #bPOM
      for i in bvars[3:]: add_var(H,i,('elem',),M.f2(gd.exy,M.attr(i))) #other SFM variables
      trs=c_[trs,array([M.f3(M.attr(i)).astype('float32') for i in svars[:ntr_icm]]).transpose([1,2,0])] #ICM state variables
      if p.flag['CLAM']==1: add_var(H,'clam',('elem','five'), array([M.f2(gd.exy,i) for i in M.clam]).T)   #CLAM   
      if p.flag['SAV']==1: add_var(H,'sav',('elem','four'), array([M.f2(gd.exy,i) for i in M.sav]).T); add_var(H,'EP', ('elem',), M.f2(gd.exy,M.EP))#SAV

   #write hotstart
   for i in ['tr_nd','tr_nd0']: add_var(H,i,('node','nVert','ntracers'),trs) 
   add_var(H,'tr_el',('elem','nVert','ntracers'),gd.interp_node_to_elem(trs).astype('float32'))
   H.save('hotstart.nc',zlib=True); trs=None

#---------------------------------------------------------------------------------------
#elev2D.th.nc: Ocean boundary nodes are fixed 
#---------------------------------------------------------------------------------------
fname='elev2D.th.nc'
if p.flag[fname]==1:
   print('writing '+fname)
  
   #get elev data (do not change boundary) 
   fpt=(M.elev_time>=p.StartT)*(M.elev_time<=(p.EndT+2)); mti=M.elev_time[fpt]
   sindp=near_pts(gd.xy[gd.iobn[0]],M.elev_xy[:])
   elev=M.elev_data[:,fpt][sindp].T

   #write elev2D.th.nc
   C=zdata(); C.dimname=['nOpenBndNodes','nLevels','nComponents','one','time']; C.dims=[gd.nobn[0],1,1,1,len(mti)]
   add_var(C,'time_step',('one',),array(3600.0)); add_var(C,'time',('time',),(mti-mti[0])*86400)
   add_var(C,'time_series',('time','nOpenBndNodes','nLevels','nComponents'), elev[...,None,None])
   C.save(fname,zlib=True)

#----------------------------------------------------------------------
#TEM_3D.th.nc,  TEM_nu.nc, TEM_nudge.gr3
#----------------------------------------------------------------------
for fname in ['TEM_3D.th.nc','TEM_nu.nc']:
    if p.flag[fname]==0: continue
    print('writing '+fname)
    nn,ne,ns,nvrt=gd.np,gd.ne,gd.ns,vd.nvrt; fmt=0 if fname.endswith('th.nc') else 1

    #write TEM_nudge.gr3
    if fname=='TEM_nu.nc':
       dist=abs(gd.cxy[:,None]-gd.cxy[gd.iobn[0]][None,:]).min(axis=1); x1,x2=[8e3,5e4]; 
       vi=1.15740741e-05*((dist-x2)/(x1-x2)).clip(0,1); gd.save('TEM_nudge.gr3',value=vi,outfmt='{:16.8e}') 
   
    #read hycom database
    C=read(p.hycom,1); ct=C.time/24+datenum(2000,1,1); cvar=C.water_temp
    clon=C.lon[:]; clat=C.lat[:]; cdepth=C.depth[:]; cxy=(clon[None,:]+1j*clat[:,None]).ravel()

    #get node xyz, and get interp indices
    bind=gd.iobn[0] if fmt==0 else pindex(read(fname[:4]+'nudge.gr3').z!=0)
    lx,ly=gd.lxy[bind].T; lz=abs(compute_zcor(vd.sigma[bind],gd.dp[bind])).T; nobn=bind.size
    idx=(lx[:,None]>clon[None,:]).sum(axis=1)-1; ratx=(lx-clon[idx])/(clon[idx+1]-clon[idx]); ratx=ratx[None,:]
    idy=(ly[:,None]>clat[None,:]).sum(axis=1)-1; raty=(ly-clat[idy])/(clat[idy+1]-clat[idy]); raty=raty[None,:]
    mti=arange(p.StartT,p.EndT+2.0); dt=90; irec=0 #each time, do interplation for 30 days

    #do interpolation and write netcdf file
    dns=['nOpenBndNodes' if fmt==0 else 'node','nLevels','one','time']; dms=[nobn,nvrt,1,None]
    from netCDF4 import Dataset
    fid=Dataset(fname,'w',format='NETCDF4'); fvar=fid.variables
    for i,k in zip(dns,dms): fid.createDimension(i,k)   #define dimension
    if fmt==0: fid.createDimension('nComponents',1)
    mvar='time_series' if fmt==0 else 'tracer_concentration' 
    fid.createVariable('time','float64',('time',)); fvar['time'][:]=(mti-mti[0])*86400.0
    if fmt==0:
       fid.createVariable('time_step','float64',('one',)); fvar['time_step'][0]=86400.0 
       fid.createVariable(mvar,'float32',('time','nOpenBndNodes','nLevels','nComponents'),zlib=True) 
    else:
       fid.createVariable('map_to_global_node','int',('node',)); fvar['map_to_global_node'][:]=bind+1 
       fid.createVariable(mvar,'float32',('time','node','nLevels','one'),zlib=True)

    #interp
    while(irec<mti.size):
         it1=irec; it2=min(it1+dt,mti.size); irec=irec+dt; its=arange(pindex(ct<mti[it1]).max(),pindex(ct>mti[it2-1]).min()+1)
         vs=cvar[its]; temp=[]; print('  {}: {}'.format(fname,num2date(mti[it2-1]))) #get temp data 
         for i, v in enumerate(vs): #fix nan
             if it1==0 and i==0: #compute interp indices for nan
                ins,ips=[],[]; fps=[abs(m.ravel())>1e3 for m in v]; [[ins.append(pindex(m)),ips.append(pindex(~m))] for m in fps]  
                sindns=ins; sindps=[n if n.size==0 else n[near_pts(cxy[m],cxy[n])] for m,n in zip(ins,ips)] 
             for k,cv in enumerate(v): #fix nan 
                 sindn,sindp=sindns[k],sindps[k]; cv=cv.ravel()
                 if sindp.size!=0: fp0=(abs(cv[sindn])>1e3)*(abs(cv[sindp])<1e3); cv[sindn[fp0]]=cv[sindp[fp0]] #init fix
                 fpn=abs(cv)>1e3
                 if sum(~fpn)==0:
                    v[k]=v[k-1]
                 else: 
                    fn=pindex(fpn); ip=pindex(~fpn); fp=ip[near_pts(cxy[fn],cxy[ip])]; cv[fn]=cv[fp]  #final fix

             #interp in space 
             v11=v[:,idy,idx]; v12=v[:,idy,idx+1]; v21=v[:,idy+1,idx]; v22=v[:,idy+1,idx+1]
             vi=(v11*(1-ratx)+v12*ratx)*(1-raty)+(v21*(1-ratx)+v22*ratx)*raty 
             fv=interpv(vi,cdepth,lz).astype('float32'); temp.append(fv)
         fvar[mvar][it1:it2]=interp(ct[its],array(temp),mti[it1:it2],axis=0).transpose([0,2,1])[...,None]
    fid.close(); C.close(); temp=None

#----------------------------------------------------------------------
#SAL_3D.th.nc, SAL_nu.nc, SAL_nudge.gr3
#----------------------------------------------------------------------
for fname in ['SAL_3D.th.nc','SAL_nu.nc']:
    if p.flag[fname]==0: continue
    print('writing {}'.format(fname))
    nn,ne,ns,nvrt=gd.np,gd.ne,gd.ns,vd.nvrt; fmt=0 if fname.endswith('th.nc') else 1

    #write SAL_nudge.gr3
    if fname=='SAL_nu.nc':
       dist=abs(gd.cxy[:,None]-gd.cxy[gd.iobn[0]][None,:]).min(axis=1); x1,x2=[8e3,1.35e4]; fnudge=1.15740741e-05
       vi=fnudge*((dist-x2)/(x1-x2)).clip(0,1); vi[read(p.region+'CD_canal.reg').inside(gd.lxy)]=fnudge 
       gd.save('SAL_nudge.gr3',value=vi,outfmt='{:16.8e}')

    #get salinity climatology database, and clip salinity data, and shrink datasize 
    ctime=M.csalt_time[:]; clon=M.csalt_lon[:]; clat=M.csalt_lat[:]; cdepth=M.csalt_depth[:]; csalt=M.csalt_data[:]
    fpz=cdepth<=5; fpv=csalt>31.5;  csalt[fpz[None,:,None,None]*fpv]=31.5
    fpz=(cdepth>5)*(cdepth<=10); fpv=csalt>33; csalt[fpz[None,:,None,None]*fpv]=33

    fpx=pindex((clon>=gd.lon.min())*(clon<=gd.lon.max())); fpy=pindex(clat<=gd.lat.max()) 
    ix1=fpx.min()-1; ix2=fpx.max()+2; iy=fpy.max()+2; 
    clon=clon[ix1:ix2]; clat=clat[:iy]; csalt=csalt[:,:,:iy,ix1:ix2]; cxy=(clon[None,:]+1j*clat[:,None]).ravel()

    for i in arange(ctime.size): #remove nan   
        for k in arange(cdepth.size):
            cv=csalt[i,k].ravel(); fpn=isnan(cv); sindn=pindex(fpn); sindp=pindex(~fpn)
            fp=near_pts(cxy[sindn],cxy[sindp]); cv[sindn]=cv[sindp[fp]]; csalt[i,k]=cv.reshape(csalt.shape[2:]) 

    #interpolation: get node xyz, interp indices, then interp horizontal, vertical and time
    bind=gd.iobn[0] if fmt==0 else pindex(read(fname[:4]+'nudge.gr3').z!=0)
    lx,ly=gd.lxy[bind].T; ly=ly.clip(clat.min(),90); lz=abs(compute_zcor(vd.sigma[bind],gd.dp[bind])).T; nobn=bind.size
    idx=(lx[:,None]>clon[None,:]).sum(axis=1)-1; ratx=(lx-clon[idx])/(clon[idx+1]-clon[idx]); ratx=ratx[None,None,:]
    idy=(ly[:,None]>clat[None,:]).sum(axis=1)-1; raty=(ly-clat[idy])/(clat[idy+1]-clat[idy]); raty=raty[None,None,:]
    v11=csalt[...,idy,idx]; v12=csalt[...,idy,idx+1]; v21=csalt[...,idy+1,idx]; v22=csalt[...,idy+1,idx+1]
    v0=(v11*(1-ratx)+v12*ratx)*(1-raty)+(v21*(1-ratx)+v22*ratx)*raty #horizontal
    v1=interpv(v0.transpose([1,2,0]),cdepth,lz).transpose([2,1,0]) #vertical

    ctime=r_[ctime[-1]-365,ctime,ctime[0]+365]; v1=r_[v1[-1:,...],v1,v1[:1,...]]
    mti=arange(p.StartT,p.EndT+2); doy=array([i-datenum(num2date(i).year,1,1) for i in mti])
    salt=interp(ctime,v1,doy,axis=0)

    #for CD_canal
    if fmt==1:
       fpb=read(p.region+'CD_canal.reg').inside(gd.lxy[bind])
       salt[:,fpb]=interp(M.ET21_time,M.ET21_salt,mti)[:,None,None]

    #write netcdf
    dns=['nOpenBndNodes' if fmt==0 else 'node','nLevels','one','time']; dms=[nobn,nvrt,1,mti.size]
    if fmt==0: dns.append('nComponents'); dms.append(1)
    C=zdata(); C.dimname=dns; C.dims=dms
    add_var(C,'time',('time',),(mti-mti[0])*86400) #add time 
    if fmt==0:
       add_var(C,'time_step',('one',),array(86400.0))
       add_var(C,'time_series',('time','nOpenBndNodes','nLevels','nComponents'),salt[...,None])
    else:
       add_var(C,'map_to_global_node',('node',),bind+1)
       add_var(C,'tracer_concentration',('time','node','nLevels','one'),salt[...,None])
    C.save(fname,zlib=True); C,salt=None,None

#----------------------------------------------------------------------
#ICM_3D.th.nc
#----------------------------------------------------------------------
fname='ICM_3D.th.nc'
if p.flag[fname]!=0:
   print('writing '+fname)

   #inputs
   svars=array(['PB1','PB2','PB3','RPOC','LPOC','DOC','RPON','LPON','DON','NH4','NO3', 'RPOP','LPOP','DOP','PO4','COD','DO','SRPOC','SRPON','SRPOP','PIP'])
   svals=array([0.01, 0.01, 0.01,  0.01, 0.01,  0.01, 0.001, 0.001, 0.001, 0,    0,     1e-4,  1e-4, 1e-4,  0,   0,     0,   0,      0,     0,       0  ])
   if p.flag['ICM'] in [2,20]: svars=svars[:-4]; svals=svals[:-4]
   ctime=M.NEFSC_time; cdepth=M.NEFSC_depth  #saved database
   S={}; S['DO']=M.NEFSC_DO; S['NH4']=M.NEFSC_NH4; S['NO3']=M.NEFSC_NO3; S['PO4']=M.NEFSC_PO4

   #interp
   nobn=gd.nobn[0]; bind=gd.iobn[0]; lz=abs(compute_zcor(vd.sigma[bind],gd.dp[bind])).clip(cdepth.min(),cdepth.max()).T
   mti=arange(p.StartT,p.EndT+2); doy=array([i-datenum(num2date(i).year,1,1)+1 for i in mti])
   nt=len(mti); nvrt=vd.nvrt; trs=[]
   for m,[svar, sv] in enumerate(zip(svars,svals)):
       if svar not in S: trs.append((ones([nvrt,nobn,ctime.size])*sv).astype('float32')); continue #constant values 
       cv=interpv(tile(S[svar][...,None],nobn).transpose([1,2,0]),cdepth,lz).astype('float32'); trs.append(cv)
   trs=array(trs).transpose([3,2,1,0]); ctime=r_[ctime[-1]-365,ctime,ctime[0]+365]; trs=r_[trs[-1:,...],trs,trs[:1,...]]
   trs=interp(ctime,trs,doy,axis=0).astype('float32')

   #write ICM_3D.th.nc
   C=zdata(); C.dimname=['nOpenBndNodes', 'nLevels', 'nComponents', 'one', 'time']; C.dims=[nobn,nvrt,len(svars),1,nt]
   add_var(C,'time_step',('one',),array(86400.0));  add_var(C,'time',('time',),(mti-mti[0])*86400)
   add_var(C,'time_series',('time','nOpenBndNodes','nLevels','nComponents'),trs)
   C.save(fname,zlib=True); C,trs=None,None

#----------------------------------------------------------------------
#ICM_nu.nc, ICM_nudge.gr3
#----------------------------------------------------------------------
fname='ICM_nu.nc'
if p.flag[fname]!=0:
   print('writing '+fname)

   #write TEM_nudge.gr3
   fnudge=1.15740741e-05; vi=zeros(gd.np); vi[read(p.region+'ICM_nudge.reg').inside(gd.xy)]=fnudge
   gd.save('ICM_nudge.gr3',value=vi,outfmt='{:16.8e}')

   #inputs
   vnames=['PB1',    'PB2',   'PB3','RPOC','LPOC','DOC','RPON','LPON','DON','NH4', 'NO3',  'RPOP','LPOP','DOP', 'PO4', 'COD','DOX','SRPOC','SRPON','SRPOP','PIP']
   svars =['CHLA',   'CHLA',  'CHLA','POC','POC', 'DOC', 'PON','PON', 'DON','NH4', 'NO3',  'POP',  'POP', 'DOP','PO4', 'DON','DO' ,'POC',   'PON', 'POP',  'POP']
   rats  =[0.5/13.3, 0.5/16.6, 0.0,  0.1,   0.9,   1.0,  0.1,   0.9,   1.0,  1.0,   1.0,   0.1,    0.9,   1.0,   1.0,   0.0,  1.0,  0.0,    0.0,    0.0,   0.0]
   if p.flag['ICM'] in [2,20]: vnames,svars,rats=vnames[:-4],svars[:-4],rats[:-4]

   #interp
   fpt=(M.CB74_time>=p.StartT)*(M.CB74_time<p.EndT+2); mti=M.CB74_time[fpt]; cz=M.CB74_depth
   bind=pindex(read('ICM_nudge.gr3').z!=0); lz=abs(compute_zcor(vd.sigma[bind],gd.dp[bind])).T.clip(cz.min(),cz.max())
   nobn=bind.size; ntr=len(svars); nt=len(mti); nvrt=vd.nvrt
   mys=array([M.attr('CB74_'+svar)[:,fpt]*rat for svar, rat in zip(svars,rats)]).transpose([1,0,2]) #interp in time

   #interp vertically, and write ICM_nu.nc 
   from netCDF4 import Dataset; dt=30; irec=0
   fid=Dataset(fname,'w',format='NETCDF4'); fvar=fid.variables; fdim=fid.dimensions    
   dimname=['node', 'nLevels', 'ntracers', 'time']; dims=[nobn,nvrt,ntr,None]
   for i,k in zip(dimname,dims): fid.createDimension(i,k)   #define dimension
   fid.createVariable('time','float64',('time',)); fid.createVariable('map_to_global_node','int',('node',))
   fid.createVariable('tracer_concentration','float32',('time','node','nLevels','ntracers'),zlib=True)
   while(irec<nt): #put variable values
       fpt=arange(irec,min(irec+dt,nt)); trs=zeros([nvrt,nobn,ntr,len(fpt)],'float32')*nan; irec=irec+dt
       for k, zi in enumerate(lz):
           for i in arange(cz.size-1): trs[k,(zi>=cz[i])*(zi<cz[i+1])]=mys[i,:,fpt].T[None,...]
       for k in arange(nvrt-1)[::-1]: fpn=isnan(trs[k]); trs[k][fpn]=trs[k+1][fpn]  #remove nan
       fvar['tracer_concentration'][fpt]=trs.transpose([3,1,0,2])
   fvar['time'][:]=(mti-mti[0])*86400.0; fvar['map_to_global_node'][:]=bind+1
   fid.close(); trs,mys=None,None

#----------------------------------------------------------------------
#ICM_sflux.th.nc: atmosheric loading
#----------------------------------------------------------------------
fname='ICM_sflux.th.nc'
if p.flag[fname]!=0:
   print('writing '+fname)

   #Input
   svars=['NA','NA','NA','NA','NA','NA','ORGN','ORGN','ORGN','NH4','NO3','ORGP','ORGP','ORGP','PO4','NA','NA', 'NA','ORGN','ORGP','NA']
   srats=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  1.0,   0.0,   0.0,   1.0,  1.0,  1.0,   0.0,   0.0,   1.0,  0.0, 0.0, 0.0,  0.0,   0.0,   0.0]
   if p.flag['ICM'] in [2,20]: svars,rats=svars[:-4],srats[:-4]
   S=read(p.atmdep,1); tm=S.time.max(); ym=num2date(tm).year  #atm depostion database

   #interp
   sindp=read(p.region+'estuary.reg').inside(gd.xy,prj=['epsg:26918','epsg:4326']); sinda=near_pts(gd.xy[sindp],c_[S.x,S.y])
   mti=arange(p.StartT,p.EndT+2); mti=array([datenum(ym,1,1)+i-datenum(num2date(i).year,1,1) if i>tm else i for i in mti])
   nt=len(mti); ntr=len(svars); trs=zeros([ntr,gd.np,nt],'float32')
   for m, [svar,srat] in enumerate(zip(svars,srats)):
       if srat==0: continue
       trs[m,sindp]=(interp(S.time,S.attr(svar),mti,axis=1)*1e-3*srat/S.area[:,None])[sinda]

   #write netcdf
   C=zdata(); C.dimname=['node','time','ntr','one']; C.dims=[gd.np,nt,ntr,1]
   add_var(C,'time_step',('one',),array([86400.0])); 
   add_var(C,'time_series',('time','ntr','node'), trs.transpose([2,0,1]))
   C.save(fname,zlib=True); C,S,trs=None,None,None
 
#----------------------------------------------------------------------
#ICM_param.nc: write ICM model spatially varying parameters
#----------------------------------------------------------------------
fname='ICM_param.nc'
if p.flag[fname]!=0:
   print('writing '+fname)

   #read grid, get get region indices
   gd=loadz(p.grd).hgrid; gd.compute_ctr()
   sindr=read(p.region+'stream_head.reg').inside(gd.xy)  #upper river region
   sindo=read(p.region+'estuary_ext.reg').inside(gd.xy,fmt=2)
   sindoo=read(p.region+'estuary.reg').inside(gd.xy,fmt=2,prj=['epsg:26918','epsg:4326'])
   fpz0=abs(gd.dp)<1.0        #depth<1m
   fpz=(gd.dp>=1)*(gd.dp<=11) #depth in [1,11]

   #settling velocity: WSP and WSPn
   ntr=21 if p.flag['ICM'] in [1,10] else 17; WSP=zeros([ntr,gd.np]); WSPn=zeros([ntr,gd.np])

   #POC
   for i in [3,4]:
      vm=[0.1, 0.5];   WSP[i]=vm[1];  WSP[i,fpz0]=vm[0];   WSP[i,fpz]=vm[0]+(gd.dp[fpz]-1.0)*(vm[1]-vm[0])/10
      vm=[0.02, 0.2];  WSPn[i]=vm[1]; WSPn[i,fpz0]=vm[0];  WSPn[i,fpz]=vm[0]+(gd.dp[fpz]-1.0)*(vm[1]-vm[0])/10

      vi=gd.dp[sindo]/10; vi[vi>0.5]=0.5; vi[vi<0]=0; WSP[i,sindo]=vi #ocean
      vi=gd.dp[sindo]/50; vi[vi>0.2]=0.2; vi[vi<0]=0;
      fp=(gd.dp[sindo]<=20)*(vi>0.02); vi[fp]=0.02; WSPn[i,sindo]=vi #ocean

   #for PON
   for i in [6,7]:
       vm=[0.1, 0.5];   WSP[i]=vm[1];  WSP[i,fpz0]=vm[0];  WSP[i,fpz]=vm[0]+(gd.dp[fpz]-1.0)*(vm[1]-vm[0])/10
       vm=[0.05, 0.35]; WSPn[i]=vm[1]; WSPn[i,fpz0]=vm[0]; WSPn[i,fpz]=vm[0]+(gd.dp[fpz]-1.0)*(vm[1]-vm[0])/10
       vi=gd.dp[sindo]/30; vi[vi>1.0]=1.0; vi[vi<0]=0; WSP[i,sindo]=vi; WSPn[i,sindo]=vi  #ocean
   
   #for POP
   for i in [11,12]:
       vm=[0.1,  0.2]; WSP[i] =vm[1];  WSP[i,fpz0] =vm[0];  WSP[i,fpz] =vm[0]+(gd.dp[fpz]-1.0)*(vm[1]-vm[0])/10
       vm=[0.01, 0.1]; WSPn[i]=vm[1]; WSPn[i,fpz0]=vm[0]; WSPn[i,fpz]=vm[0]+(gd.dp[fpz]-1.0)*(vm[1]-vm[0])/10
       vi=gd.dp[sindo]/30; vi[vi>1.0]=1.0; vi[vi<0]=0; WSP[i,sindo]=vi; WSPn[i,sindo]=vi  #ocean

   #zero velocity @upper river
   WSP[:,sindr]=0;    WSPn[:,sindr]=0

   #for KC0,KN0,KP0 in the ocean
   KC0=zeros([3,gd.np]); KC0[2]=0.05;  KC0[2,sindoo]=0.01
   KN0=zeros([3,gd.np]); KN0[2]=0.075; KN0[2,sindoo]=0.01
   KP0=zeros([3,gd.np]); KP0[2]=0.2;  KP0[2,sindoo]=0.01

   #write ICM_param
   S=zdata(); S.WSP=WSP; S.WSPn=WSPn; S.KC0=KC0; S.KN0=KN0; S.KP0=KP0; S.save(fname)

#----------------------------------------
#add ICM modules: SAV, CLAM, wetland
#----------------------------------------
if p.flag[fname]!=0 and (1 in [p.flag['CLAM'], p.flag['SAV'],p.flag['WET']]):
   if p.flag['CLAM']==1:
      S.cpatch0=M.f2(gd.exy,M.cpatch0).astype('int') 
      for m in ['cFc','cMTB']: S.attr(m,array([M.f2(gd.exy,i) for i in M.attr(m)]))
   if p.flag['SAV']==1: S.spatch0=M.f2(gd.exy,M.spatch0).astype('int'); S.sFc=M.f2(gd.exy,M.sFc)
   if p.flag['WET']==1: S.vpatch0=M.f2(gd.exy,M.vpatch0).astype('int'); S.vAw=M.f2(gd.exy,M.vAw)
   S.save(fname)

#----------------------------------------------------------------------
#WW3-frocing: boundary condition for WWM model
#----------------------------------------------------------------------
if p.flag['WWM']==1:
   print('writing WWM forcing files '); bfiles=[]; bname='bndfiles.dat'
   for yy in arange(num2date(p.StartT).year,num2date(p.EndT).year+1):
     for mm in arange(1,13):
         for svar in ['dir','fp','hs','spr','t02']:
             sname='WW3-GLOB-30M_{:04d}{:02d}_{}.nc'.format(yy,mm,svar); fname=p.WW3+'/{}/{}'.format(yy,sname)
             if not fexist(fname): continue
             if fexist(sname): os.remove(sname)
             os.symlink(relpath(fname),sname); bfiles.append(sname)
   [os.remove(bname) if fexist(bname) else None]; fid=open(bname,'w+'); fid.write('\n'.join(bfiles)); fid.close()

#----------------------------------------------------------------------
#sflux
#----------------------------------------------------------------------
fname='sflux'
if p.flag[fname]==1:
   print('writing '+fname)
   tm=[datenum(1979,1,1),datenum(2025,3,31)] #time range of database
   tdir='./sflux'; os.rmdir(tdir) if fexist(tdir) else None; os.mkdir(tdir) #mkdir sflux dir.

   #make links 
   mts=arange(max(tm[0],p.StartT-2),min(tm[1],p.EndT+2)); svars=['air','rad','prc']
   for irec,mt in enumerate(mts):
       for m,svar in enumerate(svars):
           t=num2date(mt); fname='{}/sflux_{}.{:04d}_{:02d}_{:02d}.nc'.format(p.sflux,svar,t.year,t.month,t.day)
           if len(mts)<10000:
              sname='{}/sflux_{}_1.{:04d}.nc'.format(tdir,svar,irec+1)
           else:
              sname='{}/sflux_{}_1.{:05d}.nc'.format(tdir,svar,irec+1)
           os.symlink('../'+relpath(fname),sname)
           if m==0 and irec%365==0: print('    sflux: {:04d}-{:02d}-{}'.format(t.year,t.month,t.day))
   fid=open('{}/sflux_inputs.txt'.format(tdir),'w+'); fid.write('&sflux_inputs\n   \n/'); fid.close()

#----------------------------------------------------------------------
#source.nc
#----------------------------------------------------------------------
fname='source.nc'

#common inputs
svars=['temp','salt','clay1','clay2','silt','sand','PB1','PB2','PB3','RPOC','LPOC','DOC',
       'RPON','LPON','DON','NH4','NO23','RPOP','LPOP','DOP','PO4','COD','DO','SRPOC','SRPON','SRPOP','PIP']

#partitioning coefficients for each river basin
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
fPB2  = [0.9,       0.9,          0.9,      0.9,       0.9,    0.9,   0.9,           0.9,       0.9]
fPB3  = [0.01,      0.01,         0.01,     0.01,      0.01,   0.01,  0.01,          0.01,      0.01]
#fPB2  = [0.9,       0.9,          0.1,      0.9,       0.9,    0.9,   0.9,           0.9,       0.1]
#fPB3  = [0.01,      0.01,         0.9,      0.01,      0.01,   0.01,  0.01,          0.01,      0.9]

#--------------------------------------------------------------------------
#Phase-7 WSM loading
#--------------------------------------------------------------------------
if p.flag[fname]==1 and ('p7' in p.source): 
   print('writing '+fname)
   #read data
   C=read(p.source,1)
   if hasattr(C,'wsm'): [C.attr(i,C.wsm.attr(i)) for i in C.wsm.attr()];  [C.attr('sho_'+i,C.sho.attr(i)) for i in C.sho.attr()] #old format

   #find source element
   ie=gd.ie(); sinde=unique(ie[near_pts(r_[C.sxy,C.sho_sxy],gd.exy[ie])]) #sinde is source element
   eid=near_pts(C.sxy,gd.exy[sinde]); C.eid=array([[eid[i] for i in sid] for sid in C.sid],'O')

   #-----------------------------------
   #watershed loading
   #-----------------------------------
   fpn=sort(pindex((C.time>=p.StartT)*(C.time<(p.EndT+2)))); mtime=C.time[fpn]
   nt,ntr,nsource=len(mtime),len(svars),len(sinde); flow=zeros([nt,nsource]); trs=zeros([nt,ntr,nsource])
   sflow=C.flow[C.river=='susquehanna'].sum(axis=0)[fpn] #susquehanna river flow

   for id, sname in enumerate(C.sname): #each watershed
       river=C.river[id];  rid=rivers.index(river) if (river in rivers) else 0
       for n,[rat,eid] in enumerate(zip(C.srat[id],C.eid[id])): #each element
          #flow, temp, salt
          flowi=C.flow[id,fpn]*rat
          temp=C.wtmp[id,fpn]*flowi*86.4;  DOX=C.doxx[id,fpn]*flowi*86.4 #treat the units as mg/L
          COD=zeros(nt); salt=zeros(nt);
          NH4=C.nh4x[id,fpn]*rat; NO3=C.no3x[id,fpn]*rat; PO4=C.po4x[id,fpn]*rat
          NH4[NH4<0]=0; NO3[NO3<0]=0; PO4[PO4<0]=0; temp[temp<0]=0; DOX[DOX<0]=0

          #compute algal carbon (mg[C]/L * flow)
          PB1=c2chl[0]*(fPB1[rid]*C.chla[id,fpn]*flowi*86.4); PB1[PB1<0]=0
          PB2=c2chl[1]*(fPB2[rid]*C.chla[id,fpn]*flowi*86.4); PB2[PB2<0]=0
          PB3=c2chl[2]*(fPB3[rid]*C.chla[id,fpn]*flowi*86.4); PB3[PB3<0]=0

          #compute nitrogen
          ORGN=C.orgn[id,fpn]*rat-n2c[0]*PB1-n2c[1]*PB2-n2c[2]*PB3; ORGN[ORGN<0]=0
          DON=fDON[rid]*ORGN;   PON=ORGN-DON
          RPON=fRPON[rid]*PON;  LPON=fLPON[rid]*PON;  SRPON=PON-RPON-LPON

          #compute carbon
          ORGC=c2n[rid]*(DON+PON); DOC=fDOC[rid]*ORGC; POC=ORGC-DOC
          RPOC=fRPOC[rid]*POC; LPOC=fLPOC[rid]*POC; SRPOC=POC-RPOC-LPOC

          #compute phosphorus
          ORGP=C.orgp[id,fpn]*rat+C.pipx[id,fpn]*rat-p2c[0]*PB1-p2c[1]*PB2-p2c[2]*PB3; ORGP[ORGP<0]=0
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
          sand=C.sand[id,fpn]*rat; silt=C.silt[id,fpn]*rat; clay1=0.0*C.clay[id,fpn]*rat; clay2=1.0*C.clay[id,fpn]*rat

          #add load
          SRPOC[SRPOC<0]=0; SRPON[SRPON<0]=0; SRPOP[SRPOP<0]=0 #make SRPOM conc not negative
          flow[:,eid]=flow[:,eid]+flowi
          trs[:,:,eid]=trs[:,:,eid]+c_[temp,salt,clay1,clay2,silt,sand,PB1,PB2,PB3,RPOC,LPOC,DOC,
                                   RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4,COD,DOX,SRPOC,SRPON,SRPOP,PIP]

   #--------------------------------------------------------------------------
   #add SHO input
   #--------------------------------------------------------------------------
   evars0=['clay1','clay2','silt', 'sand', 'RPOC', 'SRPOC', 'RPON', 'SRPON', 'RPOP', 'SRPOP', 'PIP']
   evars= ['clay' ,'clay' ,'silt', 'sand', 'RPOC', 'G3OC',  'RPON', 'G3ON',  'RPOP', 'G3OP',  'PIP']
   fps=sort(pindex((C.sho_time>=p.StartT)*(C.sho_time<(p.EndT+2)))); fpt=near_pts(C.sho_time[fps],mtime)
   for i, id in enumerate(near_pts(C.sho_sxy,gd.exy[sinde])):
       for evar0,evar in zip(evars0,evars):
           rat=0.5 if evar=='clay' else 1
           itr=svars.index(evar0); trs[fpt,itr,id]=trs[fpt,itr,id]+rat*C.attr('sho_'+evar)[i,fps]

   #--------------------------------------------------------------------------
   #preparing the model input
   #--------------------------------------------------------------------------
   fpn=flow.mean(axis=0)!=0; flow=flow[:,fpn]; trs=trs[...,fpn]; sinde=sinde[fpn]; nsource=len(sinde) #remove zero flows

   #compute conc.
   flow[flow<1.e-8]=1.e-8; trs=trs/(86.4*flow[:,None,:])
   tri=trs[:,0];  tri[tri>40]=40; trs[:,0]=tri   #limit temp value
   tri=trs[:,-5]; tri[tri>20]=20; trs[:,-5]=tri  #limit DO value
   trs[:,2:6]=trs[:,2:6]/1e3                     #convert sediment concentration to g/L

   #check modules
   if p.flag['SED']==0: ntr=ntr-4; tid=r_[arange(2),arange(6,27)]; trs=trs[:,tid,:] #remove SED setup
   if p.flag['ICM']==0: ntr=ntr-21; trs=trs[:,:-21,:]                               #remove ICM setup
   if p.flag['ICM'] in [2,20]: ntr=ntr-4; trs=trs[:,:-4,:]                          #remove ICM's CBPO setup

   #write source.nc
   C=zdata();C.dimname=['nsources','nsinks','ntracers','time_msource','time_vsource','time_vsink','one']; C.dims=[nsource,0,ntr,nt,nt,0,1]
   for i in ['time_step_vsource','time_step_vsink','time_step_msource']: add_var(C,i,('one',),array(86400.0))
   add_var(C,'source_elem',('nsources',),sinde+1); add_var(C,'vsource',('time_vsource','nsources'),flow.astype('float32'))
   add_var(C,'msource',('time_msource','ntracers','nsources'),trs.astype('float32'))
   C.save('source.nc',zlib=True); C,S,flow,trs=None,None,None,None 

#--------------------------------------------------------------------------
#Phase-6 WSM loading
#--------------------------------------------------------------------------
if p.flag[fname]==1 and ('p6' in p.source): 
   print('writing '+fname)
   #--------------------------------------------------------------------------
   #input for NPS and PS loading
   #--------------------------------------------------------------------------
   #read data
   C=loadz(p.source); nps,ps,sho=C.nps,C.ps,C.sho
   lexy=zeros([gd.ne,2]); lexy[gd.fp3]=gd.lxy[gd.elnode[gd.fp3,:3]].mean(axis=1); lexy[gd.fp4]=gd.lxy[gd.elnode[gd.fp4]].mean(axis=1)

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

   #pre-proc
   mtime=arange(p.StartT,p.EndT); nt,ntr,nsource=len(mtime),len(svars),len(C.sinde)
   fpn=sort(nonzero((nps.time>=mtime[0])*(nps.time<=mtime[-1]))[0])  #time flag for NPS
   fpp=sort(nonzero((ps.time>=mtime[0])*(ps.time<=mtime[-1]))[0])    #time flag for PS
   fps=sort(nonzero((sho.time>=mtime[0])*(sho.time<=mtime[-1]))[0])  #time flag for SHO
   sinde=near_pts(C.xy,lexy) #grid elements

   flow=zeros([nt,nsource]); trs=zeros([nt,ntr,nsource])

   #collect sources from NPS
   nps.rivers=array([C.rivers[list(C.RiverSeg).index(i[7:] if len(i.strip())==20 else i.strip())] for i in nps.sname])
   sflow=nps.flow[nps.rivers=='susquehanna'].sum(axis=0)[fpn] #susquehanna river flow

   for k,[rats,ids] in enumerate(zip(C.nrat,C.nid)):   #each grid elements
       if len(ids)==0: continue                        #no NPS into this elements
       for n,[rat,id] in enumerate(zip(rats,ids)):     #each nps
           river=nps.rivers[id]; rid=rivers.index(river) #river basin

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

   #--------------------------------------------------------------------------
   #preparing the model input
   #--------------------------------------------------------------------------
   #compute conc.
   flow[flow<1.e-8]=1.e-8; trs=trs/(86.4*flow[:,None,:])
   tri=trs[:,0]; tri[tri>40]=40; trs[:,0]=tri #limit temp value
   tri=trs[:,-5]; tri[tri>20]=20; trs[:,-5]=tri  #limit DO value
   trs[:,2:6]=trs[:,2:6]/1e3 #convert sediment concentration to g/L

   #check modules
   if p.flag['SED']==0: ntr=ntr-4; tid=r_[arange(2),arange(6,27)]; trs=trs[:,tid,:] #remove SED setup
   if p.flag['ICM']==0: ntr=ntr-21; trs=trs[:,:-21,:] #remove ICM setup
   if p.flag['ICM'] in [2,20]: ntr=ntr-4; trs=trs[:,:-4,:] #remove ICM's CBPO setup

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
   WriteNC('source.nc',nd,zlib=True)

#----------------------------------------------------------------------
#change grid projection for WWM
#----------------------------------------------------------------------
if p.flag['WWM']==1: os.remove('hgrid.gr3'); os.symlink('hgrid.ll','hgrid.gr3')
