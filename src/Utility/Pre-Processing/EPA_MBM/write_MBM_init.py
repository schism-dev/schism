#!/usr/bin/env python3
from pylib import *
close("all")

#read grid
gd,vd=grd('grid.npz',2)
C=zdata(); C.gd=gd; C.vd=vd
#C=read('MBM_init.npz')

#----------------------------------------------------
#ICM
#----------------------------------------------------
#read all data
C0=ReadNC('hotstart.nc_174240',1)

#variables to be saved
bvars=['bPOC','bPON','bPOP','btemp','bstc','bSTR','bThp','bTox','bNH4','bNH4s','bNO3','bPO4','bH2S','bCH4']
ivars=['PB1','PB2','PB3','RPOC','LPOC','DOC','RPON','LPON','DON','NH4','NO3','RPOP','LPOP','DOP','PO4','COD','DOX','SRPOC','SRPON','SRPOP','PIP']

#get variable values for ICM
for i in bvars: C.attr(i,C0.attr(i).astype('float32').T)
for i,k in enumerate(ivars): C.attr(k,C0.tr_nd.value[...,i+2].astype('float32').T)
C0.close()

#----------------------------------------------------
#SED
#----------------------------------------------------
#read all data
C0=ReadNC('hotstart.nc_SED',1)

#variables to be saved
svars=['SED3D_dp','SED3D_rough','SED3D_bed','SED3D_bedfrac']

#get variable values for ICM
for i in svars: C.attr(i,squeeze(C0.attr(i).astype('float32')).T)
C0.close()

#----------------------------------------------------
#clam
#----------------------------------------------------
gd0=grd('grid_LR.npz')

#read init data
C0=ReadNC('hotstart.nc_clam',1)
C.clam=array([gd0.interp(gd.exy,value=i).astype('float32') for i in C0.clam.T]) #save clam
C0.close()

#read param
C0=ReadNC('ICM_param.nc_clam',1)
C.cpatch0=gd0.interp(gd.exy,value=C0.cpatch0).astype('float32') #save patch
for i in ['cFc','cMTB']: C.attr(i, array([gd0.interp(gd.exy,value=k).astype('float32') for k in C0.attr(i)]))
C0.close()

#----------------------------------------------------
#sav
#----------------------------------------------------
gd0=grd('grid_LR.npz')

#read init data
C0=ReadNC('hotstart.nc_sav',1)
C.sav=array([gd0.interp(gd.exy,value=i).astype('float32') for i in C0.sav.T]) #save sav
C.EP=gd0.interp(gd.exy,value=C0.EP).astype('float32') #save EP
C0.close()

#read param
C0=ReadNC('ICM_param.nc_sav',1)
C.spatch0=gd0.interp(gd.exy,value=C0.spatch0).astype('float32') #save patch
C.sFc=gd0.interp(gd.exy,value=C0.sFc).astype('float32') #save patch
C0.close()

#----------------------------------------------------
#marsh
#----------------------------------------------------
gd0=grd('grid_LR.npz')

#read param
C0=ReadNC('ICM_param.nc_marsh',1)
C.vpatch0=gd0.interp(gd.exy,value=C0.vpatch0).astype('float32') #save patch
C.vAw=gd0.interp(gd.exy,value=C0.vAw).astype('float32') #save patch
C0.close()

#----------------------------------------------------
#rough.gr3, and watertype
#----------------------------------------------------
C.rough=read('rough.gr3').interp(gd.xy)
C.watertype=read('watertype.gr3').interp(gd.xy).astype('int')

#----------------------------------------------------
#for temp
#----------------------------------------------------
C0=read('temp.npz',1)
mti=C0.time+datenum(1991,1,1); doy=array([i-datenum(num2date(i).year,1,1) for i in mti]); mts=[]; mys=[]
for month in arange(1,13):
    year=2000; t0=datenum(year,1,1); t1=datenum(year,month,1)-t0; t2=datenum(year,month+1,1)-t0
    fp=(doy>=t1)*(doy<t2); mts.append(doy[fp].mean())
    mys.append([gd.interp(i).astype('float32') for i in C0.temp_elem.value[fp].mean(axis=0)]) 
C.temp_time=array(mts); C.temp_data=array(mys)

#----------------------------------------------------
#for salinity
#----------------------------------------------------
#model data
C0=read('salt.npz',1)
mti=C0.time+datenum(1991,1,1); mts=[]; mys=[]
for year in arange(num2date(mti[0]).year,num2date(mti[-1]).year+1):
    print(year)
    for month in arange(1,13):
        t1=datenum(year,month,1); t2=datenum(year,month+1,1)
        fp=(mti>=t1)*(mti<t2); mts.append(mti[fp].mean())
        mys.append([gd.interp(i).astype('float32') for i in C0.salt_elem.value[fp].mean(axis=0)]) 
C.salt_time=array(mts); C.salt_data=array(mys)

#climatology data
C0=read('ChesBay_Salinity_climatology.nc',1)
C.csalt_lon=C0.lon[:]; C.csalt_lat=C0.lat[:]; C.csalt_depth=C0.depth[:]
C.csalt_time=C0.time[:]; C.csalt_data=C0.salt[:].astype('float32')

#----------------------------------------------------
#for elev
#----------------------------------------------------
stations=[8656483,8557380]
dt=1/24 #time step of elevation, to be consistent with dt in elev2D.th.nc
fc=0.5  #cutoff frequency fc for low-pass filter

#read database
E=loadz('/sciclone/data10/wangzg/CBP/database/noaa_elev_msl.npz') 
S=loadz('/sciclone/data10/wangzg/CBP/database/elev_adjust_v0.npz')   

elev=[]; mti=arange(datenum(1979,1,2),datenum(2025,3,31),dt)
for m,station in enumerate(stations):
    fp=E.station==station; cti,cyi=E.time[fp],E.elev[fp]; #shift by 0.35 hr
    eyi=array([interp(cti+S.pha[i,m]/(4*pi)-0.01458,cyi*S.amp[i,m]*S.weight[i,m],mti) for i in arange(gd.nobn[0])])
    elev.append(eyi)
elev=array(elev).sum(axis=0).T
elev=1.2*elev.astype('float32') #adjust elev.
C.elev_time=mti; C.elev_xy=gd.xy[gd.iobn[0]]; C.elev_data=elev.T

#----------------------------------------------------
#ICM in ocean and bay mouth
#----------------------------------------------------
#ocean
C0=read('NEFSC_climatology.npz',1)
C.NEFSC_time=C0.time[:]*1.0; C.NEFSC_depth=C0.depth[:]*1.0
C.NEFSC_DO=C0.DO[:]; C.NEFSC_NH4=C0.NH4[:]; C.NEFSC_NO3=C0.NO23[:]
C.NEFSC_PO4=C0.PO4[:]; C.NEFSC_SiO4=C0.SiO4[:]

#bay mouth
C0=read('CBP_WQData.npz',1) #input for get obs. data @CB7.4
svars=['CHLA','PC', 'DOC', 'PN',  'NH4F','NO23F','PP', 'PO4F','DO','DON','DOP']
rvars=['CHLA','POC','DOC', 'PON', 'NH4', 'NO3',  'POP','PO4', 'DO','DON','DOP']
svms =[ 50,    5,    5,     0.5,   0.1,  0.3,     0.2,  0.05,  15,  0.4,  0.4] #maximum
mti=arange(datenum(1988,1,1),datenum(2024,9,30)); mdepth=arange(18)*1.0 #time and depth range
for svar, rvar, svm in zip(svars,rvars,svms):
    fp=(C0.station=='CB7.4')*(C0.var==svar)*(C0.data<svm); ots,oys,ods=C0.time[fp],C0.data[fp],C0.depth[fp] #all data
    uts=unique(ots); uys=[] 
    for ti in uts: #interp in vertical
        fp=abs(ots-ti)<=1; od,oy=ods[fp],oys[fp] #get data around time
        od,sind=unique(od,return_index=True); oy=oy[sind] #unique depth
        fy=ones(len(mdepth))*oy if len(od)==1 else interp(od,oy,mdepth); uys.append(fy)
    data=interp(uts,uys,mti,axis=0).astype('float32'); C.attr('CB74_'+rvar,data.T)
fpc=(mti>=datenum(1996,1,1))*(mti<=datenum(2007,1,1)); C.CB74_DOC[:,fpc]=11.25*C.CB74_DON[:,fpc] #fix DOC gap
C.CB74_time=mti; C.CB74_depth=mdepth

#----------------------------------------------------
#salinity @ET2.1
#----------------------------------------------------
C0=read('CBP_WQData.npz',1)
fp=(C0.var=='SALINITY')*(C0.station=='ET2.1')*(C0.layer=='S')
C.ET21_time,C.ET21_salt=sort_all(C0.time[fp],C0.data[fp])

#----------------------------------------------------
#station.bp
#----------------------------------------------------
C0=read('station_CB.bp'); C.station_xyz=C0.xyz; C.station_name=C0.station

#----------------------------------------------------
#hydro_veg module
#----------------------------------------------------
gd0=grd('grid_sav.npz'); C0=read('sav_region.npz',1)
ev=zeros(gd0.ne); ev[C0.isav]=1; pv=gd0.interp(ev,fmt=2); 
ec=zeros(gd0.ne); ec[C0.isav]=C0.fc; pc=gd0.interp(ec)

sindp=pindex(gd0.interp(gd.xy,pv)>0); fc=gd0.interp(gd.xy,pc)
C.veg_h=zeros(gd.np); C.veg_D=zeros(gd.np); C.veg_N=zeros(gd.np); C.veg_cd=ones(gd.np)
C.veg_h[sindp]=0.5; C.veg_D[sindp]=0.01; C.veg_N[sindp]=(50*fc[sindp]).astype('int')

#----------------------------------------------------
#save data
#----------------------------------------------------
C.save('MBM_init')
