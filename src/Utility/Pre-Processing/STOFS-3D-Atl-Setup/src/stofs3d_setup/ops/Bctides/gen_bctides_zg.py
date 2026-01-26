#!/usr/bin/env python3
#generate bctides.in's tide harmonics
from pylib import *

#---------------------------------------------------------------------
#input
#---------------------------------------------------------------------
tnames=['S2','P1','N2','K1','M2','Q1','K2','O1']
StartT=[2022,1,1,0]  #year,month,day,hour
nday=171    #number of days
ibnds=[1,]           #order of open boundaries (starts from 1)
flags=[[5,5,4,4],]   #SCHISM bnd flags for each boundary
Z0=0.0               #add Z0 constant if Z0!=0.0

grd='./hgrid.ll'       #hgrid.ll (includes bndinfo), or grid.npz (include lon,lat)
bdir=r'/sciclone/data10/wangzg/FES2014'   #FES2014 database 

#---------------------------------------------------------------------
#read bndinfo, amp, freq, nodal factor and tear
#---------------------------------------------------------------------
#read grid information
if grd.endswith('.npz'):
   gd=loadz(grd).hgrid; gd.x=gd.lon; gd.y=gd.lat
else: 
   gd=read_schism_hgrid(grd)

#get tidal amplitude and frequency
amp=[]; freq=[]; ts=loadz('{}/tide_fac_const/tide_fac_const.npz'.format(bdir))
for tname in tnames: 
    sind=nonzero(ts.name==tname.upper())[0][0]
    amp.append(ts.amp[sind]); freq.append(ts.freq[sind])

#get nodal factor
tdir='{}/tide_fac_improved'.format(bdir)
fid=open('./tide_fac.in','w+'); fid.write('{}\n{} {} {} {}\n0\n'.format(nday,*StartT[::-1])); fid.close()
print('ifort -o tide_fac_improved {}/tf_main.f90 {}/tf_selfe.f90; ./tide_fac_improved < tide_fac.in'.format(tdir,tdir))
os.system('ifort -o tide_fac_improved {}/tf_main.f90 {}/tf_selfe.f90; ./tide_fac_improved < tide_fac.in'.format(tdir,tdir))

nodal=[]; tear=[]  #read nodal factor 
for tname in tnames:
    lines=[i for i in open('./tide_fac.out','r').readlines() if len(i.split())==3]
    line=[i for i in lines if i.strip().startswith(tname.upper())][0]
    nodal.append(float(line.strip().split()[1]))
    tear.append(float(line.strip().split()[2]))
os.system('rm tide_fac_improved tide_fac.in tide_fac.out')

#---------------------------------------------------------------------
#write bctides.in
#---------------------------------------------------------------------
fid=open('bctides.in','w+')
fid.write('!{:02}/{:02}/{:4} {:02}:00:00 UTC\n'.format(*array(StartT)[array([1,2,0,3])]))

#tidal potential, and frequency
fid.write(' {:d}  50.000 !number of earth tidal potential, cut-off depth for applying tidal potential\n'.format(len(tnames)))
for i,tname in enumerate(tnames): fid.write('{}\n{} {:<.6f}  {:<.9e}  {:7.5f}  {:.2f}\n'.format(tname,tname[1],amp[i],freq[i],nodal[i],tear[i])) 
fid.write('{} !nbfr\n'.format(len(tnames)+int(Z0!=0)))
if Z0!=0: fid.write('Z0\n  0.0 1.0 0.0\n')
for i,tname in enumerate(tnames): fid.write('{}\n  {:<.9e}  {:7.5f}  {:.2f}\n'.format(tname,freq[i],nodal[i],tear[i])) 
fid.write('{} !nope\n'.format(gd.nob))

#write tidal harmonic for each boundary
for i,ibnd in enumerate(ibnds):
    fstr='{} '+'{} '*len(flags[i])+'!ocean\n'
    fid.write(fstr.format(gd.nobn[ibnd-1],*flags[i]))

    #get bundary lon&lat, and interpolate for the amp and pha
    nobn=gd.nobn[ibnd-1]; sind=gd.iobn[ibnd-1]; xi=mod(gd.x[sind]+360,360); yi=gd.y[sind]; ap=[]
    for m,tname in enumerate(tnames): 
        print('compute amp. and pha. for tide {} of boundary: {}'.format(tname,i+1))
        api=[]
        for n in arange(3):  
            if n==0: fname='{}/fes2014b_elevations_extrapolated/ocean_tide_extrapolated/{}.nc'.format(bdir,tname.lower()); an,pn='amplitude','phase'
            if n==1: fname='{}/eastward_velocity/{}.nc'.format(bdir,tname.lower()); an,pn='Ua','Ug'
            if n==2: fname='{}/northward_velocity/{}.nc'.format(bdir,tname.lower());an,pn='Va','Vg'
            C=ReadNC(fname,1); lon=array(C.variables['lon'][:]); lat=array(C.variables['lat'][:])
            amp0=array(C.variables[an][:])/100; pha0=array(C.variables[pn][:]); C.close()
            fpn=pha0<0; pha0[fpn]=pha0[fpn]+360
            dxs=unique(diff(lon)); dys=unique(diff(lat))
            if len(dxs)!=1 or len(dys)!=1: sys.exit('{}: lon,lat not uniform specified'.format(fname))
            dx=dxs[0]; dy=dys[0]

            #get interp index (todo, if lon&lat not uniformal interval, do loop to find the right index)
            idx=floor((xi-lon[0])/dx).astype('int'); sind=nonzero((lon[idx]-xi)>0)[0]; idx[sind]=idx[sind]-1
            idy=floor((yi-lat[0])/dy).astype('int'); sind=nonzero((lat[idy]-yi)>0)[0]; idy[sind]=idy[sind]-1 
            xrat=(xi-lon[idx])/(lon[idx+1]-lon[idx]); yrat=(yi-lat[idy])/(lat[idy+1]-lat[idy])
            if sum((xrat>1)|(xrat<0)|(yrat>1)|(yrat<0))!=0: sys.exit('xrat or yrat >1 or <0')

            #interp for amp,pha
            apii=[]
            for k in arange(2): 
                if k==0: v0=c_[amp0[idy,idx],amp0[idy,idx+1],amp0[idy+1,idx],amp0[idy+1,idx+1]].T; vm=100
                if k==1: v0=c_[pha0[idy,idx],pha0[idy,idx+1],pha0[idy+1,idx],pha0[idy+1,idx+1]].T; vm=370
                vmax=v0.max(axis=0); vmin=v0.min(axis=0)
                if k==1: #deal with phase jump
                   for kk in nonzero(abs(vmax-vmin)>180)[0]: 
                       fpn=abs(v0[:,kk]-vmax[kk])>180; v0[fpn,kk]=v0[fpn,kk]+360
                v1=v0[0]*(1-xrat)+v0[1]*xrat; v2=v0[2]*(1-xrat)+v0[3]*xrat; apiii=v1*(1-yrat)+v2*yrat
                sind=nonzero((vmax>vm)*(vmin<=vm)*(vmin>=0))[0]; apiii[sind]=vmin[sind]
                if sum((vmax>vm)*((vmin>vm)|(vmin<0)))!=0: sys.exit('All junks for amp or pha')
                apii.append(apiii)
            api.append(apii)
        ap.append(api)
    ap=array(ap).transpose([1,0,3,2])

    #write tidal amp and pha for elev 
    if Z0!=0: fid.write('Z0\n'); [fid.write('{} 0.0\n'.format(Z0)) for i in arange(nobn)]
    for m,tname in enumerate(tnames):  
        fid.write('{}\n'.format(tname.lower()))
        for k in arange(nobn):
            fid.write('{:8.6f} {:.6f}\n'.format(*ap[0,m,k]))

    #write tidal amp and pha for uv
    if Z0!=0: fid.write('Z0\n'); [fid.write('0.0 0.0 0.0 0.0\n') for i in arange(nobn)]
    for m,tname in enumerate(tnames): 
        fid.write('{}\n'.format(tname.lower()))
        for k in arange(nobn):
            fid.write('{:8.6f} {:.6f} {:8.6f} {:.6f}\n'.format(*ap[1,m,k],*ap[2,m,k]))
fid.close()
