#!/usr/bin/env python3
from pylib import *
close('all')

#----------------------------------------------------------------------
#inputs
#----------------------------------------------------------------------
grd='MBM_v9b.2dm'                               #SMS grid
prj='epsg:26918'                                #projection of SMS grid
bxy=[-76.485900,-75.063539,34.599849,38.628769] #lon&lat of bnd nodes
headers=('etopo1','ncei19_CB','crm_3arcs','BlueTopo','cb_ll','cb_bay_dem_v3.1_ll') #DEM headers
regions=['RappChannel.reg','CD_canal.reg']      #regions of setting mini. depth
rvalues=[16,10]                                 #mini. depth in each region
reg_max='DEM_maximum_v9b.reg'                   #use maximum depth of all DEMs in region
bdir=r'/sciclone/data10/wangzg/CBP/grid/dem'    #directory of DEM data
sdir='/sciclone/data10/wangzg/CBP/database'     #ChesBay setup dir

#----------------------------------------------------------------------
#get hgrid.gr3, hgrid.ll and vgrid.in
#----------------------------------------------------------------------
#read 2dm, do projection, compute bndfile, and split bad quads
print('converting *2dm to *.gr3')
gd=sms2grd(grd); gd.lon,gd.lat=gd.proj(prj,fmt=1); gd.x,gd.y=gd.proj(prj,'epsg:26918',fmt=1) #projection 
gd.split_quads(angle_min=60,angle_max=120) #split bad quads
bxy=array(bxy); bx,by=proj_pts(bxy[:2],bxy[2:],'epsg:4326','epsg:26918'); gd.compute_bnd([*bx,*by]) #define boundary

#load bathymetry
print('loading bathymetry to grid')
fpm=read(sdir+'/region/'+reg_max).inside(gd.xy,fmt=1)
for header in headers: 
    DEMs=array([i for i in os.listdir(bdir) if i.startswith(header)])
    for n,DEM in enumerate(DEMs): 
        print('   DEM: {}, {}/{}'.format(DEM,n+1,len(DEMs)))
        dp,ip=load_bathymetry(gd.lon,gd.lat,'{}/{}'.format(bdir,DEM),fmt=1); dp=-dp #depth from new DEM
        fpz=fpm[ip]*(gd.dp[ip]>dp); dp[fpz]=gd.dp[ip[fpz]]  #impose maximum DEM depth in region
        gd.dp[ip]=dp
gd.dp=gd.dp-0.258 #from navd to ngvd (msl)

#minimum depth in regions, etc.  (todo: regions)
print("apply mininum depth")
gd.dp[gd.dp<0.5]=0.5
for region,dpm in zip(regions,rvalues):
    fname=sdir+'/region/'+region
    if not os.path.exists(fname): continue
    gd.z[read(fname).inside(gd.lxy,1)*(gd.dp<dpm)]=dpm

#save hgrid
gd.save('hgrid.gr3',fmt=1); gd.save('hgrid.ll',fmt=1,xy=1) #; gd.x,gd.y=px,py; gd.save('hgrid_bk.gr3',fmt=1)

#get vgrid
print('generating vgrid.in')
code='ifx -Bstatic -O2 -CB -o gen_vqs gen_vqs.f90 schism_geometry.f90'
v=zeros(gd.np); v[read(sdir+'/region/estuary.reg').inside(gd.lxy)]=1; gd.save('estuary.gr3',value=v)
os.system('ln -sf {}/vgrid/*.f90 ./; {}; ./gen_vqs'.format(sdir,code)) #gen vgrid.in
read('vgrid.in').write_vgrid() #change vgrid format

#clean folders
#os.system('mv hgrid.gr3 hgrid.ll; mv hgrid_bk.gr3 hgrid.gr3') #rename hgrid
os.system('rm bad_quad.bp estuary.gr3 fort.* gen_vqs* schism_geometry* vgrid_master_*')
save_schism_grid()
