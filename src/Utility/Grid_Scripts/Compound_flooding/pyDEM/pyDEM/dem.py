#!/usr/bin/env pythonw
#this is the python library used to extracted river network from DEM
import richdem as rd
from osgeo import gdal
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, MultiPoint
from shapely.ops import split as ssplit

from pylib import *

class dem(object):
    def __init__(self):
        pass

    def read_data(self,path,svar=None):
        '''
        read data of *.npz format
        '''
        if svar is None:
            #read data
            data=loadz(path)

            #get variable name
            svars=list(data.__dict__.keys())
            if 'VINFO' in svars: svars.remove('VINFO')

            #put variables in self
            for svar in svars:
                exec('self.{}=data.{}'.format(svar,svar))
        else:
            fid=load(path,allow_pickle=True)
            data=fid[svar]
            fid.close()
            return data[()]

    def save_data(self,fname,svars):
        '''
        save attributes of svars in *.npz format
        '''
        if not isinstance(svars,list): svars=[svars,]

        S=zdata();
        for svar in svars:
            exec('S.{}=self.{}'.format(svar,svar))
        save_npz(fname,S)

    def read_demdata(self,outname='dem',data=None,save_bnd=True):
        '''
        read dem raw data from info.name or data(global)
        also save values on subdomain boundary for later combination
        '''
        ym,xm=self.info.ds
        #read data directly
        if data is None:
            if self.info.name.endswith('.asc'):
                self.data=loadtxt(self.info.name,skiprows=self.info.skiprows,max_rows=self.info.max_rows,usecols=self.info.usecols)
                self.__dict__[outname]=self.__dict__.pop('data')
            elif self.info.name.endswith('.tif'):
                #breakpoint()
                self.data = rd.LoadGDAL(self.info.name, no_data=-9999)
                self.__dict__[outname]=self.__dict__.pop('data')
            else:
                sys.exit('wrong format')
        else:
            #extract data
            if self.info.is_subdomain:
                iy,ix=self.info.ixy_global
            else:
                iy=0; ix=0
            eind=[iy,iy+ym,ix,ix+xm]
            exec('self.{}=data[{}:{},{}:{}].copy()'.format(outname,*eind))

        #convert nan to nodata
        self.remove_nan(outname)

        if type(self.dem) is not rd.rdarray:
            self.dem=rd.rdarray(self.dem, no_data=self.info.nodata)
        #only include data in depth_limit
        exec('fp=(self.{}<self.info.depth_limit[0])|(self.{}>self.info.depth_limit[1]); self.{}[fp]=self.info.nodata'.format(outname,outname,outname))

        #add nodata bnd
        sind00=nonzero(self.dem.ravel()==self.info.nodata)[0]
        if len(sind00)!=0:
            ns=min(int(1e6),sind00.size); nsub=int(sind00.size/ns); sindn=[];
            for i in arange(nsub):
                if i==(nsub-1):
                    sind0=sind00[arange(i*ns,sind00.size)]
                else:
                    sind0=sind00[arange(i*ns,(i+1)*ns)]
                #get all neighbors
                iy0,ix0=unravel_index(sind0,self.info.ds)
                yind=r_[iy0,  iy0-1, iy0,  iy0+1, iy0-1,iy0-1, iy0+1, iy0+1]
                xind=r_[ix0+1,ix0,   ix0-1,ix0,   ix0+1,ix0-1, ix0-1, ix0+1]
                #true neighbors
                fpt=nonzero((xind>=0)*(xind<=(xm-1))*(yind>=0)*(yind<=(ym-1)))[0]
                sind=ravel_multi_index([yind[fpt],xind[fpt]],self.info.ds)

                #nodata bnd
                fpn=nonzero(self.dem.ravel()[sind]!=self.info.nodata)[0]
                sindn.extend(unique(sind[fpn]))
            sindn=unique(array(sindn).astype('int'))
            self.info.sind_bnd=unique(r_[self.info.sind_bnd,sindn])
            self.info.sind_bnd_nodata=unique(r_[self.info.sind_bnd_nodata,sindn])

        #save bnd values
        if hasattr(self,'dem'):
            self.info.dem_bnd=self.dem.ravel()[self.info.sind_bnd].copy()
            if hasattr(self,'dir'): self.info.dir_bnd=self.dir.ravel()[self.info.sind_bnd].copy()

            #exclude nodata
            fpn=self.info.dem_bnd!=self.info.nodata
            self.info.sind_bnd=self.info.sind_bnd[fpn].copy()
            self.info.dem_bnd=self.info.dem_bnd[fpn].copy()
            if hasattr(self,'dir'): self.info.dir_bnd=self.info.dir_bnd[fpn].copy()


    def collect_subdomain_data(self,name='dem',outname='dem',collect_bndinfo=True):
        '''
        collect subdomain values to reconstruct global domain data
        name: attribute to be combined
        outname: name of combined global variable
        '''

        #initialize
        exec('self.{}=ones(self.info.ds)*self.info.nodata'.format(outname))

        #colloect data and also bnd dir and dem
        sind_bnd_local=[]; dem_bnd_local=[]; dir_bnd_local=[]; sind_bnd_nodata=[]
        for i in arange(self.info.nsubdomain):
            iy,ix=self.domains[i].info.ixy_global; ym,xm=self.domains[i].info.ds
            eind=[iy,iy+ym,ix,ix+xm]
            exec('self.{}[{}:{},{}:{}]=self.domains[{}].{}'.format(outname,*eind,i,name))

            #bnd info
            if collect_bndinfo:
                biy,bix=unravel_index(self.domains[i].info.sind_bnd, self.domains[i].info.ds)
                sind_bnd=ravel_multi_index([biy+iy,bix+ix],self.info.ds)
                sind_bnd_local.extend(sind_bnd)
                if hasattr(self.domains[i].info,'dir_bnd'): dir_bnd_local.extend(self.domains[i].info.dir_bnd)
                if hasattr(self.domains[i].info,'dem_bnd'): dem_bnd_local.extend(self.domains[i].info.dem_bnd)
                if len(self.domains[i].info.sind_bnd_nodata)!=0:
                    biy,bix=unravel_index(self.domains[i].info.sind_bnd_nodata, self.domains[i].info.ds)
                    sind_bnd=ravel_multi_index([biy+iy,bix+ix],self.info.ds)
                    sind_bnd_nodata.extend(sind_bnd)

        #assign bnd value
        if collect_bndinfo:
            self.info.sind_bnd_local=array(sind_bnd_local)
            self.info.dir_bnd_local=array(dir_bnd_local)
            self.info.dem_bnd_local=array(dem_bnd_local)
            self.info.sind_bnd_nodata=unique(r_[self.info.sind_bnd_nodata,sind_bnd_nodata])

    def collapse_subdomain(self,svars=['dem','dir']):
        '''
        this function is used to remove outer zone of subdomain
        update on "info"
        '''

        if not self.info.is_subdomain: return

        #information of subdomain
        ym0,xm0=self.info.ds;     iy0,ix0=self.info.ixy_global
        ym,xm=self.info.ds_local; iy,ix=self.info.ixy_local
        yll,xll,dxy=self.info.header[2:].copy()

        #extract the domain data
        for svar in svars:
            if not hasattr(self,svar): continue
            exec('self.{}=self.{}[{}:{},{}:{}].copy()'.format(svar,svar,iy,iy+ym,ix,ix+xm))

        #update subdomain info
        self.info.header=[ym,xm,yll-iy*dxy,xll+ix*dxy,dxy]
        self.info.extent=[xll+ix*dxy,xll+ix*dxy+(xm-1)*dxy,yll-iy*dxy-(ym-1)*dxy,yll-iy*dxy]
        self.info.ds=[ym,xm]
        self.info.skiprows=int(self.info.skiprows+iy)
        self.info.usecols=self.info.usecols[(ix0+ix):(ix0+ix+xm)]
        self.info.max_rows=int(ym)

        #relate to global domain
        self.info.ixy_local=[0,0]
        self.info.ds_lcoal=[ym,xm]
        self.info.ixy_global=[iy0+iy,ix0+ix]

        #local boundary index
        biy=r_[zeros(xm),ones(xm)*(ym-1),arange(1,ym-1),arange(1,ym-1)].astype('int')
        bix=r_[arange(xm),arange(xm),zeros(ym-2),ones(ym-2)*(xm-1)].astype('int')
        sind_bnd=ravel_multi_index([biy,bix],self.info.ds)
        if hasattr(self,'dem'): dem_bnd=self.dem.ravel()[sind_bnd].copy()
        if hasattr(self,'dir'): dir_bnd=self.dir.ravel()[sind_bnd].copy()

        #if there were bnd inside new domain
        if len(self.info.sind_bnd)!=0:
            biy,bix=unravel_index(self.info.sind_bnd,[ym0,xm0])
            fpb=nonzero((biy>iy)*(biy<(iy+ym-1))*(bix>ix)*(bix<(ix+xm-1)))[0]
            if len(fpb)!=0:
                sind_bnd0=ravel_multi_index([biy[fpb]-iy,bix[fpb]-ix], self.info.ds)
                sind_bnd=r_[sind_bnd,sind_bnd0]
                if hasattr(self,'dem'): dem_bnd=r_[dem_bnd,self.info.dem_bnd[fpb]]
                if hasattr(self,'dir'): dir_bnd=r_[dir_bnd,self.info.dir_bnd[fpb]]
        self.info.sind_bnd=sind_bnd; self.info.dem_bnd=dem_bnd; self.info.dir_bnd=dir_bnd

        sind_bnd_nodata=array([]).astype('int')
        #for sind_bnd_nodata
        if len(self.info.sind_bnd_nodata)!=0:
            biy,bix=unravel_index(self.info.sind_bnd_nodata,[ym0,xm0])
            fpb=nonzero((biy>=iy)*(biy<=(iy+ym-1))*(bix>=ix)*(bix<=(ix+xm-1)))[0]
            if len(fpb)!=0:
                sind_bnd_nodata0=ravel_multi_index([biy[fpb]-iy,bix[fpb]-ix], self.info.ds)
                sind_bnd_nodata=r_[sind_bnd_nodata,sind_bnd_nodata0]
        self.info.sind_bnd_nodata=sind_bnd_nodata

        #exclude nodata bnd
        if hasattr(self,'dem'):
            fpn=self.info.dem_bnd!=self.info.nodata
            #update bnd info if there is nodata
            self.info.sind_bnd=self.info.sind_bnd[fpn].copy()
            self.info.dem_bnd=self.info.dem_bnd[fpn].copy()
            if hasattr(self,'dir'): self.info.dir_bnd=self.info.dir_bnd[fpn].copy()


    def read_deminfo(self,name=None,subdomain_size=5e6,offset=0,depth_limit=[-1e5,1e4],header=None,sdir='.'):
        '''
        read header info, and divide the data to subdomains based on subdomain_size
        work for *asc, add other format later
        '''

        #read header
        self.info=zdata()
        if header is None:
            if name.endswith('.asc'):
                with open(name,'r') as fid:
                    xm=int(fid.readline().split()[1])
                    ym=int(fid.readline().split()[1])
                    xll=float(fid.readline().split()[1])
                    yll=float(fid.readline().split()[1])
                    dxy=float(fid.readline().split()[1])
                    nval=fid.readline().split()[1]
                    if nval.lower() in ('nan','-nan'):
                        nodata=nan
                    else:
                        nodata=float(nval)

                #change yll to upper left corner
                yll=yll+(ym-1)*dxy
                skiprows=0; skipcols=0

            elif name.endswith('.tif'):
                fid = gdal.Open(name)
                xm = fid.RasterXSize
                ym = fid.RasterYSize
                xmin, xres, _, ymax, _, yres  = fid.GetGeoTransform()
                xmax = xmin + (xm-1) * xres
                #ymin = ymax + (ym-1) * yres
                xll = xmin
                yll = ymax
                dxy = xres
                #get nodata
                band = fid.GetRasterBand(1)
                nodata = band.GetNoDataValue()

        else:
           xll,xllm,yllm,yll,dxy,xm,ym,nodata=header.diminfo
           skiprows=header.skiprows; skipcols=header.skipcols
           name=header.name
           depth_limit=header.depth_limit

        #process info
        nsize=ym*xm;
        self.info.name=name
        self.info.header=[ym,xm,yll,xll,dxy]
        self.info.extent=[xll,xll+(xm-1)*dxy,yll-(ym-1)*dxy,yll]
        self.info.nodata=nodata
        self.info.ds=[ym,xm];
        self.info.skiprows=int(6+skiprows)
        self.info.usecols=arange(xm).astype('int')+int(skipcols)
        self.info.max_rows=int(ym)
        self.info.depth_limit=depth_limit

        #divide the domain into subdomains
        nsub=max(around(nsize/subdomain_size),1);
        nx=max(int(round(sqrt(nsub))),1); ny=max(int(nsub/nx),1);
        dy0=int(floor(ym/ny)); dx0=int(floor(xm/nx))

        self.info.is_subdomain=False
        self.info.nsubdomain=nx*ny;
        self.info.subdomain_shape=[ny,nx]
        self.info.subdomain_size=[dy0,dx0]
        self.domains=[]

        #boundary
        biy=r_[zeros(xm),ones(xm)*(ym-1),arange(1,ym-1),arange(1,ym-1)].astype('int')
        bix=r_[arange(xm),arange(xm),zeros(ym-2),ones(ym-2)*(xm-1)].astype('int')
        self.info.sind_bnd=ravel_multi_index([biy,bix],self.info.ds)
        self.info.sind_bnd_nodata=array([]).astype('int')

        #calcuate subdomain info
        for i in arange(ny):
            #subdomain index
            if ny==1:
                iy=0; dy=dy0; eiy=0; edy=dy0 #eiy,edy is local index and shape
            else:
                if i==0:
                    iy=0;  dy=dy0+offset;  eiy=0; edy=dy0
                elif i==(ny-1):
                    iy=i*dy0-offset;  dy=ym-iy;  eiy=offset; edy=ym-i*dy0
                else:
                    iy=i*dy0-offset;  dy=dy0+2*offset;   eiy=offset; edy=dy0

            for k in arange(nx):
                #subdomain index
                if nx==1:
                    ix=0;  dx=dx0;  eix=0; edx=dx0 #eiy,edy is local index and shape
                else:
                    if k==0:
                        ix=0;  dx=dx0+offset;  eix=0; edx=dx0
                    elif k==(nx-1):
                        ix=k*dx0-offset; dx=xm-ix;  eix=offset; edx=xm-k*dx0
                    else:
                        ix=k*dx0-offset; dx=dx0+2*offset; eix=offset; edx=dx0

                #save subdomain info
                sinfo=zdata();
                sinfo.name=name
                sinfo.header=[dy,dx,yll-iy*dxy,xll+ix*dxy,dxy]
                sinfo.extent=[xll+ix*dxy,xll+ix*dxy+(dx-1)*dxy,yll-iy*dxy-(dy-1)*dxy,yll-iy*dxy]
                sinfo.nodata=nodata
                sinfo.ds=[dy,dx]
                sinfo.skiprows=int(6+iy+skiprows)
                sinfo.usecols=arange(ix,ix+dx).astype('int')+int(skipcols)
                sinfo.max_rows=int(dy)
                sinfo.depth_limit=depth_limit

                #go global
                sinfo.is_subdomain=True
                sinfo.ixy_local=[eiy,eix]
                sinfo.ds_local=[edy,edx]
                sinfo.ixy_global=[iy,ix]
                sinfo.ds_global=[ym,xm]

                #local boundary index
                biy=r_[zeros(dx),ones(dx)*(dy-1),arange(1,dy-1),arange(1,dy-1)].astype('int')
                bix=r_[arange(dx),arange(dx),zeros(dy-2),ones(dy-2)*(dx-1)].astype('int')
                sinfo.sind_bnd=ravel_multi_index([biy,bix],sinfo.ds)
                sinfo.sind_bnd_nodata=array([]).astype('int')

                #save info
                sdata=dem(); sdata.info=sinfo
                self.domains.append(sdata)

    def read_files_info(self,names,ids,sname='deminfo',sdir='.',depth_limit=[-1e4,1e4],plot_domain=False):
        #read information about multiple dem files, and remove the overlap zone

        #read global information of files
        diminfo0=[]; nfile=len(names)
        for name in names:
            #read header
            if name.endswith('.asc'):
                with open(name,'r') as fid:
                    xm=int(fid.readline().split()[1])
                    ym=int(fid.readline().split()[1])
                    xll=float(fid.readline().split()[1])
                    yll=float(fid.readline().split()[1])
                    dxy=float(fid.readline().split()[1])
                    nval=fid.readline().split()[1]
                    if nval.lower() in ('nan','-nan'):
                        nodata=nan
                    else:
                        nodata=float(nval)
                #change yll to upper left corner
                yll=yll+(ym-1)*dxy
            elif name.endswith('.npz'):
                sinfo=loadz(name,svars=['info',]).info
                ym,xm=sinfo.ds; yll=sinfo.yll; xll=sinfo.xll; dxy=sinfo.dxy; nodata=sinfo.nodata

            elif name.endswith('.tif') or name.endswith('.TIF'):
                fid = gdal.Open(name)
                xm = fid.RasterXSize
                ym = fid.RasterYSize
                xmin, xres, _, ymax, _, yres  = fid.GetGeoTransform()
                xmax = xmin + (xm-1) * xres
                ymin = ymax + (ym-1) * yres
                xll = xmin
                yll = ymax
                dxy = xres
                #get nodata
                band = fid.GetRasterBand(1)
                nodata = band.GetNoDataValue()

            #process info
            diminfo0.append([xll,xll+(xm-1)*dxy,yll-(ym-1)*dxy,yll,dxy,xm,ym,nodata])

        if nfile==1:
            #if there is only one dem, return
            self.headers=[]
            header=zdata()
            header.diminfo=diminfo0[0]
            header.name=names[0]
            header.sname0=sname
            if ids is None:
               header.id=None
               header.sname=sname
            else:
               header.id=ids[0]
               header.sname='{}/{}_{}'.format(sdir,sname,ids[0])
            header.depth_limit=depth_limit
            header.nbs=[]
            header.names=names
            header.ids=ids
            header.skipcols=0
            header.skiprows=0
            self.headers.append(header)
        else:
            #for multiple files

            #find neighboring domains
            nbs=[]; slims=[]
            for i in arange(nfile):
                x1,x2,y1,y2,dxy=diminfo0[i][:5]

                #neighboring domain
                nb=[]
                cdxy=1.5*dxy; cx1=x1-cdxy; cx2=x2+cdxy; cy1=y1-cdxy; cy2=y2+cdxy
                for m in arange(nfile):
                    if m==i: continue
                    sx1,sx2,sy1,sy2=diminfo0[m][:4]

                    fx10=(sx1>x1)*(sx1<x2);fx20=(sx2>x1)*(sx2<x2);fy10=(sy1>y1)*(sy1<y2);fy20=(sy2>y1)*(sy2<y2)
                    fx1=(sx1>cx1)*(sx1<cx2);fx2=(sx2>cx1)*(sx2<cx2);fy1=(sy1>cy1)*(sy1<cy2);fy2=(sy2>cy1)*(sy2<cy2)

                    #if no overlap, skip
                    if not (fx1|fx2)*(fy1|fy2): continue
                    if (sum([fx1*fy1,fx2*fy1,fx1*fy2,fx2*fy2])<2)*(sum([fx10*fy10,fx20*fy10,fx10*fy20,fx20*fy20])==0): continue

                    #add overlap info
                    nb.append(m)
                nbs.append(nb);

                #find which domain corner resides
                #****************************************************************************
                #this part need mannaul input, which dem file matters in the overlapping zone
                #****************************************************************************
                snb=[None,None,None,None]
                xc=[x1,x2,x2,x1]; yc=[y1,y1,y2,y2]
                for m in arange(len(nb)):
                    if int(ids[i])>int(ids[nb[m]]): continue
                    sx1,sx2,sy1,sy2=diminfo0[nb[m]][:4]
                    for k in arange(4):
                        xi=xc[k]; yi=yc[k]
                        if (xi>sx1)*(xi<sx2)*(yi>sy1)*(yi<sy2): snb[k]=nb[m]

                #calculate new xy limit. There are eight situations
                slim=[-99999,99999,-99999,99999]
                #bottom side
                if snb[0]!=None and snb[1]!=None:
                    slim[2]=max(slim[2],diminfo0[snb[0]][3],diminfo0[snb[1]][3])

                #right side
                if snb[1]!=None and snb[2]!=None:
                    slim[1]=min(slim[1],diminfo0[snb[1]][0],diminfo0[snb[2]][0])

                #upper side
                if snb[2]!=None and snb[3]!=None:
                    slim[3]=min(slim[3],diminfo0[snb[2]][2],diminfo0[snb[3]][2])

                #left side
                if snb[0]!=None and snb[3]!=None:
                    slim[0]=max(slim[0],diminfo0[snb[0]][1],diminfo0[snb[3]][1])

                #this part needed to be add manually if results are not correct
                if sum(array(snb)!=None)==1:
                   #lower left corner
                   if snb[0]!=None:
                      slim[2]=max(slim[2],diminfo0[snb[0]][3])

                   #lower right corner
                   if snb[1]!=None:
                      slim[1]=min(slim[1],diminfo0[snb[1]][0])

                   #upper right corner
                   if snb[2]!=None:
                      if ids[i] in ['08',]:
                        slim[3]=min(slim[3],diminfo0[snb[2]][2])
                      else:
                        slim[1]=min(slim[1],diminfo0[snb[2]][0])

                   #upper left corner
                   if snb[3]!=None:
                      slim[3]=min(slim[3],diminfo0[snb[3]][2])
                slims.append(slim)

            #update new xy limits
            skiprows=[]; skipcols=[]; diminfo=[]
            for i in arange(nfile):
                x1,x2,y1,y2,dxy,xm,ym,nodata=diminfo0[i][:8]
                xi=x1+dxy*arange(xm); xind=nonzero((xi>slims[i][0])*(xi<slims[i][1]))[0]
                yi=y2-dxy*arange(ym); yind=nonzero((yi>slims[i][2])*(yi<slims[i][3]))[0]
                ix1=xind.min(); ix2=xind.max(); iy1=yind.min(); iy2=yind.max();
                ym=(iy2-iy1+1); xm=(ix2-ix1+1); y2=y2-iy1*dxy; x1=x1+ix1*dxy
                diminfo.append([x1,x1+(xm-1)*dxy,y2-(ym-1)*dxy,y2,dxy,xm,ym,nodata])
                skiprows.append(iy1); skipcols.append(ix1)

            #reorgnize info
            self.headers=[]
            for i in arange(nfile):
                header=zdata()
                header.diminfo=diminfo[i]
                header.id=ids[i]
                header.name=names[i]
                header.sname0=sname
                header.sname='{}/{}_{}'.format(sdir,sname,ids[i])
                header.depth_limit=depth_limit
                header.names=names
                header.ids=ids
                header.nbs=nbs[i]
                header.skipcols=skipcols[i]
                header.skiprows=skiprows[i]
                self.headers.append(header)

            # #plot domains before and after removing overlapping zone
            # if plot_domain:
            #     colors=['r','g','b','k','m']
            #     figure(figsize=[12,10])
            #     subplot(2,1,1)
            #     for i in arange(nfile):
            #         x1,x2,y1,y2=diminfo0[i][:4]
            #         xi=array([x1,x2,x2,x1,x1]); yi=array([y1,y1,y2,y2,y1])
            #         plot(xi,yi,'-',color=colors[mod(i,5)],alpha=0.6)
            #         text((x1+x2)/2,(y1+y2)/2,ids[i],color=colors[mod(i,5)],alpha=0.6)
            #     title('domains of each dem files: original')

            #     subplot(2,1,2)
            #     for i in arange(nfile):
            #         x1,x2,y1,y2=diminfo[i][:4]
            #         xi=array([x1,x2,x2,x1,x1]); yi=array([y1,y1,y2,y2,y1])
            #         xi=array([x1,x2,x2,x1,x1]); yi=array([y1,y1,y2,y2,y1])
            #         plot(xi,yi,'-',color=colors[mod(i,5)],alpha=0.6)
            #         text((x1+x2)/2,(y1+y2)/2,ids[i],color=colors[mod(i,5)],alpha=0.6)
            #     title('domains of each dem files: removing overlapping zone')
            #     gcf().tight_layout()
            #     savefig('{}/dem_domains_overlap'.format(sdir))
            #     close()

        #save domain information
        if len(names)!=1:
            self.save_data('{}/{}'.format(sdir,sname),'headers')

    def compute_river(self,seg=None,sind=None,acc_limit=1e4,area_limit=None,rat_prj=1.11e5,nodata=None,apply_mask=False,msg=False):
        '''
        compute river network for watersheds
        seg is segment number, sind is catchment indices.
        area_limit: compute river for cell with watershed area>area_limit
        rat_prj: conversion ratio between projection to meter (default is for lat_long)
        '''

        #compute acc_limit
        if area_limit is not None:
           dxy=self.info.header[4];
           acc_limit=area_limit/(dxy*rat_prj)**2

        #compute river mouths
        if sind is not None:
            sind0=sind.copy(); seg0=self.seg.ravel()[sind0]

        if seg is not None: #may rewrite sind
            if hasattr(seg,'__len__'):
                seg=array(seg)
            else:
                seg=array([seg,])

            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]
            seg_in,ia,ib=intersect1d(seg,seg0,return_indices=True)
            seg0=seg0[ib]; sind0=sind0[ib];

        if (sind is None) and (seg is None):
            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]

        #pre-define varibles, will update in the loop
        if nodata is None: nodata=self.info.nodata
        sind=sind0; slen=len(sind);
        num=arange(slen).astype('int'); pind0=None;

        print('computing river network')

        #each pts cooresponding to a river network
        self.rivers=[[] for i in arange(len(sind0))];

        while len(sind)!=0:
            #search the largest stream
            sind_list=self.search_upstream(sind,ireturn=7,acc_limit=acc_limit,msg=msg)

            #add rivers
            for i in arange(slen):
                if pind0 is not None: self.rivers[num[i]].append(pind0[i])
                self.rivers[num[i]].extend(sind_list[i])
                self.rivers[num[i]].append(self.info.nodata)

            #create new pind for search
            pind=[]; pnum=[];
            for i in arange(slen):
                pind.extend(sind_list[i])
                pnum.extend(ones(len(sind_list[i]))*num[i])
            pind=array(pind); pnum=array(pnum).astype('int')

            #create new sind
            sind,pnum_ind=self.search_upstream(pind,ireturn=5,acc_limit=acc_limit,msg=msg)
            slen=len(sind); num=pnum[pnum_ind];  pind0=pind[pnum_ind]

        #format
        inum=[]
        for i in arange(len(sind0)):
            river=array(self.rivers[i]).astype('int')
            #apply mask
            if apply_mask:
                fpm=self.mask.ravel()[river]
                river[~fpm]=self.info.nodata

                #only one nodata in between
                nind=nonzero(river==self.info.nodata)[0]; dind=diff(nind);
                fpd=nonzero(dind==1)[0]+1; nind=nind[fpd]

                fpn=ones(len(river)).astype('bool'); fpn[nind]=False;
                river=river[fpn]

            self.rivers[i]=river
            if len(self.rivers[i])>3: inum.append(i)
        self.rivers=array(self.rivers, dtype='object'); inum=array(inum)

        if len(self.rivers)==0: return
        if len(inum)==0: self.rivers=[]; return
        #exclude rivers with len<3
        self.rivers=self.rivers[inum]

        #add river information to info
        self.info.rseg=seg0[inum]
        self.info.rind=sind0[inum]

    def compute_boundary(self,seg=None,sind=None,acc_limit=None,level_max=100,msg=True):
        '''
        comtpue boundary for pts(sind) or segments(seg); seg overite sind
        '''

        #pre-define variables
        ds=self.info.ds; ym,xm=ds

        #compute river mouths
        if sind is not None:
            sind0=sind.copy(); seg0=self.seg.ravel()[sind0]

        if seg is not None: #may rewrite sind0
            if hasattr(seg,'__len__'):
                seg=array(seg)
            else:
                seg=array([seg,])

            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]
            seg_in,ia,ib=intersect1d(seg,seg0,return_indices=True)
            seg0=seg0[ib]; sind0=sind0[ib];

        if (sind is None) and (seg is None):
            sind0=nonzero(self.dir.ravel()==0)[0]; seg0=self.seg.ravel()[sind0]

        #compute watershed bnd
        sind=sind0.copy()

        #exclude pts with acc<acc_limit
        if acc_limit is not None:
            fpa=self.acc.ravel()[sind]>=acc_limit
            sind=sind[fpa]; seg=seg0[fpa]

        #---------------method 2-----------------------------------------------
        #find boundary pts
        fps=arange(len(sind)).astype('int'); seg=self.seg.ravel()[sind]
        while len(fps)!=0:
            #exclude pts already on domain bounday
            iy,ix=unravel_index(sind[fps],ds)
            fpb=(iy==(ym-1))|(ix==(xm-1))|(iy==0)|(ix==0)
            fps=fps[~fpb]

            #check the other three directions
            sind_next=r_[sind[fps]-xm,sind[fps]-xm-1,sind[fps]-1]
            seg_next=self.seg.ravel()[sind_next]

            #has next cells with same seg
            fp=seg_next==tile(seg[fps],3); sind_next[~fp]=sys.maxsize
            fpt=(fp).reshape([3,len(fps)]).sum(axis=0)!=0

            #reset
            sind[fps[fpt]]=sind_next.reshape([3,len(fps)]).min(axis=0)[fpt]
            fps=fps[fpt]

        #search boundary
        self.boundary=self.search_boundary(sind,level_max=level_max,msg=msg)

        #add boundary info
        #self.info.bseg=seg

    def compute_watershed(self,msg=False,ireturn=0):
        '''
        acc_max: compute acc and segment number for all watersheds
        '''

        if not hasattr(self,'dir'): sys.exit('dir not exist')

        #all the catchments()
        sind0=nonzero(self.dir.ravel()==0)[0];
        if len(sind0)==0: return

        print('---------compute watershed------------------------------------')
        #initialize acc
        self.search_upstream(sind0,ireturn=3,level_max=100,msg=msg)

        if ireturn==0:
            print('---------sort watershed number--------------------------------')
            #reorder segment number
            acc0=self.acc.ravel()[sind0] 
            ind=flipud(argsort(acc0)) #get index of acc from high to low
            sind=sind0[ind]
            seg=arange(len(sind)).astype('int')+1
            self.search_upstream(sind,ireturn=3,seg=seg,level_max=100,msg=msg)

            #add external acc (for subdomain)
            if hasattr(self.info,'sind_ext'):
                if len(self.info.sind_ext)!=0:
                    slen=len(self.info.sind_ext)
                    sind_list=self.search_downstream(self.info.sind_ext,ireturn=2)
                    sind_all=[]; acc_all=[]
                    for i in arange(slen):
                        sind_all.extend(sind_list[i])
                        acc_all.extend(ones(len(sind_list[i]))*self.info.acc_ext[i])
                    sind_all=array(sind_all); acc_all=array(acc_all)

                    #organize acc.
                    sind_unique=unique(sind_all); acc_unique=zeros(len(sind_unique)).astype('int')
                    for i in arange(slen):
                        sindc,iA,iB=intersect1d(sind_unique,sind_list[i],return_indices=True)
                        acc_unique[iA]=acc_unique[iA]+self.info.acc_ext[i]
                    self.acc.ravel()[sind_unique]=self.acc.ravel()[sind_unique]+acc_unique
        else:
            #save boundary acc
            sind0=sort(sind0)
            if len(setdiff1d(sind0,self.info.sind_bnd))!=0: 
                sys.exit('sind0-sind_bnd!=0, there still have flats/depressions in interior area.')
            self.info.sind0=sind0; self.info.acc0=self.acc.ravel()[sind0]
            #sind0 is a subset of self.info.sind_bnd and must have unique values.
            sindc,iA,iB=intersect1d(sind0,self.info.sind_bnd,return_indices=True)
            if not array_equal(sindc,sind0): 
                sys.exit('sind0!=sind_unique, sind0 is not unique')
            self.info.dem0=self.info.dem_bnd[iB].copy()

            #save seg depth
            sind=self.search_downstream(self.info.sind_bnd,ireturn=1,msg=msg)
            sind_unique,fpu=unique(sind,return_inverse=True)
            sindc,iA,iB=intersect1d(sind_unique,self.info.sind_bnd,return_indices=True)
            if not array_equal(sindc,sind_unique): sys.exit('sindc!=sind_unique')
            self.info.dem_bnd0=self.info.dem_bnd[iB][fpu].copy()

        return

    def compute_dir(self,data='dem',outname='dir',subdomain_size=1e5,zlimit=0,method=0):
        '''
        when diff(dem)<zlimit, regard them the same elevation
        method=0: don't assign dir to flat cells;
        method=1: assign dir to flat cells
        method=2: don't compute dir on boundary
        returns: dir == 0 (flats or depressions)
                     == -1 (nodata)
                     == [1, 2, 4, 8, 16, 32, 64, 128] flow direction
        '''

        #dem data
        if isinstance(data,str): #assume data is dem
            #change dem dimension
            self.dem=self.dem.ravel()
            dem0=self.dem
        else:
            dem0=data.ravel()

        #pre-calculation
        ds=self.info.ds; ym,xm=ds; nsize=prod(ds); nlim=subdomain_size
        #dir=zeros(nsize).astype('int32'); nodata=self.info.nodata
        dir=zeros(nsize, dtype=int32); nodata=self.info.nodata
        nsubdomain=int(max(round(nsize/nlim),1))

        offsets_all=array([1,-xm,-1,xm, -xm+1,-xm-1,xm-1,xm+1])
        ndirs_all=array([1, 64, 16, 4, 128, 32, 8,  2])
        #pdirs_all=array([16, 4,  1, 64, 8,   2, 128,32])

        if method==2:
            nloop=1
        else:
            nloop=3
        #calculate dir for each subdomain
        for it in arange(nloop):
            if it==0:
                nsub=nsubdomain #for subdomain
            else:
                nsub=4 #for 4 sides and 4 corners

            #for each section
            for i in arange(nsub):
                if it==0:
                    #get index for subdomain
                    if i==(nsub-1):
                        sind0=arange(i*nlim,nsize).astype('int')
                    else:
                        sind0=arange(i*nlim,(i+1)*nlim).astype('int')

                    #exclude pts on boundary
                    iy,ix=unravel_index(sind0,ds)
                    fp=(ix>0)*(ix<(xm-1))*(iy>0)*(iy<(ym-1)); sind0=sind0[fp]

                    flag=arange(8)
                elif it==1:
                    #get index for each side
                    if i==0:
                        sind0=arange(1,xm-1);
                        flag=array([0,2,3,6,7 ])
                    elif i==1:
                        sind0=arange(1,xm-1)+(ym-1)*xm;
                        flag=array([0,1,2,4,5])
                    elif i==2:
                        sind0=arange(1,ym-1)*xm
                        flag=array([0,1,3,4,7])
                    elif i==3:
                        sind0=arange(1,ym-1)*xm+(xm-1)
                        flag=array([1,2,3,5,6])
                    dir[sind0]=0
                elif it==2:
                    #get index for each corner
                    if i==0:
                        sind0=array([0])
                        flag=array([0,3,7])
                    elif i==1:
                        sind0=array([xm-1])
                        flag=array([2,3,6])
                    elif i==2:
                        sind0=array([(ym-1)*xm])
                        flag=array([0,1,4])
                    elif i==3:
                        sind0=array([nsize-1])
                        flag=array([1,2,5])
                    dir[sind0]=0

                #begin compute dir-----------------------------------------------
                #define offsets for each sid-e
                offsets=offsets_all[flag]
                ndirs=ndirs_all[flag]
                #pdirs=pdirs_all[flag]

                #exclude nodata
                sdem=dem0[sind0];
                fp=sdem!=nodata; dir[sind0[~fp]]=-1
                sind=sind0[fp]; dem=sdem[fp];  slen=sum(fp)
                if len(sind)==0: continue

                #find the neighbor with min. depth
                #sindi - index of the neighbors
                #ndem - dem of the neighbor with min. depth
                #nind - index of the neighbor with min. depth
                #ndir - direction of the neighbor with min.depth
                ndem=ones(slen)*nodata; nind=zeros(slen).astype('int');  ndir=nind.copy() # pdir=nind.copy()
                for m in arange(len(offsets)):
                    sindi=sind+offsets[m];
                    fpo=sindi>=nsize; sindi[fpo]=nsize-1
                    demi=dem0[sindi]

                    #check values
                    fpm1=(demi!=nodata)*(ndem==nodata)
                    fpm2=(demi!=nodata)*(ndem!=nodata)*(demi<ndem)
                    fpm=fpm1|fpm2

                    #assign new ndem if ndem==nodata and demi!=nodata, and if demi<ndem
                    nind[fpm]=sindi[fpm]
                    ndem[fpm]=demi[fpm]
                    ndir[fpm]=ndirs[m]
                    #pdir[fpm]=pdirs[m]

                #assign dir when dem>ndem (focal cell toward neighbor)
                fpm=(ndem!=nodata)*(dem>(ndem+zlimit))
                dir[sind[fpm]]=ndir[fpm]

                #when dem=ndem (this is a flat), assign dir first if method==1
                if method==1:
                    fpm=(ndem!=nodata)*(abs(dem-ndem)<=zlimit)*(sind<nind)
                    dir[sind[fpm]]=ndir[fpm]

        #make sure dir=-1 for nodata
        dir[dem0==nodata]=-1

        #reshape
        if isinstance(data,str): self.dem=self.dem.reshape(ds) #assume data is dem
        #dir=dir.reshape(ds).astype('int32')
        dir=dir.reshape(ds)
        exec('self.{}=dir'.format(outname))

        #add bnd info
        self.info.dir_bnd=dir.ravel()[self.info.sind_bnd].copy()

    def change_bnd_dir(self):
        '''
        let dir=0 if dir on the boundary is outflow
        '''

        ds=self.info.ds; ym,xm=ds

        #bnd info
        sind_bnd=self.info.sind_bnd
        dir_bnd=self.dir.ravel()[sind_bnd]
        iy,ix=unravel_index(sind_bnd,ds)

        #dirs on four sides
        dirs=[[32,64,128],[2,4,8],[8,16,32],[1,2,128]]
        #inds=[arange(1,xm-1),xm+arange(1,xm-1),2*xm+arange(0,ym-2),2*xm+(ym-2)+arange(0,ym-2)]
        for i in arange(4):
            if i==0:
                fpr=(iy==0)*(ix>0)*(ix<(xm-1))
            elif i==1:
                fpr=(iy==(ym-1))*(ix>0)*(ix<(xm-1))
            elif i==2:
                fpr=(ix==0)*(iy>0)*(iy<(ym-1))
            elif i==3:
                fpr=(ix==(xm-1))*(iy>0)*(iy<(ym-1))

            fpr=nonzero(fpr)[0]; diri=dir_bnd[fpr]
            fpd=nonzero((diri==dirs[i][0])|(diri==dirs[i][1])|(diri==dirs[i][2]))[0]
            self.dir.ravel()[sind_bnd[fpr[fpd]]]=0

        #dirs on four corners
        dirs=[[1,2,4],[4,8,16],[1,128,64],[16,32,64]]
        #inds=[0,xm-1,xm,2*xm-1]
        for i in arange(4):
            if i==0:
                fpr=(iy==0)*(ix==0)
            elif i==1:
                fpr=(iy==0)*(ix==(xm-1))
            elif i==2:
                fpr=(iy==(ym-1))*(ix==0)
            elif i==3:
                fpr=(iy==(ym-1))*(ix==(xm-1))
            fpr=nonzero(fpr)[0]; diri=dir_bnd[fpr]
            if len(diri)==0: continue
            if (diri!=dirs[i][0])*(diri!=dirs[i][1])*(diri!=dirs[i][2]):
                self.dir.ravel()[sind_bnd[fpr]]=0

    def fill_depression(self,level_max=100,method=0,msg=True):
        '''
        method=0: resolve all flats
        method=1: don't resolve flats on boundary (dir=0)
        method=2: resolve data when only self.info.dem_bnd_local is available, and no self.dem
        '''

        #pre-define variables
        ds=self.info.ds; ym,xm=ds

        #initialize
        #get all cells within the flat(dir==0)
        sind0=nonzero(self.dir.ravel()==0)[0]; iy,ix=unravel_index(sind0,ds)
        fp=(ix>0)*(ix<(xm-1))*(iy>0)*(iy<(ym-1))
        sind0=setdiff1d(sind0[fp],self.info.sind_bnd_nodata) #excluding cells on the boundary
        if sind0.size==0: return

        #resolve flats (after this process, there are still cells with dir == 0, these are depressions)
        if method==0:
            print('---------resolve DEM flats------------------------------------')
            # This search_upsteram returns if flat cells have incoming flow (inum != 0) or not (inum == 0)
            isum=self.search_upstream(sind0,ireturn=2,msg=msg)
            self.resolve_flat(sind0[isum==0])

        #resolve flats for single cells
        # sind0=nonzero(self.dir.ravel()==0)[0]
        # if method==1:
        #     iy,ix=unravel_index(sind0,ds)
        #     fp=(ix>0)*(ix<(xm-1))*(iy>0)*(iy<(ym-1)); sind0=sind0[fp]
        # isum=self.search_upstream(sind0,ireturn=2)
        # self.search_flat(sind0[isum==0],ireturn=6)

        #identify depression catchment
        sind0=nonzero(self.dir.ravel()==0)[0]; iy,ix=unravel_index(sind0,ds)
        fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1)); sind0=setdiff1d(sind0[fp],self.info.sind_bnd_nodata)
        if sind0.size==0: return

        #breakpoint()
        #search and mark each depression
        print('---------identify and mark depressions------------------------')
        slen=len(sind0); seg0=arange(slen)+1; h0=self.dem.ravel()[sind0];
        self.search_upstream(sind0,ireturn=3,seg=seg0,level_max=level_max,msg=msg)
        #breakpoint()
        #At this point, self.seg is created

        #get indices for each depression; here seg is numbering, not segment number
        print('---------save all depression points---------------------------')
        sind_segs=self.search_upstream(sind0,ireturn=9,seg=arange(slen),acc_calc=True,level_max=level_max,msg=msg)
        ns0=array([len(i) for i in sind_segs])

        #get depression boundary
        print('---------compute depression boundary--------------------------')
        self.compute_boundary(sind=sind0,msg=msg)

        #loop to reverse dir along streams
        print('---------fill depressions-------------------------------------')
        ids=arange(slen); iflag=0
        while len(sind0)!=0:
            iflag=iflag+1
            print('fill depression: iloop={}, ndep={}'.format(iflag,slen))
            # slen=len(ids); h0=h00[ids]; sind0=sind00[ids]; ns0=ns00[ids]
            # ids=arange(slen)
            #-------------------------------------------------------------------------
            #seg0,     ns0,     h0,      sind0,       sind_bnd,   h_bnd
            #seg_min,  ns_min,  h0_min,  sind0_min,   sind_min,   h_min  h_bnd_min  dir_min
            #-------------------------------------------------------------------------

            #exclude nonboundary indices
            sind_bnd_all=[]; ids_bnd=[]; id1=0; id2=0;
            for i in arange(slen):
                sindi=unique(self.boundary[i]); id2=id1+len(sindi)
                sind_bnd_all.extend(sindi)
                ids_bnd.append(arange(id1,id2).astype('int')); id1=id1+len(sindi)
            sind_bnd_all=array(sind_bnd_all);ids_bnd=array(ids_bnd)
            flag_bnd_all=self.search_flat(sind_bnd_all,ireturn=10,msg=msg)
            for i in arange(slen):
                self.boundary[i]=sind_bnd_all[ids_bnd[i][flag_bnd_all[ids_bnd[i]]]]

            #recalcuate bnd for seg with false boundary
            len_seg=array([len(i) for i in sind_segs]);
            len_bnd=array([len(i) for i in self.boundary])
            idz=nonzero((len_bnd==0)*(len_seg!=0))[0]
            if len(idz)!=0:
                #save boundary first
                boundary=self.boundary.copy(); delattr(self,'boundary')
                self.compute_boundary(sind=sind0[idz],msg=msg);
                if len(idz)==1:
                    boundary[idz[0]]=self.boundary[0]
                else:
                    boundary[idz]=self.boundary

                self.boundary=boundary; boundary=None

            #rearange the boundary index
            sind_bnd_all=[]; ids_bnd=[]; id1=0; id2=0;
            for i in arange(slen):
                sindi=self.boundary[i]; id2=id1+len(sindi)
                sind_bnd_all.extend(sindi)
                ids_bnd.append(arange(id1,id2).astype('int')); id1=id1+len(sindi)
            sind_bnd_all=array(sind_bnd_all);ids_bnd=array(ids_bnd);

            #find all the neighboring indices with minimum depth
            sind_min_all,h_min_all,dir_min_all=self.search_flat(sind_bnd_all,ireturn=8,msg=msg)

            #minimum for each seg
            h_min=array([h_min_all[id].min() for id in ids_bnd])
            mind=array([nonzero(h_min_all[id]==hi)[0][0] for id,hi in zip(ids_bnd,h_min)])
            sind_bnd=array([sind_bnd_all[id][mi] for id,mi in zip(ids_bnd,mind)]); h_bnd=self.dem.ravel()[sind_bnd]
            sind_min=array([sind_min_all[id][mi] for id,mi in zip(ids_bnd,mind)])
            dir_min=array([dir_min_all[id][mi] for id,mi in zip(ids_bnd,mind)])

            #get s0_min,h0_min,sind0_min
            seg_min=self.seg.ravel()[sind_min]; fps=nonzero(seg_min!=0)[0];

            ind_sort=argsort(seg_min[fps]);
            seg_1, ind_unique=unique(seg_min[fps[ind_sort]],return_inverse=True)
            seg_2,iA,iB=intersect1d(seg_1,seg0,return_indices=True)
            ids_min=zeros(slen).astype('int'); ids_min[fps[ind_sort]]=ids[iB][ind_unique]

            ns_min=zeros(slen).astype('int');    ns_min[fps]=ns0[ids_min[fps]]
            h0_min=-1e5*ones(slen);              h0_min[fps]=h0[ids_min[fps]]
            sind0_min=-ones(slen).astype('int'); sind0_min[fps]=sind0[ids_min[fps]]
            h_bnd_min=zeros(slen);               h_bnd_min[fps]=h_bnd[ids_min[fps]]

            #get stream from head to catchment
            fp1=seg_min==0
            fp2=(seg_min!=0)*(h0_min<h0)
            fp3=(seg_min!=0)*(h0_min==h0)*(ns_min>ns0)
            fp4=(seg_min!=0)*(h0_min==h0)*(ns_min==ns0)*(h_bnd>h_bnd_min)
            fp5=(seg_min!=0)*(h0_min==h0)*(ns_min==ns0)*(h_bnd==h_bnd_min)*(sind0>sind0_min)
            fph=nonzero(fp1|fp2|fp3|fp4|fp5)[0]
            sind_head=sind_bnd[fph]
            if len(fph)==0:
                #sind0 is on the boundary bounded by nodata
                #local
                self.info.sind_bnd=r_[self.info.sind_bnd,sind0]
                self.info.sind_bnd_nodata=r_[self.info.sind_bnd_nodata, sind0]
                if hasattr(self.info,'dir_bnd'): self.info.dir_bnd=r_[self.info.dir_bnd,self.dir.ravel()[sind0]]
                if hasattr(self.info,'dem_bnd'): self.info.dem_bnd=r_[self.info.dem_bnd,self.dem.ravel()[sind0]]
                break
            sind_streams=self.search_downstream(sind_head,ireturn=2,msg=msg)

            #get all stream index and its dir
            sind=[]; dir=[];
            for i in arange(len(fph)):
                id=fph[i];
                #S.seg.ravel()[sind_segs[id]]=seg_min[id]; seg0[id]=seg_min[id]
                sind_stream=sind_streams[i]
                dir_stream=r_[dir_min[id],self.dir.ravel()[sind_stream][:-1]]
                sind.extend(sind_stream); dir.extend(dir_stream)
            sind=array(sind); dir0=array(dir).copy(); dir=zeros(len(dir0)).astype('int')

            #reverse dir
            dir_0=[128,64,32,16,8,4,2,1]; dir_inv=[8,4,2,1,128,64,32,16]
            for i in arange(8):
                fpr=dir0==dir_0[i]; dir[fpr]=dir_inv[i]
            self.dir.ravel()[sind]=dir;

            #----build the linkage--------------------------------------------------
            s_0=seg0.copy(); d_0=seg_min.copy(); d_0[setdiff1d(arange(slen),fph)]=-1
            while True:
                fpz=nonzero((d_0!=0)*(d_0!=-1))[0]
                if len(fpz)==0: break

                #assign new value of d_0
                s1=d_0[fpz].copy();

                ind_sort=argsort(s1);
                s2,ind_unique=unique(s1[ind_sort],return_inverse=True)
                s3,iA,iB=intersect1d(s2,s_0,return_indices=True)
                d_0[fpz[ind_sort]]=d_0[iB[ind_unique]]

                #assign new value fo s_0
                s_0[fpz]=s1;

            #--------------------------
            seg_stream=seg0[fph[d_0[fph]==-1]] #seg that flows to another seg (!=0)
            seg_tmp,iA,sid=intersect1d(seg_stream,seg0,return_indices=True)

            seg_target=s_0[fph[d_0[fph]==-1]]; tid=zeros(len(seg_target)).astype('int');
            ind_sort=argsort(seg_target); seg_tmp, ind_unique=unique(seg_target[ind_sort],return_inverse=True)
            seg_tmp,iA,iB=intersect1d(seg_target,seg0,return_indices=True);
            tid[ind_sort]=iB[ind_unique]
            #----------------------------------------------------------------------

            #reset variables-----
            ids_stream=ids[fph]; ids_left=setdiff1d(ids,ids_stream)

            #collect boundary
            for si,ti in zip(sid,tid):
                self.seg.ravel()[sind_segs[si]]=seg0[ti]
                self.boundary[ti]=r_[self.boundary[ti],self.boundary[si]]
                sind_segs[ti]=r_[sind_segs[ti],sind_segs[si]]
                ns0[ti]=ns0[ti]+ns0[si]

            #set other seg number to zeros
            seg_stream2=seg0[fph[d_0[fph]==0]] #seg that flows to another seg (!=0)
            seg_tmp,iA,sid2=intersect1d(seg_stream2,seg0,return_indices=True)
            for si in sid2:
                self.seg.ravel()[sind_segs[si]]=0

            #update variables
            sind0=sind0[ids_left]; seg0=seg0[ids_left]; h0=h0[ids_left]; ns0=ns0[ids_left]
            self.boundary=self.boundary[ids_left]; sind_segs=sind_segs[ids_left]
            slen=len(sind0); ids=arange(slen)
        #clean
        delattr(self,'seg');delattr(self,'boundary')

    def fill_depression_global(self,level_max=100,msg=False):
        '''
        method=0: resolve depression without DEM data
        '''
        print('--------------------------------------------------------------')
        print('---------work on global domain depression---------------------')
        print('--------------------------------------------------------------')

        #pre-define variables
        ds=self.info.ds; ym,xm=ds; nodata=self.info.nodata

        #this part assign dir directly if the original dir of boundary cells is not zero
        if len(self.info.sind_bnd_local)==0: return
        sind0=self.info.sind_bnd_local; iy,ix=unravel_index(sind0,ds);
        fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1)); sind0=sind0[fp]
        dem0=self.info.dem_bnd_local[fp]; dir0=self.info.dir_bnd_local[fp]; #original
        dir=self.dir.ravel()[sind0] #at present

        #exclude nodata bnd pts
        if len(self.info.sind_bnd_nodata)!=0:
            sind_bnd_nodata,iA,iB=intersect1d(sind0,self.info.sind_bnd_nodata,return_indices=True)
            fp=setdiff1d(arange(len(sind0)),iA)
            sind0=sind0[fp]; dem0=dem0[fp]; dir0=dir0[fp]; dir=dir[fp]

        if len(sind0)==0: return

        print('---------assign dir directly if original dir is not zero------')
        #put flow into another seg if its original dir is not zero, it doesn't form loops
        fpz=nonzero((dir0!=0)*(dir==0))[0];
        self.dir.ravel()[sind0[fpz]]=dir0[fpz]

        #excude pts that forms loops
        sind_list,flag_loop=self.search_downstream(sind0[fpz],ireturn=3,msg=msg)
        fpl=nonzero(flag_loop==1)[0]; sindl=sind0[fpz[fpl]]; deml=dem0[fpz[fpl]]
        if len(fpl)!=0:
            sind_list=array([intersect1d(sindl,i) for i in sind_list[fpl]])
            sind_all=[]; ids=[]; id1=0; id2=0;
            for i in arange(len(fpl)):
                sindi=unique(sind_list[i]); id2=id1+len(sindi)
                sind_all.extend(sindi)
                ids.append(arange(id1,id2).astype('int')); id1=id1+len(sindi)
            sind_all=array(sind_all)
            sindu,fpu=unique(sind_all,return_inverse=True); sindc,iA,iB=intersect1d(sindu,sindl,return_indices=True)
            if not array_equal(sindu,sindc): sys.exit('sindc!=sindu')
            dem_all=deml[iB][fpu]
            #assign dir=0 for the minimum depth pts
            sind_min=[]
            for idi in ids:
                sind_min.append(sind_all[idi[nonzero(dem_all[idi]==min(dem_all[idi]))[0][0]]])
            sind_min=array(sind_min)
            self.dir.ravel()[sind_min]=0

        #if dem exists, use fill_depression directly
        if hasattr(self,'dem'):
            #for bnd
            if hasattr(self.info,'sind_bnd_local'):
                self.info.sind_bnd=self.info.sind_bnd_local.copy();
                delattr(self.info,'sind_bnd_local')
            self.info.dir_bnd=self.dir.ravel()[self.info.sind_bnd]
            self.info.dem_bnd=self.dem.ravel()[self.info.sind_bnd]

            self.fill_depression(method=1)
            return

        # # S.save_data('S8_m2',['dir','info'])
        # # return

        #----------------------------------------------------------------------
        #fill depression in global domain without dem (only dem on boundary and sind0 are known)
        #----------------------------------------------------------------------

        #compute all segment number
        fp=self.dir.ravel()[self.info.sind_bnd_local]==0
        sind00=sort(self.info.sind_bnd_local[fp]); seg00=arange(len(sind00)).astype('int')+1
        self.search_upstream(sind00,ireturn=3,seg=seg00,level_max=level_max,msg=msg)
        #get h00
        tmp,iA,iB=intersect1d(sind00,self.info.sind_bnd_local,return_indices=True); h00=self.info.dem_bnd_local[iB]

        #identify all depressions
        iy,ix=unravel_index(self.info.sind_bnd_local,ds)
        fp=(iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1))*(self.dir.ravel()[self.info.sind_bnd_local]==0)
        sind0=self.info.sind_bnd_local[fp]; sind0=setdiff1d(sind0,self.info.sind_bnd_nodata); slen=len(sind0)

        #get indices for each depression
        print('---------save all depression points: global---------------------')
        sind_segs=self.search_upstream(sind0,ireturn=9,acc_calc=True,level_max=level_max,msg=msg)

        print('---------assign dir to neighboring segs if possible-----------')
        iloop=0
        while slen!=0:
            iloop=iloop+1; print("assign dir: iloop={}, ndep={}".format(iloop,slen))
            #get h0,seg0
            tmp,iA,iB=intersect1d(sind0,self.info.sind_bnd_local,return_indices=True); h0=self.info.dem_bnd_local[iB].copy()
            tmp,iA,iB=intersect1d(sind0,sind00,return_indices=True); seg0=seg00[iB].copy()

            #--------------------------------------------------------------------------
            #get sind_min, dir_min and h_min
            #--------------------------------------------------------------------------
            #index
            iy0,ix0=unravel_index(sind0,ds)
            yind=r_[iy0,  iy0-1, iy0,  iy0+1, iy0-1,iy0-1, iy0+1, iy0+1]
            xind=r_[ix0+1,ix0,   ix0-1,ix0,   ix0+1,ix0-1, ix0-1, ix0+1]
            #true neighbors
            fpt=nonzero((xind>=0)*(xind<xm)*(yind>=0)*(yind<ym))[0]
            sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds);
            fptt=self.dir.ravel()[sind_true]!=-1; fpt=fpt[fptt]; sind_true=sind_true[fptt]

            #find neighbor segment depth
            fpn=nonzero(self.seg.ravel()[sind_true]!=tile(seg0,8)[fpt])[0];
            sindu,fpu=unique(sind_true[fpn],return_inverse=True);
            seguu,fpuu=unique(self.seg.ravel()[sindu],return_inverse=True)
            tmp,iA,iB=intersect1d(seguu,seg00,return_indices=True); demn=h00[iB][fpuu][fpu]

            #find seg_min, and dir_min
            dem_all=1e6*ones([8,slen]); dem_all.ravel()[fpt[fpn]]=demn;  iy_min=argmin(dem_all,axis=0);
            sind_all=-ones([8,slen]).astype('int'); sind_all.ravel()[fpt[fpn]]=sind_true[fpn]
            dir_all=tile(array([1,64,16,4,128,32,8,2]),slen).reshape([slen,8]).T

            h_b=dem_all[iy_min,arange(slen)];
            sindb=sind_all[iy_min,arange(slen)]; #on the segment boundary
            dir_min=dir_all[iy_min,arange(slen)];
            seg_min=self.seg.ravel()[sindb]

            #find sind_min, h_min, id_min (refer to index of sind0)
            sind_min=-ones(slen).astype('int'); id_min=sind_min.copy(); h_min=int(1e6)*ones(slen)
            fpn=nonzero(seg_min!=-1)[0]; seg_minu,fpu=unique(seg_min[fpn],return_inverse=True)
            tmp,iA,iB=intersect1d(seg_minu,seg00,return_indices=True);
            sind_min[fpn]=sind00[iB][fpu].copy(); h_min[fpn]=h00[iB][fpu].copy()

            tmp,iA,iB=intersect1d(seg_minu,seg0,return_indices=True);
            id_tmp=-ones(len(seg_minu)); id_tmp[iA]=iB; id_min[fpn]=id_tmp[fpu].copy()

            #--------------------------------------------------------------------------
            #compare h0 and h_min, and assign dir
            #--------------------------------------------------------------------------
            fpm=nonzero(((h0>h_min)*(sind_min!=-1))|((h0==h_min)*(sind_min!=-1)*(sind0<sind_min)))[0]
            if len(fpm)==0: fpm=nonzero((h0==h_min)*(sind_min!=-1)*(sind0>sind_min))[0]
            if len(fpm)==0: break

            #assign dir, and break loops in case
            self.dir.ravel()[sind0[fpm]]=dir_min[fpm]
            sind_list, flag_loop=self.search_downstream(sind0[fpm],ireturn=3,msg=msg)
            fpl=fpm[nonzero(flag_loop==1)[0]]; self.dir.ravel()[sind0[fpl]]=0
            fpr=fpm[nonzero(flag_loop!=1)[0]]
            if len(fpr)==0: break

            #to find seg_min_final, for mofidy
            seg_min[setdiff1d(arange(slen),fpm)]=-1;  seg_min[fpl]=-1; seg_min_final=seg_min.copy()
            fpi=nonzero(seg_min!=-1)[0];
            while True:
                seg_min_unique,fpu=unique(seg_min_final[fpi],return_inverse=True)
                segc,iA,iB=intersect1d(seg0[fpi],seg_min_unique,return_indices=True)
                if len(segc)==0:break
                seg_min_unique[iB]=seg_min[fpi[iA]]
                seg_min_final[fpi]=seg_min_unique[fpu]

            #to find id_min_final
            id_min[setdiff1d(arange(slen),fpm)]=-1; id_min[fpl]=-1;
            fpi=nonzero(id_min!=-1)[0]; id_min_final=id_min.copy()
            #because all of fpi are linked to another seg that is not in fpi itself. Therefore, id_min and fpi has no intersection
            while True:
                id_min_unique,fpu=unique(id_min_final[fpi],return_inverse=True)
                idc,iA,iB=intersect1d(fpi,id_min_unique,return_indices=True)
                if len(idc)==0: break
                id_min_unique[iB]=id_min[idc]
                id_min_final[fpi]=id_min_unique[fpu]

            #assign new seg number
            for i in fpr:
                self.seg.ravel()[sind_segs[i]]=seg_min_final[i]
                if id_min[i]!=-1: sind_segs[id_min_final[i]]=r_[sind_segs[id_min_final[i]],sind_segs[i]]

            #update sind0, sind_segs
            ids=setdiff1d(arange(slen),fpr)
            sind0=sind0[ids]; slen=len(sind0)
            sind_segs=sind_segs[ids]

        #--------------------------------------------------------------------------
        #find the shortest route to other lower segments
        #sind0,sind_segs,slen
        #--------------------------------------------------------------------------
        print('---------fill depression in global domain---------------------')

        #compute boundary
        print('---------save all depression bounary points: global-----------')
        self.compute_boundary(sind=sind0,msg=msg)

        iloop=0
        while slen!=0:
            iloop=iloop+1; print("fill depression: iloop={}, ndep={}".format(iloop,slen))
            #get h0,seg0
            tmp,iA,iB=intersect1d(sind0,self.info.sind_bnd_local,return_indices=True); h0=self.info.dem_bnd_local[iB].copy()
            tmp,iA,iB=intersect1d(sind0,sind00,return_indices=True); seg0=seg00[iB].copy()

            #exclude nonboundary indices
            sind_bnd_all=[]; ids=[]; id1=0; id2=0;
            for i in arange(slen):
                nsind=len(self.boundary[i])
                id2=id1+nsind
                sind_bnd_all.extend(self.boundary[i])
                ids.append(arange(id1,id2).astype('int')); id1=id1+nsind
            sind_bnd_all=array(sind_bnd_all); ids=array(ids)
            flag_bnd_all=self.search_flat(sind_bnd_all,ireturn=10,msg=msg)
            for i in arange(slen):
                self.boundary[i]=sind_bnd_all[ids[i][flag_bnd_all[ids[i]]]]

            #recalcuate bnd for seg with false boundary; caused by some segs resides some other segs
            len_seg=array([len(i) for i in sind_segs]);
            len_bnd=array([len(i) for i in self.boundary])
            idz=nonzero((len_bnd==0)*(len_seg!=0))[0]
            if len(idz)!=0:
                #save boundary first
                boundary=self.boundary.copy(); delattr(self,'boundary')
                self.compute_boundary(sind=sind0[idz],msg=msg);
                if len(idz)==1:
                    boundary[idz[0]]=self.boundary[0]
                else:
                    boundary[idz]=self.boundary
                self.boundary=boundary; boundary=None

            #rearrange boundary index
            sind_bnd_all=[]; h0_all=[]; ids=[]; id1=0; id2=0;
            for i in arange(slen):
                nsind=len(self.boundary[i])
                id2=id1+nsind
                sind_bnd_all.extend(self.boundary[i])
                h0_all.extend(ones(nsind)*h0[i])
                ids.append(arange(id1,id2).astype('int')); id1=id1+nsind
            sind_bnd_all=array(sind_bnd_all); h0_all=array(h0_all); ids=array(ids)

            #length of routes
            len_stream=self.search_downstream(sind_bnd_all,ireturn=4,level_max=50,msg=msg)

            #sort by length
            for i in arange(slen):
                id=ids[i]
                ind_sort=argsort(len_stream[id])
                len_stream[id]=len_stream[id][ind_sort]
                sind_bnd_all[id]=sind_bnd_all[id][ind_sort]

            #neighboring segments
            iy0,ix0=unravel_index(sind_bnd_all,ds)
            yind=r_[iy0,  iy0-1, iy0,  iy0+1, iy0-1,iy0-1, iy0+1, iy0+1]
            xind=r_[ix0+1,ix0,   ix0-1,ix0,   ix0+1,ix0-1, ix0-1, ix0+1]
            #true neighbors
            fpt=nonzero((xind>=0)*(xind<xm)*(yind>=0)*(yind<ym))[0]
            sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds);
            fptt=nonzero((self.dir.ravel()[sind_true]!=-1)*(self.seg.ravel()[sind_true]!=self.seg.ravel()[tile(sind_bnd_all,8)[fpt]]))[0]
            fpt=fpt[fptt]; sind_true=sind_true[fptt];

            #to find h_true
            seg_true=self.seg.ravel()[sind_true]
            seg_true_unique,fpu=unique(seg_true,return_inverse=True)
            tmp,iA,iB=intersect1d(seg_true_unique,seg00,return_indices=True); h_true=h00[iB][fpu]

            #compare segment depth
            flag_true=h_true<tile(h0_all,8)[fpt]
            if sum(flag_true)==0:
                seg_true_unique,fpu=unique(seg_true,return_inverse=True)
                segc,iA,iB=intersect1d(seg_true_unique,seg0,return_indices=True)
                seg_true_unique[iA]=nodata
                flag_true=seg_true_unique[fpu]!=nodata
                #if flag_true is still all false,exclude the large seg
                if sum(flag_true)==0:
                    sid=nonzero(len_seg==max(len_seg))[0][0]; seg0_left=setdiff1d(seg0,seg0[sid])
                    seg_true_unique,fpu=unique(seg_true,return_inverse=True)
                    segc,iA,iB=intersect1d(seg_true_unique,seg0_left,return_indices=True)
                    seg_true_unique[iA]=nodata
                    flag_true=seg_true_unique[fpu]!=nodata

            #create flag
            flag_all=zeros([8,len(sind_bnd_all)])==1; flag_all.ravel()[fpt]=flag_true
            seg_all=-ones([8,len(sind_bnd_all)]).astype('int'); seg_all.ravel()[fpt]=seg_true
            #sind_all=-ones([8,len(sind_bnd_all)]).astype('int'); sind_all.ravel()[fpt]=sind_true

            flagn=zeros(len(sind_bnd_all))==1;
            dirn=-ones(len(sind_bnd_all)).astype('int');
            segn=dirn.copy()

            dir0=array([16,4,1,64,8,2,128,32]).astype('int')
            for i in arange(8):
                fp=(flag_all[i])*(~flagn)
                flagn[fp]=True
                dirn[fp]=dir0[i]
                segn[fp]=seg_all[i][fp]
                #sindn[fp]=sind_all[i][fp]

            #find the outlet for each segment
            fpo=[]; sind_min=[]; dir_min=[]; seg_min=[]
            for i in arange(slen):
                id=ids[i]
                fp=nonzero(flagn[id])[0]
                if len(fp)==0: continue

                #save information about the outlet
                fpo.append(i)
                sind_min.append(sind_bnd_all[id][fp[0]])
                dir_min.append(dirn[id][fp[0]])
                seg_min.append(segn[id][fp[0]])
            fpo=array(fpo); sind_min=array(sind_min); dir_min=array(dir_min); seg_min=array(seg_min)

            #save index and dir along the shortest route
            sind_streams=self.search_downstream(sind_min,ireturn=2,msg=msg)
            sind_all=[]; dir_all=[]
            for i in arange(len(fpo)):
                sind_stream=sind_streams[i]
                sind_all.extend(sind_stream)
                dir_all.extend(r_[dir_min[i], self.dir.ravel()[sind_stream][:-1]])
            sind_all=array(sind_all); dir_all=array(dir_all); dir_all0=dir_all.copy()

            #reverse dir along the shortest route
            dir_0=[128,64,32,16,8,4,2,1]; dir_inv=[8,4,2,1,128,64,32,16]
            for i in arange(8):
                fpr=dir_all0==dir_0[i]; dir_all[fpr]=dir_inv[i]
            self.dir.ravel()[sind_all]=dir_all;

            #to find seg_min_final and id_min_final
            seg_min_final=seg_min.copy()
            while True:
                segc,iA,iB=intersect1d(seg0[fpo],seg_min_final,return_indices=True)
                if len(segc)==0: break
                seg_min_final[iB]=seg_min[iA]

            seg_min_final_unique,fpu=unique(seg_min_final,return_inverse=True)
            segc,iA,iB=intersect1d(seg0,seg_min_final_unique,return_indices=True)
            id_min_final_unique=-ones(len(seg_min_final_unique)).astype('int')
            id_min_final_unique[iB]=iA; id_min_final=id_min_final_unique[fpu]

            #assign sind_segs, S.boundary
            for i in arange(len(fpo)):
                self.seg.ravel()[sind_segs[fpo[i]]]=seg_min_final[i]
                if id_min_final[i]!=-1:
                    sind_segs[id_min_final[i]]=r_[sind_segs[id_min_final[i]],sind_segs[fpo[i]]]
                    self.boundary[id_min_final[i]]=r_[self.boundary[id_min_final[i]],self.boundary[fpo[i]]]

            #update sind0, sind_segs, S.boundary, slen
            id_left=setdiff1d(arange(slen),fpo)
            sind0=sind0[id_left]; slen=len(sind0)
            sind_segs=sind_segs[id_left]
            self.boundary=self.boundary[id_left]

        #update boundary information
        tmp,iA,iB=intersect1d(self.info.sind_bnd_local,self.info.sind_bnd_nodata,return_indices=True)
        if not array_equal(tmp,self.info.sind_bnd_nodata): sys.exit('not equal: tmp,self.info.sind_bnd_nodata')
        flag_nodata=zeros(self.info.sind_bnd_local.size).astype('int'); flag_nodata[iA]=1

        iy,ix=unravel_index(self.info.sind_bnd_local,ds)
        fp=~((iy>0)*(iy<(ym-1))*(ix>0)*(ix<(xm-1))*(flag_nodata==0));
        self.info.sind_bnd=self.info.sind_bnd_local[fp]
        self.info.dem_bnd=self.info.dem_bnd_local[fp]
        self.info.dir_bnd=self.dir.ravel()[self.info.sind_bnd]

        #delete attributes
        delattr(self.info,'sind_bnd_local'); delattr(self.info,'dem_bnd_local'); delattr(self.info,'dir_bnd_local');
        delattr(self,'seg'); delattr(self,'boundary')

    def resolve_flat(self,sind0=None,zlimit=0):
        '''
        resolve flat based on following paper
        An Efficient Assignment of Drainage Direction Over Flat Surfaces In
        Raster Digital Elevation Models". Computers & Geosciences. doi:10.1016/j.cageo.2013.01.009"
        '''

        #pre-define variables
        ds=self.info.ds; ym,xm=ds; nsize=ym*xm

        #find indices with (dir==0)
        if sind0 is None:
            sind0=nonzero(self.dir.ravel()==0)[0]

        #exclude boundary pts
        iy,ix=unravel_index(sind0,ds)
        fp=(ix>0)*(ix<(xm-1))*(iy>0)*(iy<(ym-1)); sind0=sind0[fp]
        if sind0.size==0: return

        #return cell indices with neighbors having same elevation (isum > 0)
        isum=self.search_flat(sind0,ireturn=1,zlimit=zlimit) #indices of flat
        fp=isum>0; sind0=sind0[fp]
        if sind0.size==0: return

        #obtain flat cell indices
        sind_flat=self.search_flat(sind0,ireturn=3,zlimit=zlimit); #sind_flat=unique(sind_flat)
        dir_flat0=self.dir.ravel()[sind_flat]

        #find the low edges
        fpl=dir_flat0!=0; sind_low=sind_flat[fpl]

        #find the high edges
        fp=dir_flat0==0; sind_high=self.search_flat(sind_flat[fp],ireturn=2,zlimit=zlimit);

        #initialize new dem, and assing dem begining from high edges
        if len(sind_high)!=0:
            dem_high=self.search_flat(sind_high,ireturn=4,zlimit=zlimit)[sind_flat]
        else:
            dem_high=zeros(len(sind_flat)).astype('int')

        #initialize new dem, and assing dem begining from low edges
        if len(sind_low)!=0:
            dem_low=self.search_flat(sind_low,ireturn=5,zlimit=zlimit)[sind_flat]
        else:
            dem_low=zeros(len(sind_flat)).astype('int')

        #combine dem_high and dem_low
        dem_new=2*dem_low+dem_high

        #calcuate dir for sind_flat
        dir_flat=self.search_flat(sind_flat,ireturn=9,zlimit=zlimit,dem_flat=dem_new)
        dir_flat[fpl]=dir_flat0[fpl]; self.dir.ravel()[sind_flat]=dir_flat

    def search_flat(self,sind0,ireturn=0,zlimit=0,wlevel=0,level=0,level_max=100,dem_flat=None,method=0,msg=True):
        '''
        ireturn=0: return one-level nearby indices with same elevation
        ireturn=1: return how many neighboring indices with same elevation
        ireturn=2: return high edges of flat
        ireturn=3: return all the nearby indice with same elevation (all the cells in flat)
        ireturn=4: calculate elevation from high edges
        ireturn=5: calculate elevation from low edges
        ireturn=6: assign dir to neighboring cells that is with dir=0 & without upstream cells
        ireturn=7: return neighboring indices in seg==0 & with minimum depth.
        ireturn=8: return neighboring indices in neighboring seg & with minimum depth &dir
                   method=1 means no dem available
        ireturn=9: compute and return dir in sind_flat region based on dem_flat
        ireturn=10: return flags (True for seg-boundary cells; False for non-seg-boundary cells)
        '''

        #pre-define variables
        slen=len(sind0); ds=self.info.ds; ym,xm=ds; nsize=ym*xm
        if hasattr(self,'dem'): dem0=self.dem.ravel()[sind0]

        #convert index
        iy0,ix0=unravel_index(sind0,ds)

        #construct maxtrix for neighbor
        #yind=r_[iy0-1,iy0-1,iy0-1,iy0,iy0+1,iy0+1,iy0+1,iy0]
        #xind=r_[ix0+1,ix0,ix0-1,ix0-1,ix0-1,ix0,ix0+1,ix0+1]
        yind=r_[iy0,  iy0-1, iy0,  iy0+1, iy0-1,iy0-1, iy0+1, iy0+1]
        xind=r_[ix0+1,ix0,   ix0-1,ix0,   ix0+1,ix0-1, ix0-1, ix0+1]

        #true neighbors
        fpt=nonzero((xind>=0)*(xind<xm)*(yind>=0)*(yind<ym))[0]
        sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds);
        fptt=self.dir.ravel()[sind_true]!=-1; fpt=fpt[fptt]; sind_true=sind_true[fptt]
        if hasattr(self,'dem'): dem_true=self.dem.ravel()[sind_true]

        if ireturn==6:
            #get neighboring minimum dem and indices
            fp=dem_true==self.info.nodata;  dem_true[fp]=1e5
            dem=ones([8,slen])*1e5; dem.ravel()[fpt]=dem_true
            iy_min=argmin(dem,axis=0); dem_min=dem[iy_min,arange(slen)]

            #modify sind0's dir and dem
            #dir0=array([128,64,32,16,8,4,2,1])
            dir0=array([1,64,16,4, 128,32,8,2])
            self.dir.ravel()[sind0]=dir0[iy_min]
            self.dem.ravel()[sind0]=dem_min

        if ireturn==7:
            #get neighboring minimum dem and indices in seg=0
            fp=(self.seg.ravel()[sind_true]!=0)|(dem_true==self.info.nodata); dem_true[fp]=1e5;
            dem=ones([8,slen])*1e5; dem.ravel()[fpt]=dem_true

            sind=zeros([8,slen]).astype('int'); sind.ravel()[fpt]=sind_true
            iy_min=argmin(dem,axis=0);
            sind_min=sind[iy_min,arange(slen)]; dem_min=dem[iy_min,arange(slen)];
            return sind_min, dem_min

        if ireturn==8:
            #get neighboring minimum dem and indices in seg=0
            dir=tile(array([16,4,1,64,8,2,128,32]),slen).reshape([slen,8]).T
            seg_true=self.seg.ravel()[sind_true]; seg0=self.seg.ravel()[sind0]

            if method==0:
                fp=(seg_true==tile(seg0,8)[fpt])|(dem_true==self.info.nodata); dem_true[fp]=1e5;
                dem=ones([8,slen])*1e5; dem.ravel()[fpt]=dem_true
                sind=zeros([8,slen]).astype('int'); sind.ravel()[fpt]=sind_true

                iy_min=argmin(dem,axis=0);
                sind_min=sind[iy_min,arange(slen)];
                dem_min=dem[iy_min,arange(slen)];
                dir_min=dir[iy_min,arange(slen)];
                return sind_min,dem_min,dir_min
            elif method==1:
                #find depth for sindn based on self.info.dem_bnd_local
                fpn=nonzero(seg_true!=tile(seg0,8)[fpt])[0];
                sindn=self.search_downstream(sind_true[fpn],ireturn=1,msg=True)
                sindu, fpu=unique(sindn,return_inverse=True)
                sindc,iA,iB=intersect1d(self.info.sind_bnd_local,sindu,return_indices=True)
                if not array_equal(sindu,sindc): sys.exit('sindu!=sindc: wrong')
                dem=ones([8,slen])*1e5; dem.ravel()[fpt[fpn]]=self.info.dem_bnd_local[iA][fpu]
                sind=zeros([8,slen]).astype('int'); sind.ravel()[fpt[fpn]]=sind_true[fpn]
                sindc=zeros([8,slen]).astype('int'); sindc.ravel()[fpt[fpn]]=sindn

                iy_min=argmin(dem,axis=0)
                sind_min=sind[iy_min,arange(slen)]
                dem_min=dem[iy_min,arange(slen)]
                dir_min=dir[iy_min,arange(slen)]
                sindc_min=sindc[iy_min,arange(slen)]
                return sind_min,sindc_min,dem_min,dir_min

        if ireturn==10:
            seg_next=-ones(8*slen); seg_next[fpt]=self.seg.ravel()[sind_true]
            isum=reshape(seg_next==tile(self.seg.ravel()[sind0],8),[8,slen]).sum(axis=0)
            fps=isum!=8
            return fps

        if ireturn==2:
            fph=nonzero(dem_true>(tile(dem0,8)[fpt]+zlimit))[0]
            num=zeros(slen*8); num[fpt[fph]]=1
            isum=num.reshape([8,slen]).sum(axis=0)
            fp=isum>0; sind_next=sind0[fp]
            return sind_next

        #neighbor with same elevation (sind_next)
        fps=nonzero(abs(dem_true-tile(dem0,8)[fpt])<=zlimit)[0]
        sind_next=sind_true[fps]

        #exclude low edges
        if ireturn==4:
            fpl=nonzero(self.dir.ravel()[sind_next]==0)[0]; fpl_len=len(sind_next)
            sind_next=sind_next[fpl]

        #unique indices
        sind_next,fpu,fpu_inverse=unique(sind_next,return_index=True,return_inverse=True)

        if ireturn==9:
            dir_flat=zeros(slen).astype('int')
            dem0=int(1e5)*ones(nsize).astype('int'); dem0[sind0]=dem_flat;
            dem=int(1e5)*ones([8,slen]).astype('int'); dem.ravel()[fpt[fps]]=dem0[sind_true[fps]]
            iy_min=argmin(dem,axis=0); dem_min=dem[iy_min,arange(slen)]
            dir_min=tile(array([1,64,16,4,128,32,8,2]),[slen,1]).T[iy_min,arange(slen)]
            fp=dem_flat>dem_min; dir_flat[fp]=dir_min[fp]
            return dir_flat

        if ireturn==0:
            return sind_next

        if ireturn==1:
            num=zeros(slen*8); num[fpt[fps]]=1
            isum=reshape(num,[8,slen]).sum(axis=0).astype('int')
            return isum

        if ireturn in [3,4,5]:
            if wlevel==0: #first level loop
                #init
                self.sind_next=sind0.copy()
                self.sind_list=[]
                self.flag_search=True
                #self.dem_flag=zeros(nsize).astype('int32')
                self.dem_flag=zeros(nsize, dtype=int32)

                if ireturn in [4,5]:
                    self.v0=1
                    self.vmax=None
                    #self.dem_flat=zeros(nsize).astype('int32')
                    self.dem_flat=zeros(nsize, dtype=int32)
                    self.dem_flat_save=None
                if ireturn==4: self.sind_list.append(sind0)

                #1st round of search loop from outside to inside
                while self.flag_search:
                    sind_next=self.sind_next
                    if len(sind_next)!=0:
                        self.search_flat(sind_next,ireturn=ireturn,wlevel=1,level=0,level_max=level_max)

                #2nd round of search loop from inside to outside
                if ireturn==4:
                    self.dem_flat_save=self.dem_flat.copy()
                    for i in arange(len(self.sind_list)):
                        self.vmax=self.search_flat(self.sind_list[-i-1],ireturn=ireturn,wlevel=2,level=0,level_max=level_max)

                #save results
                if ireturn==3:
                    sind=array(self.sind_list)
                else:
                    dem_flat=self.dem_flat

                #clean
                delattr(self,'sind_next'); delattr(self,'sind_list');delattr(self,'dem_flag'); delattr(self,'flag_search')
                if ireturn in [4,5]: delattr(self,'v0'); delattr(self,'vmax');delattr(self,'dem_flat');delattr(self,'dem_flat_save');

                #return results
                if ireturn==3:
                    return sind
                else:
                    return dem_flat

            elif wlevel==1: #2nd level search, for ireturn==3
                self.dem_flag[sind0]=1
                fpn=self.dem_flag[sind_next]==0; sind_next=sind_next[fpn]
                if ireturn==3: self.sind_list.extend(sind0);
                if ireturn in [4,5]: self.dem_flat[sind0]=(ones(slen)*(level+self.v0)).astype('int32')

                if level!=level_max:
                    if sum(fpn)!=0: #continue
                        self.search_flat(sind_next,ireturn=ireturn,wlevel=1,level=level+1,level_max=level_max)
                    else:
                        self.flag_search=False #reach the end
                else:
                    if sum(fpn)!=0:
                        self.sind_next=sind_next
                        if ireturn in [4,5]:
                            self.sind_list.append(sind_next)
                            self.v0=int(self.v0+level_max+1)

            elif wlevel==2: #2nd level search, for ireturn==4
                # self.dem_flag[sind0]=0
                vmax=self.dem_flat[sind0].copy()
                #define the sind_next using self.dem_flat_save
                fpn=nonzero(tile(vmax,8)[fpt[fps[fpu]]]<self.dem_flat_save[sind_next])[0];
                fpn_len=len(sind_next); sind_next=sind_next[fpn]

                if level!=level_max:
                    if len(fpn)!=0: #continue
                        vmax_next=self.search_flat(sind_next,ireturn=ireturn,wlevel=2,level=level+1,level_max=level_max)
                    else:
                        self.flag_search=False #reach the end
                else:
                    if sum(fpn)!=0:
                        vmax_next=self.vmax

                #modify dem_flat
                if sum(fpn)!=0:
                    #replace vmax with vmax_next if vmax_next is larger
                    num_n=zeros(fpn_len); num_n[fpn]=vmax_next
                    num=zeros(8*slen);  num[fpt[fps[fpl]]]=num_n[fpu_inverse]

                    nmax=reshape(num,[8,slen]).max(axis=0)
                    fp=vmax<nmax; vmax[fp]=nmax[fp]

                self.dem_flat[sind0]=vmax-self.dem_flat[sind0]

                return vmax

    def search_boundary(self,sind0,wlevel=0,level=0,level_max=100,msg=True):
        '''
        search pts of boundary pts of watershed segment, sind0 are the initial index of segments
        '''

        #pre-defind variables
        slen=len(sind0); ds=self.info.ds
        # print(sind0)
        iy0,ix0=unravel_index(sind0,ds)
        seg0=self.seg.ravel()[sind0]

        #construct matrix
        yind=r_[iy0-1,iy0-1,iy0-1,iy0,iy0+1,iy0+1,iy0+1,iy0]
        xind=r_[ix0+1,ix0,ix0-1,ix0-1,ix0-1,ix0,ix0+1,ix0+1]
        sind=-ones(8*slen).astype('int'); seg=sind.copy()

        #true neighbor
        fpt=nonzero((xind>=0)*(xind<ds[1])*(yind>=0)*(yind<ds[0]))[0]
        sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds);
        fptt=self.dir.ravel()[sind_true]!=-1; fpt=fpt[fptt]; sind_true=sind_true[fptt]
        sind[fpt]=sind_true; seg[fpt]=self.seg.ravel()[sind_true]

        #indices belong to the same segment
        fps=seg!=tile(seg0,8); sind[fps]=-1

        #exclude previous cells
        # if hasattr(self,'pind'): fpn=sind==tile(self.pind,8); sind[fpn]=0
        if hasattr(self,'ppind'): fpn=sind==tile(self.ppind,8); sind[fpn]=-1

        sind=sind.reshape([8,slen])

        iy0=zeros(slen).astype('int');
        if hasattr(self,'pind'):
            pind=self.pind
        else:
            pind=0
        #find the starting search pt
        for i in arange(8):
            fp=sind[i,:]==pind
            iy0[fp]=i

        #find the next pts
        sind_next=-ones(slen).astype('int')
        ix=arange(slen).astype('int')
        for i in arange(8):
            iy0=iy0+1
            iy1=mod(iy0,8); iy2=mod(iy0+1,8)
            fpn=(sind[iy1,ix]==-1)*(sind[iy2,ix]!=-1)*(sind_next==-1)
            sind_next[fpn]=sind[iy2,ix][fpn]

        #if sind_next=0, search back
        fpb=sind_next==-1; sind_next[fpb]=sind0[fpb]

        #recursive search
        if wlevel==0:
            #init
            self.sind_next=sind0.copy()
            self.sind0=sind0.copy()
            self.pind=sind0.copy() #previous pts
            self.ppind=sind0.copy() #previous pts
            self.flag_search=True
            self.sind_list=[[i] for i in sind0]
            self.snum=arange(slen).astype('int')

            #1-level search loop
            iflag=0
            while self.flag_search:
                iflag=iflag+1;
                if msg: print('search boundary: loop={}, npt={}'.format(iflag,len(self.sind_next)))
                if len(self.sind_next)==0: break
                self.search_boundary(self.sind_next,wlevel=1,level=0,level_max=level_max)
                #print('search bnd: {}'.format(iflag*level_max))

            #save list info

            sind_list=array([array(i) for i in self.sind_list])
            if len(sind0)==1:
                slist_copy=array([0]).astype('O')
                slist_copy[0]=sind_list[0].copy()
                sind_list=slist_copy

            #clean
            delattr(self,'sind_next'); delattr(self,'sind0'); delattr(self,'flag_search')
            delattr(self,'sind_list'); delattr(self,'snum'); delattr(self,'pind');delattr(self,'ppind');

            return sind_list

        elif wlevel==1:
            fpz=sind_next!=self.sind0
            [self.sind_list[i].append(j) for i,j in zip(self.snum[~fpz],sind_next[~fpz])]
            if sum(fpz)==0: self.flag_search=False
            if (level!=level_max) and sum(fpz)!=0: #not reach the end
                self.sind0=self.sind0[fpz]
                self.snum=self.snum[fpz]
                self.ppind=self.pind[fpz]
                self.pind=self.sind_next[fpz]
                self.sind_next=sind_next[fpz]

                #save pts to the list
                [self.sind_list[i].append(j) for i,j in zip(self.snum,self.sind_next)]

                #continue search
                self.search_boundary(self.sind_next,wlevel=1,level=level+1,level_max=level_max)

    def search_upstream(self,sind0,ireturn=0,seg=None,acc_limit=None,wlevel=0,level=0,level_max=100,acc_calc=False,msg=True):
        '''
        ireturn=0: all the 1-level upstream index. if seg is not None, return seg number (seg_up) also
        ireturn=1: all the 1-level upstream index, also with flags of true neighbor and incoming flow
        ireturn=2: number of upstream indices
        ireturn=3: search all upstreams, compute acc and seg
        ireturn=4: just one-level upstream index with largest acc
        ireturn=5: all the one-level upstream index except the one with largest acc, and numbering index
        ireturn=6: uppermost upstream index along the largest river stem
        ireturn=7: save all the cells along the largest upstream river
        ireturn=8: return how many neighbors with same segment number after acc and seg are known
        ireturn=9: return segment indices.
        '''

        #--pre-define variables
        slen=len(sind0); ds=self.info.ds
        dir_in0=[8,4,2,1,128,64,32,16]

        #convert index
        iy0,ix0=unravel_index(sind0,ds)

        #compute init neighbor index
        yind=r_[iy0-1,iy0-1,iy0-1,iy0,iy0+1,iy0+1,iy0+1,iy0]
        xind=r_[ix0+1,ix0,ix0-1,ix0-1,ix0-1,ix0,ix0+1,ix0+1]
        if seg is not None: seg0=tile(seg,8)

        #true neighbor index (fpt - all 8 neighbors)
        fpt=nonzero(((xind>=0)*(xind<ds[1])*(yind>=0)*(yind<ds[0])).ravel())[0]

        #init dir
        dir_in=tile(dir_in0,[len(sind0),1]).T.ravel()

        #true neighbor (fpt - all neighbors not nodata)
        sind_true=ravel_multi_index([yind[fpt],xind[fpt]],ds)
        fptt=self.dir.ravel()[sind_true]!=-1
        fpt=fpt[fptt]
        sind_true=sind_true[fptt]
        if seg is not None: seg_true=seg0[fpt]

        #dir diff for true neighbors
        dir_diff=self.dir.ravel()[sind_true]-dir_in[fpt]

        #index of neighbors with incoming flow (fpf)
        fpf=nonzero(dir_diff==0)[0]
        sind_up=sind_true[fpf]
        if seg is not None: seg_up=seg_true[fpf]

        if ireturn==0:
            if seg is None:
                return sind_up
            else:
                return sind_up, seg_up

        if ireturn==1:
            return sind_up, fpt[fpf]

        if ireturn==2:
            num=zeros(slen*8); num[fpt[fpf]]=1
            #num_up > 0 means cells have incoming flow from upstream
            num_up=(reshape(num,[8,slen]).sum(axis=0)).astype('int')
            return num_up

        if ireturn==8:
            seg0=self.seg.ravel()[sind0]
            #method 1
            fps=nonzero(tile(seg0,8)[fpt]==self.seg.ravel()[sind_true])[0]
            num=zeros(slen*8); num[fpt[fps]]=1

            num_up=(reshape(num,[8,slen]).sum(axis=0)).astype('int')
            return num_up

        if ireturn==3:
            if wlevel==0: #first level recursive search
                #init
                if seg is None:
                    seg0=arange(slen)+1
                    #self.acc=zeros(prod(ds)).astype('int'); self.acc[self.dir.ravel()==-1]=-1
                    self.acc=zeros(prod(ds), dtype=int32); self.acc[self.dir.ravel()==-1]=-1
                else:
                    seg0=seg

                self.seg=zeros(prod(ds), dtype=int32); self.seg[self.dir.ravel()==-1]=-1
                #self.seg=zeros(prod(ds)).astype('int'); self.seg[self.dir.ravel()==-1]=-1
                self.sind_list=[]
                self.seg_list=[]
                self.flag_search=True

                self.sind_next=sind0.copy(); self.seg_next=seg0
                self.sind_list.append(sind0); self.seg_list.append(seg0)

                #1-level search loop, from downstream to upstream
                iflag=0
                while self.flag_search:
                    iflag=iflag+1
                    if msg: print('search upstream,ireturn={}: loop={}, npt={}'.format(ireturn,iflag,len(self.sind_next)))
                    if len(self.sind_next)==0: break
                    self.search_upstream(self.sind_next,ireturn=ireturn,seg=self.seg_next, wlevel=1,level=0,level_max=level_max)

                #search from upstream to downstream on the 1-level
                for i in arange(len(self.sind_list)):
                    if seg is not None: continue
                    if msg: print('search upstream backward,ireturn={}: loop={}, npt={}'.format(ireturn,len(self.sind_list)-i,len(self.sind_list[-i-1])))
                    self.search_upstream(self.sind_list[-i-1],ireturn=ireturn,seg=self.seg_list[-i-1], wlevel=1,level=0,level_max=level_max,acc_calc=True)

                #clean
                delattr(self,'sind_list');delattr(self,'seg_list'); delattr(self,'flag_search')
                delattr(self,'sind_next'); delattr(self,'seg_next')
                if seg is None: delattr(self,'seg')

                #reshape
                if hasattr(self,'acc'): self.acc=self.acc.reshape(ds)
                if seg is not None: self.seg=self.seg.reshape(ds)

                return

            elif wlevel==1: #2-level recursive search
                #assign seg number
                self.seg[sind0]=seg

                acc_up=0; acc=0
                if level!=level_max:
                    if len(sind_up)==0: #reach the end
                        self.flag_search=False
                    else:
                        acc_up=self.search_upstream(sind_up,ireturn=ireturn,seg=seg_up,wlevel=1,level=level+1,level_max=level_max,acc_calc=acc_calc)

                else:
                    if not acc_calc:
                        self.sind_next=sind_up
                        self.seg_next=seg_up
                        self.sind_list.append(sind_up)
                        self.seg_list.append(seg_up)
                    else:
                        acc_up=self.acc[sind_up]

                if acc_calc:
                    #compute acc for sind0
                    #get flags of sind_up first
                    sind_up,fptf=self.search_upstream(sind0,ireturn=1)
                    acc0=zeros(slen*8); acc0[fptf]=acc_up
                    acc=reshape(acc0,[8,slen]).sum(axis=0)+1

                    self.acc[sind0]=acc

                return acc

        #----------------------------------------------------------------------
        #code below works after self.acc is computed, except for ireturn=9
        #----------------------------------------------------------------------
        sind=-ones(8*slen).astype('int'); sind[fpt[fpf]]=sind_up

        if ireturn==9:
            sind_next=sind_up
        else:
            #get corresponding acc
            acc=zeros(slen*8).astype('int'); acc[fpt[fpf]]=self.acc.ravel()[sind_up];
            acc=acc.reshape([8,slen]); sind=sind.reshape([8,slen])

            #apply limit
            if acc_limit is not None:
                fpc=acc<acc_limit; acc[fpc]=0; sind[fpc]=-1;

            #get index for cells_next with maximum acc
            ix=arange(slen).astype('int'); iy=argmax(acc,axis=0)
            sind_next=sind[iy,ix]

        #just one level
        if ireturn==4:
            return sind_next

        #one-level upstream, all the cells except the largest one
        if ireturn==5:
            sind[iy,ix]=-1; fpn=sind!=-1;
            nind0=tile(arange(slen),8).reshape([8,slen]).astype('int')
            sind_next=sind[fpn]; nind=nind0[fpn]
            return sind_next,nind

        #get to the uppermost pts, using 2-level recursive search.
        #1-level recursive search will crash reaching the setrecursionlimit
        if ireturn in [6,7,9]:
            if wlevel==0: #first level recursive
                #init
                #if not hasattr(self,'sind_next'):
                self.sind_next=sind0.copy()
                self.sind_flag=ones(slen)
                self.flag_search=True
                self.sind_list=[[i] for i in sind0]
                self.sind_bnd=[[] for i in sind0]
                self.sind_seg=[[] for i in sind0]
                self.seg_next=arange(slen)

                #1-level search loop
                iflag=0
                while self.flag_search:
                    iflag=iflag+1
                    if msg: print('search upstream,ireturn={}: loop={}, npt={}'.format(ireturn,iflag,len(self.sind_next)))
                    if ireturn==9:
                        sind_next=self.sind_next
                    else:
                        fp=self.sind_flag>0; sind_next=self.sind_next[fp]
                    if len(sind_next)==0: break;
                    self.search_upstream(sind_next,ireturn=ireturn,seg=self.seg_next,acc_limit=acc_limit,wlevel=1,level=0,level_max=level_max,acc_calc=acc_calc)

                #search finished at this point
                sind_next=self.sind_next
                sind_list=array([array(i) for i in self.sind_list], dtype = 'object')
                # sind_bnd=array([array(i) for i in self.sind_bnd])
                sind_seg=array([sort(array(i)).astype('int') for i in self.sind_seg])

                #clean
                delattr(self,'sind_next'); delattr(self,'flag_search'); delattr(self,'sind_list')
                delattr(self,'sind_bnd'); delattr(self,'seg_next'); delattr(self,'sind_flag'); delattr(self,'sind_seg')

                if ireturn==6:
                    return sind_next
                elif ireturn==7:
                    return sind_list
                elif ireturn==9:
                    if acc_calc:
                        return sind_seg
                    else:
                        return sind_bnd

            elif wlevel==1: #second level recursive search
                fpnz=nonzero(sind_next!=-1)[0]; fpz=nonzero(sind_next==-1)[0]
                if level!=level_max:
                    if ireturn==9:
                        seg_next=seg_up
                        #collect boundary cells
                        # isum=sind.reshape([8,slen]).sum(axis=0); fpb=isum==-8
                        # [self.sind_bnd[i].append(j) for i,j in zip(seg[fpb],sind0[fpb])]
                        if acc_calc:
                            [self.sind_seg[i].append(j) for i,j in zip(seg,sind0)]
                    else:
                        seg_next=self.seg_next

                    if len(fpnz)==0:
                        #reach the end
                        if ireturn!=9: fpn=self.sind_flag>0; self.sind_next[fpn]=sind0
                        if ireturn==7:
                            ind_list=nonzero(fpn)[0]; sind_list=sind0
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]
                        self.flag_search=False
                    else:
                        if ireturn==9:
                            self.sind_next=sind_up
                            self.seg_next=seg_next
                        else:
                            #save the pts that reach the end
                            sind_next[fpz]=sind0[fpz]
                            fpn=self.sind_flag>0; self.sind_next[fpn]=sind_next; self.sind_flag[nonzero(fpn)[0][fpz]]=0
                        if ireturn==7:
                            ind_list=nonzero(fpn)[0]; sind_list=abs(sind_next)
                            [self.sind_list[i].append(j) for i,j in zip(ind_list,sind_list)]

                        #continue search the rest pts
                        self.search_upstream(sind_next[fpnz],ireturn=ireturn,seg=seg_next,acc_limit=acc_limit,wlevel=1,level=level+1,level_max=level_max,acc_calc=acc_calc)

    def search_downstream(self,sind0,ireturn=0,wlevel=0,level=0,level_max=500,msg=True):
        '''
        ireturn=0: return 1-level downstream index
        ireturn=1: most downstream index
        ireturn=2: save all the index along the downstream, return if reach dir=0.
                   return downstream network. This will fail for loops
        ireturn=3: save all the index along the downstream, return if reach dir=0 or identify loops
                   return downstream network and flag for loops (0: not loop; 1: 1st-type loop; 2: 2nd-type loop).
        ireturn=4: save length from upstream to downstream
        '''

        #pre-define variables
        slen=len(sind0); ds=self.info.ds

        #variables for each direction
        dir_out=array([128,64,32,16,8,4,2,1]).astype('int')
        offset_y=array([-1,-1,-1,0,1,1,1,0]).astype('int')
        offset_x=array([1,0,-1,-1,-1,0,1,1]).astype('int')

        #check each dir
        sdir=self.dir.ravel()[sind0]
        iy0,ix0=unravel_index(sind0,ds)
        yind=-ones(slen).astype('int'); xind=yind.copy(); sind_next=yind.copy()
        for i in arange(8):
            fp=sdir==dir_out[i]
            yind[fp]=iy0[fp]+offset_y[i]
            xind[fp]=ix0[fp]+offset_x[i]

        #true donwstream
        fpt=(xind>=0)*(xind<ds[1])*(yind>=0)*(yind<ds[0])
        if sum(fpt)!=0:
            sind_next[fpt]=ravel_multi_index([yind[fpt],xind[fpt]],ds)

        if ireturn==0:
            return sind_next

        if ireturn in [1,2,3,4]:
            if wlevel==0: #first level recrusive search
                #init
                #if not hasattr(self,'sind_next'):
                self.sind_next=sind0.copy()
                self.flag_search=True
                self.flag_next=ones(slen)==1
                self.sind_list=[[i] for i in sind0]
                self.len_stream=zeros(slen).astype('int')

                #1-level search loop
                iflag=0
                while self.flag_search:
                     iflag=iflag+1;
                     if msg: print('search downstream,ireturn={}: loop={}, npt={}'.format(ireturn,iflag,sum(self.flag_next)))

                     #to check whether there is a loop
                     if ireturn==3:
                         ns0=array([len(i) for i in self.sind_list])
                         ns=array([len(unique(i)) for i in self.sind_list])
                         fpe=nonzero(((ns0-ns)>1)*(self.flag_next))[0]
                         for i in arange(len(fpe)):
                              self.sind_list[fpe[i]]=self.sind_list[fpe[i]][:ns[fpe[i]]+2]
                         self.flag_next[fpe]=False

                     if sum(self.flag_next)==0: break;
                     self.search_downstream(self.sind_next[self.flag_next],ireturn=ireturn,wlevel=1,level=0,level_max=level_max)

                #search finished at this point
                sind_next=self.sind_next
                if ireturn==2:
                    sind_list=array([array(i)[:-1] for i in self.sind_list], dtype='object')
                elif ireturn==3:
                    sind_list=array([array(i)[:-1] for i in self.sind_list], dtype='object')
                    ns0=array([len(i) for i in sind_list])
                    ns=array([len(unique(i)) for i in sind_list])
                    fpl=array([i[0]==i[-1] for i in sind_list])*(ns0!=1)

                    flag_loop=zeros(slen).astype('int')
                    flag_loop[ns0!=ns]=2; flag_loop[fpl]=1
                elif ireturn==4:
                    len_stream=self.len_stream.copy()

                #clean
                delattr(self,'sind_next'); delattr(self,'flag_search'); delattr(self,'flag_next')
                delattr(self,'sind_list'); delattr(self,'len_stream');

                if ireturn==1:
                    return sind_next
                elif ireturn==2:
                    return sind_list
                elif ireturn==3:
                    return sind_list, flag_loop
                elif ireturn==4:
                    return len_stream

            elif wlevel==1: #second level recursive search
                fpz=nonzero(self.flag_next)[0]
                fpn=nonzero(sind_next==-1)[0]
                fpp=nonzero(sind_next!=-1)[0]

                if level!=level_max:
                    if len(fpp)==0:
                        #reach the end
                        self.sind_next[fpz]=sind0
                        if ireturn in [2,3]:
                            [self.sind_list[i].append(j) for i,j in zip(fpz,sind0)]
                        elif ireturn==4:
                            self.len_stream[fpz]=self.len_stream[fpz]+1

                        self.flag_search=False
                    else:
                        #save the pts that reach the end
                        sind_next[fpn]=sind0[fpn]
                        self.sind_next[fpz]=sind_next
                        self.flag_next[fpz[fpn]]=False
                        if ireturn in [2,3]:
                            [self.sind_list[i].append(j) for i,j in zip(fpz,sind_next)]
                        elif ireturn==4:
                            self.len_stream[fpz]=self.len_stream[fpz]+1

                        #continue search the rest pts
                        self.search_downstream(sind_next[fpp],ireturn=ireturn,wlevel=1,level=level+1,level_max=level_max)

    # def compute_extent(self,header=None):
    #     #recalculate extent from header
    #     if header is None: header=self.info.header
    #     ym,xm,yll,xll,dxy=header
    #     extent=array([xll,xll+(xm-1)*dxy,yll,yll+(ym-1)*dxy])
    #     return extent

    def remove_nan(self,name='dem',nodata=-99999):
        '''
        change nan to nodata
        '''

        fpn=isnan(self.dem)|(self.dem==self.info.nodata)
        self.dem[fpn]=nodata
        self.info.nodata=nodata

    def get_coordinates(self,sind0,header=None,nodata=None):
        '''
        example: sind0=nonzero(self.dir.ravel()==0)[0], or sind0=nonzero(self.dir==0)
        '''

        #check parameters first
        if header is None: header=self.info.header
        if nodata is None: nodata=self.info.nodata
        ym,xm,yll,xll,dxy=header; ds=[ym,xm]

        sind=sind0.copy().ravel();
        #check index
        if len(sind)==2 and hasattr(sind[0],'__len__'):
            indy,indx=sind  #sind0=[indy,indx]
            fpn=(indy==nodata)|(indx==nodata)
        else:
            fpn=sind0==nodata; sind[fpn]=0;
            indy,indx=unravel_index(sind,ds)
            indy[fpn]=nodata; indx[fpn]=nodata

        #convert index to coordinates
        xi=xll+indx*dxy; xi[fpn]=nodata
        yi=yll-indy*dxy; yi[fpn]=nodata

        #reshape
        xi=reshape(xi,sind0.shape)
        yi=reshape(yi,sind0.shape)

        return yi,xi

    def combine_shp(self,name,sdir='.',npt_sample=None,npt_smooth=None):

        #read fnames
        S=loadz('{}/{}.npz'.format(sdir,name))
        if hasattr(S,'headers'):
           snames=[i.sname for i in S.headers]
        else:
           snames=[S.info.sname,]

        #combine all shpfiles
        S0=zdata(); nrec=0; X=[]; Y=[]
        for i in arange(len(snames)):
            fname='{}.shp'.format(snames[i])
            if not os.path.exists(fname): continue

            #read each shp
            print('combining {}'.format(fname))
            S=read_shapefile_data(fname)
            S0.prj=S.prj; S0.type=S.type

            #read mth river
            for m in arange(S.nrec):
                xi=S.xy[m][:,0]; yi=S.xy[m][:,1]
                if not isnan(xi[0]): xi=r_[nan,xi]; yi=r_[nan,yi]
                if not isnan(xi[-1]): xi=r_[xi,nan]; yi=r_[yi,nan]
                nind=nonzero(isnan(xi))[0]; nsec=len(nind)-1

                #read kth section
                for k in arange(nsec):
                    i1=nind[k]+1; i2=nind[k+1]
                    xii=xi[i1:i2]; yii=yi[i1:i2]; slen=len(xii)

                    #re-sample river
                    if npt_sample is not None:
                       if slen<(2*npt_sample): continue
                       npt=int(ceil(slen/npt_sample))
                       pind=slen*linspace(0,1,npt+2)[1:-1]; pind=pind.astype('int')
                       xii=xii[pind]; yii=yii[pind]

                    #smooth
                    if npt_smooth is not None:
                       xii=smooth(xii,npt_smooth)
                       yii=smooth(yii,npt_smooth)

                    #add
                    nrec=nrec+1
                    xii=r_[nan,xii]; yii=r_[nan,yii]
                    X.extend(xii); Y.extend(yii)

        #write shapefile
        S0.xy=c_[array(X),array(Y)]
        write_shapefile_data('{}/{}'.format(sdir,name),S0)

    #def write_shapefile(self,sname,data='rivers',crs='epsg:4326',npt_smooth=None):
    def write_shapefile(self,data='rivers',crs='epsg:4326',npt_smooth=None):
        '''
        data: string name or data
        sname: shapefile name
        '''

        #get data
        SF=zdata()
        if isinstance(data, str):
            exec('SF.data=self.{}'.format(data))
        else:
            SF.data=data

        #get xy
        SF.xy=[];
        for i in arange(len(SF.data)):
            yi,xi=self.get_coordinates(SF.data[i])
            if len(xi)==0: continue
            fp=(xi==self.info.nodata)|(yi==self.info.nodata);
            xi[fp]=nan; yi[fp]=nan;
            SF.xy.append(c_[xi,yi])
        SF.xy=array(SF.xy, dtype='object')
        if SF.xy.ndim==3: SF.xy=[*SF.xy]

        if npt_smooth is not None:
            SF.xy=self.smooth_river(SF.xy,npt_smooth=npt_smooth)

        #geopandas
        total = 0
        branches = []
        for i, river in enumerate(SF.xy):
            idxs=np.squeeze(np.argwhere(np.isnan(river[:,0])))
            #breakpoint()
            try:
                #The main channel has tributaries
                n_branches = idxs.shape[0]

                intersections = []
                for j in np.arange(n_branches):
                    if j == 0:
                        coords = river[0:idxs[j]]
                        mainchannel = LineString(coords)
                    else:
                        total += 1
                        coords = river[idxs[j-1]+1:idxs[j]]
                        branch = LineString(coords)
                        branches.append({'reach_id': total, 'geometry': branch})

                        if len(mainchannel.intersection(branch).coords[:]) > 0:
                            intersections.append(mainchannel.intersection(branch))
                        #breakpoint()
                        #mainchannel = ssplit(mainchannel, intersection)
                    #df = pd.DataFrame({'longitude': coords[:,0], 'latitude': coords[:,1]})
                    #df.dropna(inplace=True)
                    #breakpoint()
                    #gdf = gpd.GeoDataFrame(df)
                    #gdf.to_file(f'branch_{i}.shp')
                
                print(total)
                channel = ssplit(mainchannel, MultiPoint(intersections))
                for branch in channel.geoms:
                    total += 1
                    branches.append({'feature_id':total, 'geometry': branch})
            except:
                #only have main channel no tributaries
                total += 1
                coords = river[0:-1]
                branches.append({'feature_id':total, 'geometry': LineString(coords)})

        #breakpoint()
        gdf = gpd.GeoDataFrame(branches)
        gdf = gdf.set_crs(crs)
        #gdf.to_file(sname)
        return gdf


        #if npt_smooth is not None:
        #    SF.xy=self.smooth_river(SF.xy,npt_smooth=npt_smooth)

        ##write shapefile
        #write_shapefile_data(sname,SF)

    def smooth_river(self,data,npt_smooth=5):

        #parameter
        ds=self.info.ds; nodata=self.info.nodata

        #for each river
        for i in arange(len(self.rivers)):
            river=self.rivers[i]; rxy=data[i].copy()
            sind_unique,fpu,npt_unique=unique(river,return_inverse=True,return_counts=True)
            npt=npt_unique[fpu]; npt[river==nodata]==-1;

            #for each section
            sids=nonzero((npt>1)|(npt==-1))[0]

            for m in arange(len(sids)):
                #get subsection indices
                if m==0:
                    id1=0;
                else:
                    id1=sids[m-1];
                id2=sids[m]
                if river[id1]==nodata: id1=id1+1
                if river[id2]!=nodata: id2=id2+1

                #modify smooth pts
                if (id2-id1)<3:
                    continue
                elif (id2-id1)<(npt_smooth+2):
                    npt_smooth_final=(id2-id1)-2
                else:
                    npt_smooth_final=npt_smooth

                #smooth river section, keep front and end pts unchanged
                srxi=smooth(rxy[id1:id2,0].copy(),npt_smooth_final); sryi=smooth(rxy[id1:id2,1].copy(),npt_smooth_final)
                rxy[(id1+1):(id2-1),0]=srxi[1:-1]; rxy[(id1+1):(id2-1),1]=sryi[1:-1];
            data[i]=rxy

        return squeeze(array(data))

    def compute_sind_ext(self,sdir='.',dis_limit=None):
        '''
        compute additional acc at the bounary, these acc are from other DEMs
        '''
        self.info.sind_ext=[]
        if len(self.info.nbs)==0: return
        if not hasattr(self.info,'sind0'): return
        if not hasattr(self.info,'dem_bnd0'): return

        #determin dis_limit
        if dis_limit is None:
            dis_limit=2*self.info.header[-1]

        #collect bnd info from other DEMs
        sind_all=[]; dem_all=[]; acc_all=[]; sx_all=[]; sy_all=[];
        for i in arange(len(self.info.nbs)):
            fname='{}/{}_{}.npz'.format(sdir,self.info.sname0,self.info.ids[self.info.nbs[i]])
            flag_loop=True
            while flag_loop:
                  try:
                     sinfo=self.read_data(fname,'info')
                     flag_loop=False
                  except:
                     time.sleep(5)
            if not hasattr(sinfo,'sind0'): continue
            sind=sinfo.sind0; yi,xi=self.get_coordinates(sind,header=sinfo.header)
            sind_all.extend(sind); dem_all.extend(sinfo.dem0); acc_all.extend(sinfo.acc0)
            sx_all.extend(xi); sy_all.extend(yi)
        sind_all=array(sind_all); dem_all=array(dem_all); acc_all=array(acc_all); sx_all=array(sx_all); sy_all=array(sy_all)
        if len(sind_all)==0: return

        #local sind_bnd
        sind=self.info.sind_bnd; dem=self.info.dem_bnd; dem0=self.info.dem_bnd0; ly,lx=self.get_coordinates(sind)

        #exclude pts too far
        fp=~((sx_all<(min(lx)-dis_limit))|(sx_all>(max(lx)+dis_limit))|(sy_all<(min(ly)-dis_limit))|(sy_all>(max(ly)+dis_limit)))
        if sum(fp)==0: return
        sind_all=sind_all[fp]; dem_all=dem_all[fp]; acc_all=acc_all[fp]; sx_all=sx_all[fp]; sy_all=sy_all[fp]

        s=zdata(); s.p0=c_[sx_all,sy_all]; s.p=c_[lx,ly]; save_npz('tmp_near',s)

        #index of nearest pts
        ids=near_pts(c_[sx_all,sy_all],c_[lx,ly])
        sdist=abs((lx[ids]-sx_all)+1j*(ly[ids]-sy_all))
        fp=(sdist<dis_limit)*(dem0[ids]<dem_all)

        acc_all_ext=acc_all[fp]; sind_all_ext=sind[ids][fp]
        sind_ext,fpu,npt_unique=unique(sind_all_ext,return_inverse=True,return_counts=True)
        acc_ext=zeros(len(sind_ext)).astype('int')
        for i in arange(len(acc_all_ext)):
            acc_ext[fpu[i]]=acc_ext[fpu[i]]+acc_all_ext[i]

        self.info.sind_ext=sind_ext; self.info.acc_ext=acc_ext

        return

    def proc_demfile(self,names=None,ids=None,sname=None,sdir='.',findex=None,depth_limit=[-1e5,1e4],subdomain_size=1e6,offset=1, ResolveFlats=False, FillDepression=True, BreachDepression=False):
        '''
        headers: include information about each DEM file
        '''

        #read file informations
        if not os.path.exists('./{}/{}.npz'.format(sdir,sname)):
            self.read_files_info(names,ids,sname=sname,sdir=sdir,depth_limit=depth_limit,plot_domain=True)
        else:
            self.read_data('./{}/{}.npz'.format(sdir,sname))

        header = self.headers[0]
        t0=time.time()
        S=dem()
        #read diminfo
        S.read_deminfo(subdomain_size=subdomain_size,offset=offset,header=header,sdir=sdir)

            #if data size is smaller than subdomain_size, read the whole data
            # if prod(S.info.ds)<subdomain_size: S.read_demdata()

        if header.name.endswith('.npz'): S.dem=S.read_data(header.name,'dem')
        S.read_demdata()

        if ResolveFlats:
            rd.ResolveFlats(S.dem, in_place=True)

        if FillDepression:
            rd.FillDepressions(S.dem, epsilon=True, in_place=True, topology='D8')

        if BreachDepression:
            rd.BreachDepressions(S.dem, in_place=True)

        S.compute_dir()
        #update bnd info
        for svar in ['sind','dir','dem']:
            if not hasattr(S.info,'{}_bnd_local'.format(svar)): continue
            exec('S.info.{}_bnd=S.info.{}_bnd_local.copy(); delattr(S.info,"{}_bnd_local")'.format(svar,svar,svar))
                #S.info.sind_bnd=S.info.sind_bnd_local.copy(); delattr(S.info,'sind_bnd_local');
                #S.info.dir_bnd=S.info.dir_bnd_local.copy(); delattr(S.info,'dir_bnd_local');
                #S.info.dem_bnd=S.info.dem_bnd_local.copy(); delattr(S.info,'dem_bnd_local');

        #save DEMs info
        S.info.names=header.names; S.info.ids=header.ids; S.info.nbs=header.nbs
        S.info.sname0=header.sname0; S.info.sname=header.sname

        ##save acc at boundary
        #S.compute_watershed(ireturn=1)

        #save information
        S.save_data(S.info.sname,['dir','info'])

        dt=time.time()-t0
        print('total time={:.1f}s: {}'.format(time.time()-t0,S.info.sname))
        sys.stdout.flush()

        '''
        #compute sind_ext
        header=self.headers[0]; fname='{}.npz'.format(header.sname)

        #read processed DEM data
        S=dem(); sinfo=S.read_data(fname,'info')
        #if len(sinfo.nbs)==0: continue
        #if hasattr(sinfo,'sind_ext'): continue
        S.read_data(fname)

        #make sure neighboring DEMs are available
        flag_loop=True
        while flag_loop:
              flag_loop=False
              for i in arange(len(S.info.nbs)):
                  fname='{}/{}_{}.npz'.format(sdir,S.info.sname0,S.info.ids[S.info.nbs[i]])
                  if not os.path.exists(fname): flag_loop=True
              time.sleep(10)

        #compute sind_ext
        S.compute_sind_ext(sdir=sdir)

        #save information
        S.save_data(header.sname,['dir','info'])
        '''

        print('finish processing {}'.format(S.info.sname))
        S=None; sys.stdout.flush()

if __name__=="__main__":
    close('all')
