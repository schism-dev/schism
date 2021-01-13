%Author: Joseph Zhang
%Generate ampl. and phases for SAL (self-attraction loading tides) from FES2014 using linear interp directly
%from nc files.
%Requires inputs: (1) hgrid.ll (depth not used)
%                 (2) *.nc (under ./fes2014a_loadtide/load_tide/)
%  Need to first download FES2014 at AVISO site
%Outputs: loadtide_[FREQ].gr3 (lon in [0,360))
clear all; close all;

const={'m2','s2','n2','k2','k1','p1','o1','q1'};

%Read in hgrid
fid=fopen('hgrid.ll','r');
char=fgetl(fid);
tmp1=str2num(fgetl(fid));
fclose(fid);

ne=fix(tmp1(1));
np=fix(tmp1(2));

fid=fopen('hgrid.ll','r');
%Change here if there are >1 'depths'
c1=textscan(fid,'%d%f%f%f',np,'headerLines',2);
fclose(fid);
fid=fopen('hgrid.ll','r');
c2=textscan(fid,'%d%d%d%d%d%d',ne,'headerLines',2+np);
fclose(fid);

xl0=double(c1{2}(:));
yl=double(c1{3}(:));
%dp=c1{4}(:);
i34=fix(c2{2}(:));

nm(1:ne,1:4)=nan;
for i=1:ne
  for j=1:i34(i)
    nm(i,j)=fix(c2{j+2}(i));
  end %for j
end %for i

%Make lon in [0,360]
xl=xl0;
indx=find(xl<0);
xl(indx)=xl(indx)+360;

%miss_value=1.e10; %junk value
for i=1:length(const)
  ncid = netcdf.open(['fes2014a_loadtide/load_tide/' const{i} '.nc'],'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid,'lat'); %1D
  lat0 = netcdf.getVar(ncid, vid); %ascending order [-90,90]
  vid=netcdf.inqVarID(ncid,'lon');
  lon0 = netcdf.getVar(ncid, vid); %ascending order [0,360)
  nx=length(lon0); ny=length(lat0);

  clear amp0 pha0 lon2;
  amp0=zeros(nx+1,ny);
  pha0=amp0;
  lon2=nan(nx+1);
  lon2(1:nx)=lon0;

  vid=netcdf.inqVarID(ncid,'amplitude'); %(nx,ny)
  amp0(1:nx,:)= netcdf.getVar(ncid, vid);
  amp0=amp0/100; %to meters 
  vid=netcdf.inqVarID(ncid,'phase');
  pha0(1:nx,:)= netcdf.getVar(ncid, vid); %degr
  netcdf.close(ncid);

  %Close lon gap near 0,360
  lon2(nx+1)=360;
  amp0(nx+1,:)=amp0(1,:);
  pha0(nx+1,:)=pha0(1,:);

  %Fill in junks
  amp0(find(abs(amp0)>1.e10))=0;
  pha0(find(abs(pha0)>1.e10))=0;

  %Make all 1D vectors
  amp=double(reshape(amp0,(nx+1)*ny,1));
  pha=double(reshape(pha0,(nx+1)*ny,1));
  count=0;
  for iy=1:ny
    for ix=1:nx+1
       count=count+1;
       lon(count)=double(lon2(ix));
       lat(count)=double(lat0(iy));
    end %for ix
  end %for iy

  ampout=griddata(lon,lat,amp,xl,yl);
  %Use nearest approach to avoid phase wrap-around
  phaout=griddata(lon,lat,pha,xl,yl,'nearest');

%  scatter(xl,yl,6,ampout*100,'filled');
%  colormap(jet(40));
%  colorbar;

  %Output
  fid=fopen(['loadtide_' upper(const{i}) '.gr3'],'w');
  fprintf(fid,'%s\n',const{i});
  fprintf(fid,'%d %d\n',ne,np);
  fprintf(fid,'%d %f %f %f %f\n',[(1:np)' xl yl ampout phaout]');
  %fprintf(fid,'%d %f %f %f %f\n',[(1:np)' xl0 yl ampout phaout]');
  fprintf(fid,'%d %d %d %d %d %d\n',[(1:ne)' i34 nm(:,:)]');
  fclose(fid);
end %for i - freq's
