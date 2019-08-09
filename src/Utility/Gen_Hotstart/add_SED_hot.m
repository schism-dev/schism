%Add SED vars into hotstart.nc from a pure hydro run
%Input: depth.mlab (1st part of hgrid.gr3 to get depth info); 
%       hotstart.nc.0 (pure hydro); some constants below
%Output: hotstart.nc
clear all; close all;

nSED=3; %# of SED classes
SED_MBEDP=3; %# of bed prop
SED_Nbed=1; %# of bed layers
ncid0 = netcdf.open('hotstart.nc.0','NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid0);

for i=0:ndims-1 %netcdf.getVar follows C convention (starting from 0)
  [dimname{i+1}, dimlen(i+1)] = netcdf.inqDim(ncid0,i);
end %for i
np=dimlen(1); ne=dimlen(2);
nvrt=dimlen(4);

%Define some SED vars
%1st prop: thickness; 2nd: age; 3rd: porosity
bed(1,1:SED_Nbed,1:ne)=5./SED_Nbed;
bed(2,1:SED_Nbed,1:ne)=0;
bed(3,1:SED_Nbed,1:ne)=0.4;
bedfrac(1,1:SED_Nbed,1:ne)=0.3; %frac for class 1
bedfrac(2,1:SED_Nbed,1:ne)=0.3; %frac for class 2
bedfrac(nSED,1:SED_Nbed,1:ne)=0.40; %frac for class 3

for i=0:nvars-1 %netcdf.getVar follows C convention (starting from 0)
  [varname{i+1},xtype(i+1),dimids0(i+1).ids,natts] = netcdf.inqVar(ncid0,i);
  var(i+1).data= netcdf.getVar(ncid0, i);
end %for i

netcdf.close(ncid0);

%Add SED conc to tracer arrays
varout(1:8)=var(1:8);
varout(10:11)=var(10:11);
varout(14:nvars)=var(14:end);
varout(9).data=zeros(dimlen(5)+nSED,dimlen(4),dimlen(2));
varout(12).data=zeros(dimlen(5)+nSED,dimlen(4),dimlen(1));
varout(13).data=zeros(dimlen(5)+nSED,dimlen(4),dimlen(1));
varout(9).data(1:2,:,:)=var(9).data;
varout(12).data(1:2,:,:)=var(12).data;
varout(13).data(1:2,:,:)=var(13).data;

varout(9).data(3:end,:,:)=0; %SED conc
varout(12).data(3:end,:,:)=0; %SED conc
varout(13).data(3:end,:,:)=0; %SED conc

%Output
ncid2 = netcdf.create('hotstart.nc','CLOBBER');
%Def 
dimlen(5)=dimlen(5)+nSED;
for i=1:ndims
  id2=netcdf.defDim(ncid2,dimname{i},dimlen(i)); %dimids are 0-based
end %for i

for i=0:nvars-1 %netcdf.getVar follows C convention (starting from 0)
  id=netcdf.defVar(ncid2,varname{i+1},xtype(i+1),dimids0(i+1).ids);
end %for i
netcdf.endDef(ncid2);

for i=0:nvars-1 %netcdf.getVar follows C convention (starting from 0)
  netcdf.putVar(ncid2,i,varout(i+1).data);
end %for i

netcdf.reDef(ncid2);

%Other SED vars
id_ntr=netcdf.defDim(ncid2,'SED_ntr',nSED);
id_MBEDP=netcdf.defDim(ncid2,'SED_MBEDP',SED_MBEDP);
id_Nbed=netcdf.defDim(ncid2,'SED_Nbed',SED_Nbed);

dims(1)=0; %np
vid_dp=netcdf.defVar(ncid2,'SED3D_dp','double',dims(1));
vid_rough=netcdf.defVar(ncid2,'SED3D_rough','double',dims(1));
dims(1)=id_MBEDP; dims(2)=id_Nbed; dims(3)=1;
vid_bed=netcdf.defVar(ncid2,'SED3D_bed','double',dims(1:3));
dims(1)=id_ntr; %dims(2)=id_Nbed; dims(3)=
vid_bedfrac=netcdf.defVar(ncid2,'SED3D_bedfrac','double',dims(1:3));
netcdf.endDef(ncid2);

tmp=load('depth.mlab');
netcdf.putVar(ncid2,vid_dp,tmp(:,end));
rough(1:np)=8.e-4;
netcdf.putVar(ncid2,vid_rough,rough);
netcdf.putVar(ncid2,vid_bed,bed);
netcdf.putVar(ncid2,vid_bedfrac,bedfrac);

netcdf.close(ncid2);
