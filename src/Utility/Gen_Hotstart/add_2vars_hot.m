%Add 2 new vars in hotstart for newer versions than July 2021 (added 2 new vars nsteps_from_cold and cumsum_eta)
%Should work with optional modules as well.
%Input: hotstart.nc.0 (old ofrmat)
%Output: hotstart.nc
clear all; close all;

ncid0 = netcdf.open('hotstart.nc.0','NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid0);

for i=0:ndims-1 %netcdf.getVar follows C convention (starting from 0)
  [dimname{i+1}, dimlen(i+1)] = netcdf.inqDim(ncid0,i);
end %for i
np=dimlen(1); ne=dimlen(2);
nvrt=dimlen(4);

for i=0:nvars-1 %netcdf.getVar follows C convention (starting from 0)
  [varname{i+1},xtype(i+1),dimids0(i+1).ids,natts] = netcdf.inqVar(ncid0,i);
  var(i+1).data= netcdf.getVar(ncid0, i);
end %for i
netcdf.close(ncid0);

%varout(1:nvars)=var(1:end);

%Output
ncid2 = netcdf.create('hotstart.nc','CLOBBER');
%Def 
for i=1:ndims
  id2=netcdf.defDim(ncid2,dimname{i},dimlen(i)); %dimids are 0-based
end %for i

for i=0:nvars-1 %netcdf.getVar follows C convention (starting from 0)
  id=netcdf.defVar(ncid2,varname{i+1},xtype(i+1),dimids0(i+1).ids);
end %for i
netcdf.endDef(ncid2);

for i=0:nvars-1 %netcdf.getVar follows C convention (starting from 0)
  netcdf.putVar(ncid2,i,var(i+1).data);
end %for i

netcdf.reDef(ncid2);

dims(1)=5; %one
vid_nstep=netcdf.defVar(ncid2,'nsteps_from_cold','NC_INT',dims(1));
dims(1)=0; %np
vid_cum=netcdf.defVar(ncid2,'cumsum_eta','NC_DOUBLE',dims(1));
netcdf.endDef(ncid2);

cumsum(1:np)=0.;
netcdf.putVar(ncid2,vid_cum,cumsum);
netcdf.putVar(ncid2,vid_nstep,0);
netcdf.close(ncid2);
