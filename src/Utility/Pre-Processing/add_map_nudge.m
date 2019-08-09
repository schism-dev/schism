%Add mapping array into *_nu.nc for new SCHISM version
%Inputs: nu0.nc (old _nu.nc)
%Outputs: out.nc
clear all; close all;

in='nu0.nc';
ncid0 = netcdf.open(in,'NC_NOWRITE');
vid=netcdf.inqVarID(ncid0,'tracer_concentration'); %input var./array name
tr= netcdf.getVar(ncid0, vid);
[ntr,nvrt,np,ntimes]=size(tr)
vid=netcdf.inqVarID(ncid0,'time');
time= double(netcdf.getVar(ncid0, vid));
netcdf.close(ncid0);

ncid2 = netcdf.create('out.nc','CLOBBER');
dims(1)=netcdf.defDim(ncid2,'node',np);
dims(2)=netcdf.defDim(ncid2,'nLevels',nvrt); %not really used
dims(3)=netcdf.defDim(ncid2,'ntracers',ntr); %not really used
dims(4)=netcdf.defDim(ncid2,'time',ntimes);

% Define new variable in the file.
timeid = netcdf.defVar(ncid2,'time','NC_DOUBLE',dims(4));
mid = netcdf.defVar(ncid2,'map_to_global_node','NC_INT',dims(1));
tid = netcdf.defVar(ncid2,'tracer_concentration','NC_FLOAT',dims([3 2 1 4]));
netcdf.endDef(ncid2)

netcdf.putVar(ncid2,timeid,double(time));
map=1:np;
netcdf.putVar(ncid2,mid,fix(map));
netcdf.putVar(ncid2,tid,single(tr));
netcdf.close(ncid2);
