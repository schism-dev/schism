function []=check_schism_nc(fname,itype)
% Authors: Joseph Zhang
% Date: April 2023
%Check SCHISM's netcdf inputs including b.c. (*th.nc), nudging (*_nu.nc) and i.c. (hotstart.nc)
%Inputs:
%       fname: nc file name (e.g. 'TEM_nu.nc')
%       itype: 1: b.c.; 2: nudging; 3: i.c.

%Check inputs
if(itype==1) %b.c.
  if(length(regexp(fname,'th.nc'))==0) 
    error(['inputs not consistent:' fname ';' num2str(itype)]);
  end
  ncid0 = netcdf.open(fname,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time_step');
  dt=double(netcdf.getVar(ncid0, vid));
  disp('Time step='); dt
  vid=netcdf.inqVarID(ncid0,'time_series');
  var= double(netcdf.getVar(ncid0, vid)); %(nComponents,nvrt,nond,ntimes)
  netcdf.close(ncid0);
  av=mean(mean(mean(mean(var))));
  allmin=min(min(min(min(var))));
  allmax=max(max(max(max(var))));
  disp('Global mean,min,max=');
  [av allmin allmax]
elseif(itype==2) %nudging
  if(length(regexp(fname,'_nu.nc'))==0) 
    error(['inputs not consistent:' fname ';' num2str(itype)]);
  end

  ncid0 = netcdf.open(fname,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time');
  time=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'map_to_global_node');
  map=netcdf.getVar(ncid0, vid);
  
%  disp('Time step='); dt
  vid=netcdf.inqVarID(ncid0,'tracer_concentration');
  var= double(netcdf.getVar(ncid0, vid)); %(ntracers,nvrt,node_nu,ntimes)
  netcdf.close(ncid0);
  av=mean(mean(mean(mean(var))));
  allmin=min(min(min(min(var))));
  allmax=max(max(max(max(var))));
  disp('Global mean,min,max=');
  [av allmin allmax]
elseif(itype==3)
  error('not done');
else
  error('Unknown itype');
end

