%Change x,y coordinates of out2d*.nc (e.g. re-projection). Run this outside outputs/.
%Example here for polar stereographic proj for Arctic grid
%Inputs: outputs/out2d*; some parameters below
%Outputs: outputs_new/out2d*  (symlink zCo*.nc and others for visIT)
clear all; close all;

r0=6378206.4; %earth radius
nstacks=[10:20]; %stacks to be processed

if(exist('outputs_new')); status=rmdir('outputs_new','s'); end;
mkdir outputs_new;

for istack=nstacks
  disp(['doing stack ' num2str(istack)]);
  status=copyfile(['outputs/out2d_' num2str(istack) '.nc'],'outputs_new');
  ncid=netcdf.open(['outputs_new/out2d_' num2str(istack) '.nc'],'WRITE');
  lonvid=netcdf.inqVarID(ncid,'SCHISM_hgrid_node_x');
  lon= netcdf.getVar(ncid, lonvid);
  latvid=netcdf.inqVarID(ncid,'SCHISM_hgrid_node_y');
  lat= netcdf.getVar(ncid, latvid);
  
  lon2=lon/180*pi;
  lat2=lat/180*pi;
  xg=r0*cos(lat2).*cos(lon2);
  yg=r0*cos(lat2).*sin(lon2);
  zg=r0*sin(lat2);
  
  xp=2*xg*r0./(r0+zg);
  yp=2*yg*r0./(r0+zg);
  
  netcdf.putVar(ncid,lonvid,xp);
  netcdf.putVar(ncid,latvid,yp);
  
  netcdf.close(ncid);
end %for istack
