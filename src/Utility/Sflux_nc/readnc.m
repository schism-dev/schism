%Read in 1 netcdf file into matlab and plot; modify as appropriate
clear all; close all;
%Warning: dimension for arrays reversed from ncdump (FORTRAN convention)

fill_in=-999; %junk value
ncid = netcdf.open('2003-09-16/00z/03091600_nmm_d01.GrbF00000.nc','NC_NOWRITE');
vid=netcdf.inqVarID(ncid,'g3_lat_0'); %input var./array name
lat = netcdf.getVar(ncid, vid); 
vid=netcdf.inqVarID(ncid,'g3_lon_1');
lon = netcdf.getVar(ncid, vid); 
vid=netcdf.inqVarID(ncid,'U_GRD_GDS3_HTGL');
uwind = netcdf.getVar(ncid, vid); 
uwind(find(uwind<fill_in+1))=nan;

vid=netcdf.inqVarID(ncid,'V_GRD_GDS3_HTGL');
vwind = netcdf.getVar(ncid, vid); 
vwind(find(vwind<fill_in+1))=nan;

vid=netcdf.inqVarID(ncid,'TMP_GDS3_SFC'); %air temp. in K
airt= netcdf.getVar(ncid, vid); 
airt(find(airt<fill_in+1))=nan;

vid=netcdf.inqVarID(ncid,'PRMSL_GDS3_MSL'); %air pressure at MSL (Pa)
pres= netcdf.getVar(ncid, vid); 
pres(find(pres<fill_in+1))=nan;

vid=netcdf.inqVarID(ncid,'PRATE_GDS3_SFC'); %precip. rate (kg/m^2/s)
prc= netcdf.getVar(ncid, vid); 
prc(find(prc<fill_in+1))=nan;

vid=netcdf.inqVarID(ncid,'LAND_GDS3_SFC'); %Land-sea mask (land=1;sea=0)
mask= netcdf.getVar(ncid, vid); 

%plot(lon,lat);
%contour(lon,lat,pres);
%colorbar;
hold on;
%Subsample
ind1=10:10:size(lon,1);
ind2=10:10:size(lon,2);
quiver(lon(ind1,ind2),lat(ind1,ind2),uwind(ind1,ind2),vwind(ind1,ind2),0,'k');
quiver(-77.9,40.5,5,0.,0,'r');

