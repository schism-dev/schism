%Process CFSR for sflux_prc*
%Reformat prc files (monthly)
clear all; close all;

fill_in=-999; %junk value
yr_now=2012;
mon_now=4;
for imon=1:6 %monthly  
  if(mon_now>12)
    mon_now=mon_now-12;
    yr_now=yr_now+1;  
  end
  yr_char=num2str(yr_now);
  mon_char=sprintf('%2.2d',mon_now);
  if(yr_now<2011) %change in name etc from 2011
    fname3=['prate.gdas.' yr_char mon_char '.grb2.nc'];
  else
    fname3=['prate.cdas1.' yr_char mon_char '.grb2.nc'];
  end
  disp(['doing ' fname3]);

  %Read
  ncid0 = netcdf.open(fname3,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time= netcdf.getVar(ncid0, vid); %in hours
  vid=netcdf.inqVarID(ncid0,'lon'); 
  lon=netcdf.getVar(ncid0,vid); %1D
  vid=netcdf.inqVarID(ncid0,'lat');
  lat=netcdf.getVar(ncid0,vid); %1D
  vid=netcdf.inqVarID(ncid0,'PRATE_L1_Avg_1');
  prate=netcdf.getVar(ncid0,vid); 
  netcdf.close(ncid0);

  n1=size(prate,1);
  n2=size(prate,2);
  ntimes=size(prate,3);

  %Put lon in range
  itmp=find(lon>180);
  lon(itmp)=lon(itmp)-360;

  %Transform lon/lat into 2D arrays; flip later if the quad is not counter-clockwise
  for j=1:n2
    lon2(1:n1,j)=lon(1:n1);
  end
  for i=1:n1
    lat2(i,1:n2)=lat(1:n2);
  end
  area=(lon2(2,1)-lon2(1,1))*(lat2(1,2)-lat2(1,1));
  if(imon==1 && area>=0)
    error('A>=0');
  end

  %Write nc file of SCHISM format (Warning: use single precision!)
  %Write out new nc file (1 month per file)

  %prc
  filenm=['sflux_prc_' yr_char mon_char '.nc'];
  ncid2 = netcdf.create(filenm,'CLOBBER');
  dims(1)=netcdf.defDim(ncid2,'lon',n1);
  dims(2)=netcdf.defDim(ncid2,'lat',n2);
  dims(3)=netcdf.defDim(ncid2,'time',ntimes);
  % Define new variable in the file.
  timeid = netcdf.defVar(ncid2,'time','float',dims(3));
  lonid = netcdf.defVar(ncid2,'lon','float',dims(1:2));
  latid = netcdf.defVar(ncid2,'lat','float',dims(1:2));
  lid = netcdf.defVar(ncid2,'prate','float',dims(1:3));

  % Leave define mode and enter data mode to write data.
  netcdf.endDef(ncid2)

  netcdf.putVar(ncid2,timeid,single(time/24.));
  % Re-enter define mode.
  netcdf.reDef(ncid2);
  % Create an attribute associated with the variable.
  netcdf.putAtt(ncid2,timeid,'base_date',int32([yr_now mon_now 1 0])); %must use int32
  netcdf.endDef(ncid2)

  netcdf.putVar(ncid2,lonid,lon2(:,n2:-1:1));
  netcdf.putVar(ncid2,latid,lat2(:,n2:-1:1));
  netcdf.putVar(ncid2,lid,prate(:,n2:-1:1,:));
  netcdf.close(ncid2);

  mon_now=mon_now+1;

  clear prate lon2 lat2 time;
end %for day
