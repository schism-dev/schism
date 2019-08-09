%Process CFSR files to sflux_air*
%Reformat air files and interpolation for prmsl (monthly files)
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
  fname=['wnd10m.cdas1.' yr_char mon_char '.grb2.nc'];
  fname2=['tmp2m.cdas1.' yr_char mon_char '.grb2.nc'];
  fname3=['q2m.cdas1.' yr_char mon_char '.grb2.nc'];
  disp(['doing ' fname ' etc']);

  %Read
  ncid0 = netcdf.open(fname,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time= netcdf.getVar(ncid0, vid); %in hours
  vid=netcdf.inqVarID(ncid0,'lon'); 
  lon=netcdf.getVar(ncid0,vid); %1D
  vid=netcdf.inqVarID(ncid0,'lat');
  lat=netcdf.getVar(ncid0,vid); %1D
  vid=netcdf.inqVarID(ncid0,'U_GRD_L103');
  uwind=netcdf.getVar(ncid0,vid); 
  vid=netcdf.inqVarID(ncid0,'V_GRD_L103');
  vwind=netcdf.getVar(ncid0,vid); 
  netcdf.close(ncid0);

  ncid0 = netcdf.open(fname2,'NC_NOWRITE'); 
  vid=netcdf.inqVarID(ncid0,'TMP_L103');
  stmp=netcdf.getVar(ncid0,vid);
  netcdf.close(ncid0);

  ncid0 = netcdf.open(fname3,'NC_NOWRITE'); 
  vid=netcdf.inqVarID(ncid0,'SPF_H_L103');
  spfh=netcdf.getVar(ncid0,vid);
  netcdf.close(ncid0);

  %Assume all variables have same dim. except prmsl
  n1=size(uwind,1);
  n2=size(uwind,2);
  ntimes=size(uwind,3);

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

  %Read in pressure
  fname=['prmsl.cdas1.' yr_char mon_char '.grb2.nc'];
  disp(['doing ' fname]);
  ncid0 = netcdf.open(fname,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time3= netcdf.getVar(ncid0, vid); %in hours
  vid=netcdf.inqVarID(ncid0,'lon');
  lon3=netcdf.getVar(ncid0,vid); %1D
  vid=netcdf.inqVarID(ncid0,'lat');
  lat3=netcdf.getVar(ncid0,vid); %1D
  vid=netcdf.inqVarID(ncid0,'PRMSL_L101');
  prmsl=netcdf.getVar(ncid0,vid);
  netcdf.close(ncid0);

  if(length(time3) ~= length(time) || ntimes ~=length(time3))
    error('prmsl time does not match rest');
  end

  %Transform lon/lat into 2D arrays
  n1_3=size(prmsl,1);
  n2_3=size(prmsl,2);
  %Put lon in range
  itmp=find(lon3>180);
  lon3(itmp)=lon3(itmp)-360;

  for j=1:n2_3
    lon3p(1:n1_3,j)=lon3(1:n1_3);
  end
  for i=1:n1_3
    lat3p(i,1:n2_3)=lat3(1:n2_3);
  end
  area=(lon3p(2,1)-lon3p(1,1))*(lat3p(1,2)-lat3p(1,1));
  if(imon==1 && area>=0)
    error('A>=0 (2)');
  end

  %Interpolate prmsl to the other grid
  for it=1:ntimes
    prmsl2(:,:,it)=griddata(double(lon3p),double(lat3p),double(prmsl(:,:,it)),double(lon2),double(lat2),'nearest');
  end %for

  %Write nc file of SCHISM format (Warning: use single precision!)
  %Write out new nc file (1 month per file)
  filenm=['sflux_air_' yr_char mon_char '.nc'];
  ncid2 = netcdf.create(filenm,'CLOBBER');
  dims(1)=netcdf.defDim(ncid2,'lon',n1);
  dims(2)=netcdf.defDim(ncid2,'lat',n2);
  dims(3)=netcdf.defDim(ncid2,'time',ntimes);

  %Define new variable in the file.
  timeid = netcdf.defVar(ncid2,'time','float',dims(3));
  lonid = netcdf.defVar(ncid2,'lon','float',dims(1:2));
  latid = netcdf.defVar(ncid2,'lat','float',dims(1:2));
  uid = netcdf.defVar(ncid2,'uwind','float',dims(1:3));
  vid = netcdf.defVar(ncid2,'vwind','float',dims(1:3));
  tid = netcdf.defVar(ncid2,'stmp','float',dims(1:3));
  hid = netcdf.defVar(ncid2,'spfh','float',dims(1:3));
  pid = netcdf.defVar(ncid2,'prmsl','float',dims(1:3));

  % Leave define mode and enter data mode to write data.
  netcdf.endDef(ncid2)

  netcdf.putVar(ncid2,timeid,single(time)/24.);
  % Re-enter define mode.
  netcdf.reDef(ncid2);
  % Create an attribute associated with the variable.
  netcdf.putAtt(ncid2,timeid,'base_date',int32([yr_now mon_now 1 0])); %must use int32
  netcdf.endDef(ncid2)

  netcdf.putVar(ncid2,lonid,lon2(:,n2:-1:1));
  netcdf.putVar(ncid2,latid,lat2(:,n2:-1:1));
  netcdf.putVar(ncid2,uid,uwind(:,n2:-1:1,:));
  netcdf.putVar(ncid2,vid,vwind(:,n2:-1:1,:));
  netcdf.putVar(ncid2,tid,stmp(:,n2:-1:1,:));
  netcdf.putVar(ncid2,hid,spfh(:,n2:-1:1,:));
  netcdf.putVar(ncid2,pid,single(prmsl2(:,n2:-1:1,:)));
  netcdf.close(ncid2);

  mon_now=mon_now+1;

  clear uwind vwind stmp spfh lon* lat* time prmsl*;
end %for day
