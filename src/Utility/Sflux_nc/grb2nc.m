%From YC
%Read in an existing netcdf files (converted from grib2) and output SCHISM nc input; need to 
% add file names below; modify other parts as approx. (e.g. start time stamp)
clear all; close all;
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)

%netcdf files converted from grib2
%Each nc file is for one month
fill_in=-999; %junk value
%filenm={}; %define file names here;
ncid = netcdf.open('wnd10m.gdas.200704.nc','NC_NOWRITE');
ncid2 = netcdf.open('prmsl.gdas.200704.nc','NC_NOWRITE');

%read wnd data
vid1=netcdf.inqVarID(ncid,'latitude'); %input var./array name
lat = netcdf.getVar(ncid, vid1); 

vid2=netcdf.inqVarID(ncid,'longitude');
lon = netcdf.getVar(ncid, vid2); 
for i=1:length(lon)
  if(lon(i)>180)
      lon(i,1)=lon(i,1)-360;
  end    
end

vid3=netcdf.inqVarID(ncid,'time');
time = netcdf.getVar(ncid, vid3); 

vid4=netcdf.inqVarID(ncid,'UGRD_10maboveground');
uwind = netcdf.getVar(ncid, vid4); 

vid5=netcdf.inqVarID(ncid,'VGRD_10maboveground');
vwind = netcdf.getVar(ncid, vid5); 

[n1]=size(lat); 
[n2]=size(lon); 
[n3]=size(time); 

lon_x=zeros(length(lon),length(lat));
lat_y=zeros(length(lon),length(lat));

for i=1:length(lon)
  for j=1:length(lat)
     lon_x(i,j)=lon(i,1);
     lat_y(i,j)=lat(j,1);
  end    
end    

%read prmsl data
vid1p=netcdf.inqVarID(ncid2,'latitude'); %input var./array name
lat_p = netcdf.getVar(ncid2, vid1p); 

vid2p=netcdf.inqVarID(ncid2,'longitude');
lon_p = netcdf.getVar(ncid2, vid2p); 
for i=1:length(lon_p)
  if(lon_p(i)>180)
      lon_p(i,1)=lon_p(i,1)-360;
  end    
end

vid3p=netcdf.inqVarID(ncid2,'time');
time_p = netcdf.getVar(ncid2, vid3p); 

vid4p=netcdf.inqVarID(ncid2,'PRMSL_meansealevel');
prmsl = netcdf.getVar(ncid2, vid4p); 

lonp_x=zeros(length(lon_p),length(lat_p));
latp_y=zeros(length(lon_p),length(lat_p));
prmsl_org=zeros(length(lon_p)*length(lat_p),3);

for i=1:length(lon_p)
  for j=1:length(lat_p)
     lonp_x(i,j)=lon_p(i,1);
     latp_y(i,j)=lat_p(j,1);
  end    
end 

%convert prmsl grid to wind field grid
n=0;
for i=1:length(lon_p)
  for j=1:length(lat_p)
   n=n+1;   
   prmsl_org(n,1)=lonp_x(i,j);
   prmsl_org(n,2)=latp_y(i,j);
  end    
end   

pres=zeros(n2(1,1), n1(1,1),length(time));

for ii=1:length(time)
  n=0;
  for i=1:length(lon_p)
    for j=1:length(lat_p)
     n=n+1;   
     prmsl_org(n,3)=prmsl(i,j,ii);
    end    
  end 
  disp(strcat('time: ',num2str(ii)));
  F = TriScatteredInterp(prmsl_org(:,1),prmsl_org(:,2),prmsl_org(:,3));
  prmsl_w = F(lon_x,lat_y);
  pres(:,:,ii)=prmsl_w;
end


  airt=zeros(n2(1,1), n1(1,1), 24);
  %pres=zeros(n2(1,1), n1(1,1), 24);
  prc=zeros(n2(1,1), n1(1,1), 24);
  
  airt(:,:,:)=300; %nan;
%  pres(:,:,:)=1.e5; %nan;
  prc(:,:,:)=0; %nan;

%airt=zeros(n2(1,1), n1(1,1), n3(1,1));
%pres=zeros(n2, n1, n3);
%prc=zeros(n2, n1, n3);

%airt(:,:,:)=300; %nan;
%pres(:,:,:)=1.e5; %nan;
%prc(:,:,:)=0; %nan;

for i=1:length(time) %length(time)
  day=floor((i-1)/24+1);  
  hour=i-1-(day-1)*24; %0-23 hours
  
  disp(strcat('Done day cal_ ',num2str(day),' hour cal_ ',num2str(hour)));
  
  
  %Deal with junks
  uwind(find(uwind<fill_in+1))=0;
  vwind(find(vwind<fill_in+1))=0;
%  airt(find(airt<fill_in+1))=300; %nan;
  pres(find(pres<fill_in+1))=1.e5; %nan;
%  prc(find(prc<fill_in+1))=0; %nan;

  %Write nc file of SCHISM format (Warning: use single precision!)
  if(i>1 && rem(i,24)==0)

  disp(strcat('Outputting Day: ',num2str(day)));
  ncid2 = netcdf.create(strcat('sflux_air_1_Apr',num2str(day),'.nc'),'CLOBBER');

  dims(1)=netcdf.defDim(ncid2,'nx_grid',n2(1,1));
  dims(2)=netcdf.defDim(ncid2,'ny_grid',n1(1,1));
  ntimes=24;
  dims(3)=netcdf.defDim(ncid2,'time',ntimes);

  % Define new variable in the file.
  timeid = netcdf.defVar(ncid2,'time','float',dims(3));
  lonid = netcdf.defVar(ncid2,'lon','float',dims(1:2));
  latid = netcdf.defVar(ncid2,'lat','float',dims(1:2));
  uid = netcdf.defVar(ncid2,'uwind','float',dims(1:3));
  vid = netcdf.defVar(ncid2,'vwind','float',dims(1:3));
  pid = netcdf.defVar(ncid2,'prmsl','float',dims(1:3)); %pressure
  tid = netcdf.defVar(ncid2,'stmp','float',dims(1:3)); %air temp.
  hid = netcdf.defVar(ncid2,'spfh','float',dims(1:3)); %humidity

  % Leave define mode and enter data mode to write data.
  netcdf.endDef(ncid2)
  
  %Write out new nc file
  times=(0:23)/24.; %hourly
  netcdf.putVar(ncid2,timeid,times);
    % Re-enter define mode.
    netcdf.reDef(ncid2);
    % Create an attribute associated with the variable.
    netcdf.putAtt(ncid2,timeid,'units',strcat('days since 2007-04-',num2str(day)));
    netcdf.putAtt(ncid2,timeid,'base_date',int32([2007 4 day 0])); %must use int32
    netcdf.endDef(ncid2)

    netcdf.putVar(ncid2,lonid,lon_x);
    netcdf.putVar(ncid2,latid,lat_y);
    netcdf.putVar(ncid2,uid,uwind(:,:,i-23:i));
    netcdf.putVar(ncid2,vid,vwind(:,:,i-23:i));
    netcdf.putVar(ncid2,pid,pres(:,:,i-23:i));
    netcdf.putVar(ncid2,tid,airt);
    spfh=0.0175873*ones(n2(1,1),n1(1,1),ntimes);
    netcdf.putVar(ncid2,hid,spfh);
    netcdf.close(ncid2);
  end  
end    
