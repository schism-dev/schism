%Generate SCHISM compliant nc files using inverse weight interpolaton
%from obs points
clear all; close all;
delete('*.nc');
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)


%Hourly wind data (from t=0)
sta=load('sta_wind_ll.bp.mlab'); %#,lon,lat
nsta=size(sta,1);
wind_sta=load('sta_wind.in'); %hourly (u,v) @ stations 1-12
pres_sta=load('sta_pres.in'); %pres @ stations 1-12
nstep=size(wind_sta,1);

%Output ll grid
n1=51; n2=52;
xmin=-77.467285; xmax=-74.637223; dx=(xmax-xmin)/(n1-1);
ymin=36.298963; ymax=39.744966; dy=(ymax-ymin)/(n2-1);
%[lon,lat]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
%[n1 n2]=size(lon);
for i=1:n1
  for j=1:n2
    lon(i,j)=xmin+dx*(i-1);
    lat(i,j)=ymin+dy*(j-1);
  end
end

%CB_bnd=load('CB_bnd.xy'); %load domain bnd; nd,lon,lat,dp
%hold on;
%plot(lon,lat,'r.');
%plot(CB_bnd(:,2),CB_bnd(:,3),'k.');
%return;

%Calculate weights
for i=1:nsta
  dis=sqrt((lon-sta(i,2)).^2+(lat-sta(i,3)).^2);
  dis(find(dis)<1.e-8)=1.e-8;
  dis_out(i,:,:)=dis;
end %for

for ix=1:n1
  for iy=1:n2
    tmp=sum(1./dis_out(:,ix,iy));
    for i=1:nsta
      tmp1(:,i)=pres_sta(:,i)/dis_out(i,ix,iy);
      tmp2(:,i)=wind_sta(:,2*i-1)/dis_out(i,ix,iy);
      tmp3(:,i)=wind_sta(:,2*i)/dis_out(i,ix,iy);
    end %for i
    pres(ix,iy,:)=sum(tmp1,2)/tmp;
    uwind(ix,iy,:)=sum(tmp2,2)/tmp;
    vwind(ix,iy,:)=sum(tmp3,2)/tmp;
  end %for iy
end %for ix

%Plot wind
%scale=0.02; %scale to fit
%rec=364*24+1;
%hold on;
%CB_bnd=load('CB_bnd.xy'); %load domain bnd; nd,lon,lat,dp
%plot(CB_bnd(:,2),CB_bnd(:,3),'k.',sta(:,2),sta(:,3),'go');
%quiver(sta(:,2),sta(:,3),wind_sta(rec,1:2:2*nsta)'*scale,...
%       wind_sta(rec,2:2:2*nsta)'*scale,0,'b');
%quiver(lon,lat,uwind(:,:,rec)*scale,vwind(:,:,rec)*scale,0,'c');
%quiver(-77.9,39.9,5*scale,0.,0,'r');
%text(-77.9,39.8,'5 m/s');
%return;

%Check
disp('Max of pres, uwind, vwind:');
[max(max(max(pres))) max(max(max(uwind))) max(max(max(vwind)))]
disp('Min of pres, uwind, vwind:');
[min(min(min(pres))) min(min(min(uwind))) min(min(min(vwind)))]

fill_in=-999; %junk value
for day=1:366 %2003
  %Write nc file of SCHISM format (Warning: use single precision!)
  %Write out new nc file (1 day per file)
  hours=(day-1)*24:(day-1)*24+23; %hourly from Jan. 1
  ntimes=length(hours);
%  disp(strcat('Outputting Day: ',num2str(day)));
  day_str=sprintf('%4.4d',day);
  filenm=strcat('sflux_air_2.',day_str,'.nc');

%IMPORTANT: use nc4 to avoid file size limit!!
%cmode = netcdf.getConstant('NETCDF4');
%cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
%ncid2 = netcdf.create(filenm,cmode);

  ncid2 = netcdf.create(filenm,'CLOBBER');
  dims(1)=netcdf.defDim(ncid2,'n1',n1);
  dims(2)=netcdf.defDim(ncid2,'n2',n2);
  dims(3)=netcdf.defDim(ncid2,'time',ntimes);

  % Define new variable in the file.
  timeid = netcdf.defVar(ncid2,'time','float',dims(3));
  latid = netcdf.defVar(ncid2,'lat','float',dims(1:2));
  lonid = netcdf.defVar(ncid2,'lon','float',dims(1:2));
  uid = netcdf.defVar(ncid2,'uwind','float',dims(1:3));
  vid = netcdf.defVar(ncid2,'vwind','float',dims(1:3));
  pid = netcdf.defVar(ncid2,'prmsl','float',dims(1:3)); %pressure
  tid = netcdf.defVar(ncid2,'stmp','float',dims(1:3)); %air temp.
  hid = netcdf.defVar(ncid2,'spfh','float',dims(1:3)); %humidity

  % Leave define mode and enter data mode to write data.
  netcdf.endDef(ncid2)

  times=hours/24.; %time from Jan. 1
  netcdf.putVar(ncid2,timeid,times);
  % Re-enter define mode.
  netcdf.reDef(ncid2);
  % Create an attribute associated with the variable.
  netcdf.putAtt(ncid2,timeid,'units','days since 2003-01-01'); %This is not needed really
  netcdf.putAtt(ncid2,timeid,'base_date',int32([2003 1 1 0])); %must use int32
  netcdf.endDef(ncid2)

  %Good idea to add a check to make sure the quad is counter-clockwise
  netcdf.putVar(ncid2,latid,lat);
  netcdf.putVar(ncid2,lonid,lon);
  netcdf.putVar(ncid2,uid,uwind(:,:,hours+1));
  netcdf.putVar(ncid2,vid,vwind(:,:,hours+1));
  netcdf.putVar(ncid2,pid,pres(:,:,hours+1));
  airt=300*ones(n1,n2,ntimes); %no air T info so we use a constant 300K
  spfh=0.0175873*ones(n1,n2,ntimes); %use a const specific humidty
  netcdf.putVar(ncid2,tid,airt);
  netcdf.putVar(ncid2,hid,spfh);
  netcdf.close(ncid2);
end %for day
