%Read in dataset "1" (larger grid) and interpolate onto grid "2" (smaller)
%This is useful when the dataset "2" is missing some days
%No interpolation in time (as SCHISM can handle variable # of timesteps in each file)
%Assume dataset "1" and setnm2 are already in SCHISM format
%Run inside sflux/; outputs: sflux_air_2*.nc

%Dimension for arrays reversed from ncdump (FORTRAN convention)
clear all; close all;
%Input names; also need to manually update "units" and "base_date" below
setnm1='sflux_air_1.'; %more coming later
nfiles=24; %# of files from "1"
setnm2='sflux_air_2.Sept16.nc'; %used to get grid 2 only
outnm='sflux_air_2.'; %more later

CB_bnd=load('CB_bnd.xy'); %load domain bnd

%Read grid "2" first
ncid = netcdf.open(setnm2,'NC_NOWRITE');
vid=netcdf.inqVarID(ncid,'lon');
lon= netcdf.getVar(ncid, vid);
vid=netcdf.inqVarID(ncid,'lat');
lat= netcdf.getVar(ncid, vid);
netcdf.close(ncid);
[nxout,nyout]=size(lon);
%Cast into column vector for interpolation
lon2=reshape(double(lon),nxout*nyout,1);
lat2=reshape(double(lat),nxout*nyout,1);

%Read dataset "1"
fill_in=1.e9; %junk value from nc files
avi_out = avifile('out.avi');
for i=1:nfiles 
  char=sprintf('%4.4d',i);
  filen=strcat(setnm1,char,'.nc');
  ncid0 = netcdf.open(filen,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time_narr= netcdf.getVar(ncid0, vid);
  vid=netcdf.inqVarID(ncid0,'lon');
  lon_narr = netcdf.getVar(ncid0, vid); 
  [n1,n2]=size(lon_narr);
  vid=netcdf.inqVarID(ncid0,'lat');
  lat_narr = netcdf.getVar(ncid0, vid); 

  vid=netcdf.inqVarID(ncid0,'uwind');
  %Time dimension is last index in uwind_narr etc.
  uwind_narr= netcdf.getVar(ncid0, vid); 
  vid=netcdf.inqVarID(ncid0,'vwind');
  vwind_narr= netcdf.getVar(ncid0, vid); 
  vid=netcdf.inqVarID(ncid0,'prmsl');
  pres_narr= netcdf.getVar(ncid0, vid); 
  vid=netcdf.inqVarID(ncid0,'stmp');
  airt_narr= netcdf.getVar(ncid0, vid); 
  vid=netcdf.inqVarID(ncid0,'spfh'); %humidity
  spfh_narr= netcdf.getVar(ncid0, vid); 

%  uwind_narr(find(abs(uwind_narr)>fill_in))=nan;
%  vwind_narr(find(abs(vwind_narr)>fill_in))=nan;
%  pres_narr(find(abs(pres_narr)>fill_in))=nan;
%  airt_narr(find(abs(airt_narr)>fill_in))=nan;
%  spfh_narr(find(abs(spfh_narr)>fill_in))=nan;

  disp(strcat('Done reading: ',filen));

  %Do interpolation in space
  for j=1:length(time_narr)
    %Time stamp
    date_narr=strcat('Sept ',num2str(i+7),': ',num2str(time_narr(j)*24),' hour UTC');

    %TriScatteredInterp: use Delauney triangular of convex hull
    %Construct interpolant first
    fun=TriScatteredInterp(reshape(double(lon_narr),n1*n2,1),...
        reshape(double(lat_narr),n1*n2,1),reshape(double(uwind_narr(:,:,j)),n1*n2,1));    
    out=fun(lon2,lat2); uwind_out(:,:,j)=reshape(single(out),nxout,nyout);
    if(sum(isnan(out)) ~=0); error('Interp. failed(1)'); end;

    fun=TriScatteredInterp(reshape(double(lon_narr),n1*n2,1),...
        reshape(double(lat_narr),n1*n2,1),reshape(double(vwind_narr(:,:,j)),n1*n2,1));    
    out=fun(lon2,lat2); vwind_out(:,:,j)=reshape(single(out),nxout,nyout);

    fun=TriScatteredInterp(reshape(double(lon_narr),n1*n2,1),...
        reshape(double(lat_narr),n1*n2,1),reshape(double(pres_narr(:,:,j)),n1*n2,1));    
    out=fun(lon2,lat2); pres_out(:,:,j)=reshape(single(out),nxout,nyout);

    fun=TriScatteredInterp(reshape(double(lon_narr),n1*n2,1),...
        reshape(double(lat_narr),n1*n2,1),reshape(double(airt_narr(:,:,j)),n1*n2,1));    
    out=fun(lon2,lat2); airt_out(:,:,j)=reshape(single(out),nxout,nyout);

    fun=TriScatteredInterp(reshape(double(lon_narr),n1*n2,1),...
        reshape(double(lat_narr),n1*n2,1),reshape(double(spfh_narr(:,:,j)),n1*n2,1));    
    out=fun(lon2,lat2); spfh_out(:,:,j)=reshape(single(out),nxout,nyout);

    %Plot
    if(2==1 && i>=9)
      %Compute bound
%      xmin=min(min(lon_narr)); 
%      ymin=min(min(lat_narr));
%      xmax=max(max(lon_narr));
%      ymax=max(max(lat_narr));
      hold on;
      plot(CB_bnd(:,1),CB_bnd(:,2),'k'); %plot bnd

      %Scalar plot
%      contour(lon_narr,lat_narr,pres_narr(:,:,j));
%      caxis([9.7e4 10.2e4]);
%      colorbar;


      %Vel. plot
      scale=0.05; %scale to fit
      %Subsample
      ind1=10:10:size(lon,1);
      ind2=10:10:size(lon,2);
      quiver(lon(ind1,ind2),lat(ind1,ind2),uwind_out(ind1,ind2,j)*scale,vwind_out(ind1,ind2,j)*scale,0,'b');
%      ind1=2:2:size(lon_narr,1);
%      ind2=2:2:size(lon_narr,2);
%      quiver(lon_narr(ind1,ind2),lat_narr(ind1,ind2),uwind_narr(ind1,ind2,j)*scale,vwind_narr(ind1,ind2,j)*scale,0,'g');
      quiver(-77.9,39.9,5*scale,0.,0,'r');
      text(-77.9,39.8,'5 m/s');
      text(-76.9,39.8,date_narr);
      %Warning: avi does not like changing in xmin etc; use fixed scale
%      axis([xmin-0.1 xmax+0.1 ymin-0.1 ymax+0.1]);
      axis([-81 -72.6 33.32 40.45]);
      xlabel('Lon'); ylabel('Lat');

      frame = getframe(gca);
      avi_out=addframe(avi_out,frame);
      clf; %clear figure
    end %plot
  end %j - time steps

  %Output new nc files
  ncid2 = netcdf.create(strcat(outnm,char,'.nc'),'CLOBBER');

  dims(1)=netcdf.defDim(ncid2,'n1',nxout);
  dims(2)=netcdf.defDim(ncid2,'n2',nyout);
  ntimes=length(time_narr);
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

  %Write out new nc file
  netcdf.putVar(ncid2,timeid,time_narr);
  % Re-enter define mode.
  netcdf.reDef(ncid2);
  % Create an attribute associated with the variable.
  netcdf.putAtt(ncid2,timeid,'units',strcat('days since 2003-09-',num2str(i+7))); %This is not needed really
  netcdf.putAtt(ncid2,timeid,'base_date',int32([2003 9 i+7 0])); %must use int32
  netcdf.endDef(ncid2)

  %Good idea to add a check to make sure the quad is counter-clockwise
  netcdf.putVar(ncid2,latid,lat);
  netcdf.putVar(ncid2,lonid,lon);
  netcdf.putVar(ncid2,uid,uwind_out);
  netcdf.putVar(ncid2,vid,vwind_out);
  netcdf.putVar(ncid2,pid,pres_out);
  netcdf.putVar(ncid2,tid,airt_out);
  netcdf.putVar(ncid2,hid,spfh_out);
  netcdf.close(ncid2);
end %for all nc files
avi_out=close(avi_out);
