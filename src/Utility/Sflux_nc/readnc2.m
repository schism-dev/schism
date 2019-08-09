%Read in our sflux netcdf files into matlab and plot; modify as appro.
% e.g., the start and end stack #; time stamp for plot; sub-sampling freq.
% in plot
clear all; close all;

scrsz = get(0,'ScreenSize'); %screen size
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

CB_bnd=load('CB_bnd.xy'); %load domain bnd

%NARR files
fill_in=1.e9; %junk value from nc files
delete('sflux.avi','f');
avi_out = avifile('sflux.avi','FPS',5);
for i=1:10 %stack # for nc files
  char=sprintf('%3.3d',i);
  filen=strcat('sflux_air_1.',char,'.nc');
  ncid0 = netcdf.open(filen,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time_narr= double(netcdf.getVar(ncid0, vid));
  %Get start time
  base=double(netcdf.getAtt(ncid0,vid,'base_date','int32'));
  
  vid=netcdf.inqVarID(ncid0,'lon');
  lon_narr = double(netcdf.getVar(ncid0, vid)); 
  vid=netcdf.inqVarID(ncid0,'lat');
  lat_narr = double(netcdf.getVar(ncid0, vid)); 

  vid=netcdf.inqVarID(ncid0,'uwind');
  %Time dimension is last index in uwind_narr etc.
  uwind_narr= double(netcdf.getVar(ncid0, vid)); 
  uwind_narr(find(abs(uwind_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'vwind');
  vwind_narr= double(netcdf.getVar(ncid0, vid)); 
  vwind_narr(find(abs(vwind_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'prmsl');
  pres_narr= double(netcdf.getVar(ncid0, vid)); 
  pres_narr(find(abs(pres_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'stmp');
  airt_narr= double(netcdf.getVar(ncid0, vid)); %air temp.
  airt_narr(find(abs(airt_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'spfh'); %humidity
  spfh_narr= double(netcdf.getVar(ncid0, vid)); %specific humidity
  spfh_narr(find(abs(spfh_narr)>fill_in))=nan;
  netcdf.close(ncid0);

  disp(strcat('Done reading: ',filen));

  for j=1:length(time_narr)
    %Time stamp for plot
    tmp=datenum(double(base(1)),double(base(2)),double(base(3)))+time_narr(j);
    date_narr=datestr(double(tmp),21); %strcat('Sept ',num2str(i),': ',num2str(time_narr(j)*24),' hour UTC');
    %Compute bound
    xmin=min(min(lon_narr)); 
    ymin=min(min(lat_narr));
    xmax=max(max(lon_narr));
    ymax=max(max(lat_narr));

    %plot domain bnd
    hold on;
    plot(CB_bnd(:,1),CB_bnd(:,2),'k'); 

    %Scalar plot
%    contour(lon_narr,lat_narr,pres_narr(:,:,j));
%    caxis([9.7e4 10.2e4]);
%    colorbar;

    %Subsample
    ind1=1:1:size(lon_narr,1);
    ind2=1:1:size(lon_narr,2);

    %Vel. plot
    scale=0.1; %scale to fit
    quiver(lon_narr(ind1,ind2),lat_narr(ind1,ind2),...
    uwind_narr(ind1,ind2,j)*scale,vwind_narr(ind1,ind2,j)*scale,0,'k');
    quiver(30,47.9,5*scale,0.,0,'r');
    text(30,48.1,'5 m/s');

    %Atmos. pressure
    pcolor(lon_narr(ind1,ind2),lat_narr(ind1,ind2),double(airt_narr(ind1,ind2,j)-273));
    caxis([0 30]);
    shading interp; %to get rid of gridline
    colormap(hsv(40));
    colorbar;

    %Warning: avi does not like changes in xmin etc; use fixed scale
    axis([27 42 41 48]);
    xlabel('Lon'); ylabel('Lat');
    title(['Air Temp (C) @' date_narr]);

%   Stop here for testing
%    return;

    frame = getframe(gcf);
    avi_out=addframe(avi_out,frame);
    clf; %clear figure
  end %j

  clear base time_narr lon_narr lat_narr uwind_narr vwind_narr pres_narr airt_narr spfh_narr;
end %for all nc files
avi_out=close(avi_out);
