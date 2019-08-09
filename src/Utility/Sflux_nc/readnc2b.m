%Same as readnc2 but also output time series at given points (nearest neighbor interp)
clear all; close all;
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)

%Input point loc
xin=[30.71574220 35.93387997 41.61301274 37.34236566 39.16675527 31.82269307 30.31455093 29.62052320];
yin=[41.12126787 41.76327275 42.19142915 45.33307922 47.21103280 46.41283803 45.94598945 44.83146673];
npt=length(xin);

%nc files
fill_in=1.e9; %junk value from nc files
iout=0; %counter for output
for i=4:11 %stack # for nc files
  char=sprintf('%2.2d',i);
  filen=strcat('sflux_air_2013',char,'.nc');
  ncid0 = netcdf.open(filen,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time_narr= double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'lon');
  lon_narr = double(netcdf.getVar(ncid0, vid)); 
  vid=netcdf.inqVarID(ncid0,'lat');
  lat_narr = double(netcdf.getVar(ncid0, vid)); 

  vid=netcdf.inqVarID(ncid0,'stmp');
  %Time dimension is last index in airt etc.
  airt= double(netcdf.getVar(ncid0, vid)); 
  %airt(find(abs(airt)>fill_in))=nan;

  disp(strcat('Done reading: ',filen));

  %Find nearest node for interpolation
  for k=1:npt
    dist2=(lon_narr-xin(k)).^2+(lat_narr-yin(k)).^2;
    [tmp,irow0]=min(dist2);
    [tmp2,icol(k)]=min(tmp);
    irow(k)=irow0(icol(k)); 
    if(iout==0)
      disp('Nearest pt='); [lon_narr(irow(k),icol(k)) lat_narr(irow(k),icol(k))]
    end
  end %for k

  %Outputs
  base=datenum(2013,i,1)-datenum(2013,4,1);
  ntime=length(time_narr);
  time_out(iout+1:iout+ntime)=time_narr(:)+base;
  for k=1:npt
    var_out(iout+1:iout+ntime,k)=airt(irow(k),icol(k),1:ntime)-273;
  end %for k
  iout=iout+ntime;

  clear time_narr lon_narr lat_narr airt dist2 irow icol tmp tmp2;
end %for all nc files

fid=fopen('readnc2b.out','w');
fprintf(fid,'%f %f %f %f %f %f %f %f %f\n',[time_out*86400; var_out']);
fclose(fid);

