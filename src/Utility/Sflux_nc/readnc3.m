%Read in an existing netcdf files, and modify and output SCHISM nc input; need to 
% add file names below; modify other parts as appropriate (e.g. start time stamp)
clear all; close all;
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)

CB_bnd=load('CB_bnd.xy'); %load domain bnd

%WRF files
%Each nc file is for a single time
fill_in=-999; %junk value
filenm={}; %define file names here; e.g. WRF
for i=1:length(filenm)
  day=floor((i-1)/24+1); 
  hour=i-1-(day-1)*24; %0-23 hours
  ncid = netcdf.open(filenm{i},'NC_NOWRITE');
  vid1=netcdf.inqVarID(ncid,'g3_lat_0'); %input var./array name
  lat = netcdf.getVar(ncid, vid1); 
  [n1,n2]=size(lat); %reversed from C
  vid2=netcdf.inqVarID(ncid,'g3_lon_1');
  lon = netcdf.getVar(ncid, vid2); 

  vid=netcdf.inqVarID(ncid,'U_GRD_GDS3_HTGL');
  %Note the loc for time dimension (last)
  uwind(:,:,hour+1) = netcdf.getVar(ncid, vid); 

  vid=netcdf.inqVarID(ncid,'V_GRD_GDS3_HTGL');
  vwind(:,:,hour+1) = netcdf.getVar(ncid, vid); 

%  if(i==1)
%    %Debug
%    for l=10:10:n2
%      for m=10:10:n1
%        if(uwind(1,m,l)>fill_in+1)
%          fprintf(fid,'%f %f %f %f\n',lon(m,l),lat(m,l),uwind(1,m,l),vwind(1,m,l));
%        end
%      end
%    end
%  end

  vid=netcdf.inqVarID(ncid,'TMP_GDS3_SFC'); %air temp. in K
  airt(:,:,hour+1) = netcdf.getVar(ncid, vid); 

  vid=netcdf.inqVarID(ncid,'PRMSL_GDS3_MSL'); %air pressure at MSL (Pa)
  pres(:,:,hour+1) = netcdf.getVar(ncid, vid); 

  vid=netcdf.inqVarID(ncid,'PRATE_GDS3_SFC'); %precip. rate (kg/m^2/s)
  prc = netcdf.getVar(ncid, vid); 

%  vid=netcdf.inqVarID(ncid,'LAND_GDS3_SFC'); %Land-sea mask (land=1;sea=0)
%  mask= netcdf.getVar(ncid, vid); 

  disp(strcat('Done reading ',filenm{i}));

  %Deal with junks
  uwind(find(uwind<fill_in+1))=0;
  vwind(find(vwind<fill_in+1))=0;
  airt(find(airt<fill_in+1))=300; %nan;
  pres(find(pres<fill_in+1))=1.e5; %nan;
  prc(find(prc<fill_in+1))=0; %nan;

  %Write nc file of SCHISM format (Warning: use single precision!)
  if(i>1 && rem(i,24)==0)
    disp(strcat('Outputting Day: ',num2str(day)));
    ncid2 = netcdf.create(strcat('sflux_air_2.Sept',num2str(day+15),'.nc'),'CLOBBER');

    dims(1)=netcdf.defDim(ncid2,'n1',n1);
    dims(2)=netcdf.defDim(ncid2,'n2',n2);
    ntimes=24;
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
    times=(0:23)/24.; %hourly
    netcdf.putVar(ncid2,timeid,times);
    % Re-enter define mode.
    netcdf.reDef(ncid2);
    % Create an attribute associated with the variable.
    netcdf.putAtt(ncid2,timeid,'units',strcat('days since 2003-09-',num2str(day+15))); %This is not needed really
    netcdf.putAtt(ncid2,timeid,'base_date',int32([2003 9 day+15 0])); %must use int32
    netcdf.endDef(ncid2)

    %Good idea to add a check to make sure the quad is counter-clockwise
    netcdf.putVar(ncid2,latid,lat);
    netcdf.putVar(ncid2,lonid,lon);
    netcdf.putVar(ncid2,uid,uwind);
    netcdf.putVar(ncid2,vid,vwind);
    netcdf.putVar(ncid2,pid,pres);
    netcdf.putVar(ncid2,tid,airt);
    spfh=0.0175873*ones(n1,n2,ntimes);
    netcdf.putVar(ncid2,hid,spfh);
    netcdf.close(ncid2);
  end %open and write daily nc file
end %for i - all nc files
