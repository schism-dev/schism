%Purpose: convert DWD wind (nc) to sflux for SCHISM
%Author: Johannes Pein (johannes.pein@uni-oldenburg.de)
%Date: Nov 2012
clear all; close all;
for d=1:31
  if d<10
    filename=['dwd2011090' num2str(d) '.nc']
  else 
    filename=['dwd201109' num2str(d) '.nc']
  end
%  cd /h/tmp/ksd/dwd
  ncid1=netcdf.open(['/h/tmp/ksd/dwd/' filename],'NC_NOWRITE');
  time=netcdf.getVar(ncid1,0); %int32!!
  lon=netcdf.getVar(ncid1,1);
  lat=netcdf.getVar(ncid1,2);
  press=netcdf.getVar(ncid1,4);
  u10=netcdf.getVar(ncid1,7);
  v10=netcdf.getVar(ncid1,8);
  t2= netcdf.getVar(ncid1,5); %air temp.
  dev2= netcdf.getVar(ncid1,6); %Dew point temp.
  netcdf.close(ncid1);

  ntime=length(time);
  dt=time(2)-time(1); %in sec
  [X,Y] = meshgrid(lon,lat);
  %define new netcdf file
%  cd /data/pein/wind
  ncid2=netcdf.create(['/data/zhangj/wind/' filename(1:11) '_4selfe.nc'],'NC_CLOBBER');
  lon_dim=netcdf.defDim(ncid2,'nx_grid',length(lon));
  lat_dim=netcdf.defDim(ncid2,'ny_grid',length(lat));
  time_dim=netcdf.defDim(ncid2,'time',ntime);

  time_var=netcdf.defVar(ncid2,'time','float',time_dim);
  lon_var=netcdf.defVar(ncid2,'lon','float',[lat_dim lon_dim]);
  lat_var=netcdf.defVar(ncid2,'lat','float',[lat_dim lon_dim]);
  press_var=netcdf.defVar(ncid2,'prmsl','float',[lat_dim lon_dim time_dim]);
  u10_var=netcdf.defVar(ncid2,'uwind','float',[lat_dim lon_dim time_dim]);
  v10_var=netcdf.defVar(ncid2,'vwind','float',[lat_dim lon_dim time_dim]);
  stmp_var=netcdf.defVar(ncid2,'stmp','float',[lat_dim lon_dim time_dim]);
  spfh_var=netcdf.defVar(ncid2,'spfh','float',[lat_dim lon_dim time_dim]);

  %und jetzt noch attribute
  netcdf.putAtt(ncid2,time_var,'long_name','Time');
  netcdf.putAtt(ncid2,time_var,'standard_name','time');
  netcdf.putAtt(ncid2,time_var,'units',['days since ' num2str(filename(4:7)) '-' num2str(filename(8:9)) '-' num2str(filename(10:11)) ]);
  a=str2double(filename(4:7));
  b=str2double(filename(8:9));
  c=str2double(filename(10:11));
  netcdf.putAtt(ncid2,time_var,'base_date',int32([a b c 0]));
  
  netcdf.putAtt(ncid2,lon_var,'long_name','Longitude');
  netcdf.putAtt(ncid2,lon_var,'standard_name','longitude');
  netcdf.putAtt(ncid2,lon_var,'units','degrees_east'); 

  netcdf.putAtt(ncid2,lat_var,'long_name','Latitude');
  netcdf.putAtt(ncid2,lat_var,'standard_name','latitude');
  netcdf.putAtt(ncid2,lat_var,'units','degrees_north'); 

  netcdf.putAtt(ncid2,press_var,'long_name','Pressure reduced to MSL');
  netcdf.putAtt(ncid2,press_var,'standard_name','air_pressure_at_sea_level');
  netcdf.putAtt(ncid2,press_var,'units','Pa');

  netcdf.putAtt(ncid2,u10_var,'long_name','Surface Eastward Air Velocity');
  netcdf.putAtt(ncid2,u10_var,'standard_name','eastward_wind');
  netcdf.putAtt(ncid2,u10_var,'units','m/s');

  netcdf.putAtt(ncid2,v10_var,'long_name','Surface Northward Air Velocity');
  netcdf.putAtt(ncid2,v10_var,'standard_name','northward_wind');
  netcdf.putAtt(ncid2,v10_var,'units','m/s');

  netcdf.putAtt(ncid2,stmp_var,'long_name','Surface Air Temperature (2m AGL)');
  netcdf.putAtt(ncid2,stmp_var,'standard_name','air_temperature');
  netcdf.putAtt(ncid2,stmp_var,'units','Celsius');

  netcdf.putAtt(ncid2,stmp_var,'long_name','Surface Specific Humidity (2m AGL)');
  netcdf.putAtt(ncid2,stmp_var,'standard_name','specific_humidity');
  netcdf.putAtt(ncid2,stmp_var,'units','1');

  netcdf.endDef(ncid2);
  %end define mode
  %write data 
%  ncid2=netcdf.open(['/data/zhangj/wind/' filename(1:11) '_4selfe.nc'],'NC_WRITE');
  netcdf.putVar(ncid2,lon_var,X);
  netcdf.putVar(ncid2,lat_var,Y);
 
  tcount=0;
  time2=0.;
  for sl=1:ntime
    start=[0 0 tcount];
    count=[length(lat) length(lon) 1];
    press2=rot90(squeeze(press(:,:,sl)));
    netcdf.putVar(ncid2,press_var,start,count,press2);
    u10_2=rot90(squeeze(u10(:,:,sl)));
    netcdf.putVar(ncid2,u10_var,start,count,u10_2);
    v10_2=rot90(squeeze(v10(:,:,sl)));
    netcdf.putVar(ncid2,v10_var,start,count,v10_2);
    t2_2=rot90(squeeze(t2(:,:,sl)));
    netcdf.putVar(ncid2,stmp_var,start,count,t2_2);
%    dev2_2=rot90(squeeze(dev2(:,:,sl)));
    spfh=single(1.e-2*ones(size(t2_2)));
    netcdf.putVar(ncid2,spfh_var,start,count,spfh);
    netcdf.putVar(ncid2,time_var,tcount,1,single(time2/86400.));

    tcount=tcount+1;
    time2=time2+single(dt);
  end %sl

  netcdf.close(ncid2);
end %d
