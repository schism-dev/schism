%A simple aproach to download HYCOM
%Go to HYCOM site and netcdf subsetting service and fill out the query. You can see the
%URL from 'NCSS Request URL'. Use this url as:
%url=... (download SSH, TS, UV separately - include lon/lat)
%urlwrite(url,'out.nc');
%To avoid timeout/large file size, try downloading several days at a time.

close all; clear all;

url='http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1?var=surf_el&north=-15.25&west=157.58&east=171.60&south=-25.88&disableLLSubset=on&disableProjSubset=on&horizStride=1&time_start=2000-01-01T00%3A00%3A00Z&time_end=2000-02-01T00%3A00%3A00Z&timeStride=1&vertCoord=&addLatLon=true&accept=netcdf';
urlwrite(url,'SSH_1.nc');

url='http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1?var=salinity&var=water_temp&north=-15.25&west=157.58&east=171.60&south=-25.88&disableLLSubset=on&disableProjSubset=on&horizStride=1&time_start=2000-01-01T00%3A00%3A00Z&time_end=2000-02-01T00%3A00%3A00Z&timeStride=1&vertCoord=&accept=netcdf';
urlwrite(url,'TS_1.nc');

url='http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1?var=water_u&var=water_v&north=-15.25&west=157.58&east=171.60&south=-25.88&disableLLSubset=on&disableProjSubset=on&horizStride=1&time_start=2000-01-01T00%3A00%3A00Z&time_end=2000-02-01T00%3A00%3A00Z&timeStride=1&vertCoord=&accept=netcdf';
urlwrite(url,'UV_1.nc');
