#!/bin/bash
##lookback (days)
#usace
wget -O BON.Flow-Out.Ave.1Hour.txt --no-check-certificate "https://www.nwd-wc.usace.army.mil/dd/common/web_service/webexec/csv?id=BON.Flow-Out.Ave.1Hour.1Hour.CBT-REV:units=kcfs&lookback=10&timeformat=%d-%b-%Y%20%H:%M"

#usgs stations
start_date=`date --date="3 day ago" +%F`
echo $start_date
end_date=$(date "+%F")
echo $end_date
wget -O usgs_flow_14211720.csv --no-check-certificate "https://waterservices.usgs.gov/nwis/iv/?sites=14211720&parameterCd=00060&startDT=${start_date}T08:00:00-07:00&endDT=${end_date}T08:00:00-07:00&siteStatus=all&format=rdb"
wget -O usgs_flow_14220500.csv  --no-check-certificate "https://waterservices.usgs.gov/nwis/iv/?sites=14220500&parameterCd=00060&startDT=${start_date}T08:00:00-07:00&endDT=${end_date}T08:00:00-07:00&siteStatus=all&format=rdb"
wget -O usgs_flow_14243000.csv --no-check-certificate "https://waterservices.usgs.gov/nwis/iv/?sites=14243000&parameterCd=00060&startDT=${start_date}T08:00:00-07:00&endDT=${end_date}T08:00:00-07:00&siteStatus=all&format=rdb"

#fraser river
wget -O fraser_river_08MF005_hourly.csv --no-check-certificate "https://dd.meteo.gc.ca/hydrometric/csv/BC/hourly/BC_08MF005_hourly_hydrometric.csv"
