#!/bin/bash
##lookback (days)
wget -O BON.Temp-Water.Inst.1Hour.txt --no-check-certificate "https://www.nwd-wc.usace.army.mil/dd/common/web_service/webexec/csv?id=BON.Temp-Water.Inst.1Hour.0.GOES-REV:units=F&lookback=10&timeformat=%d-%b-%Y%20%H:%M"

wget -O WRNO.Temp-Water.Inst.1Hour.txt --no-check-certificate "https://www.nwd-wc.usace.army.mil/dd/common/web_service/webexec/csv?id=WRNO.Temp-Water.Inst.1Hour.0.GOES-REV:units=F&lookback=10&timeformat=%d-%b-%Y%20%H:%M"

#usgs stations
start_date=`date --date="3 day ago" +%F`
echo $start_date
end_date=$(date "+%F")
echo $end_date
wget -O usgs_temp_14211720.csv --no-check-certificate "https://waterservices.usgs.gov/nwis/iv/?sites=14211720&parameterCd=00010&startDT=${start_date}T08:00:00-07:00&endDT=${end_date}T08:00:00-07:00&siteStatus=all&format=rdb"
