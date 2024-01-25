#!/bin/csh
set idate=`/bin/date --date='0 days ago' +%Y-%m-%d`
echo $idate
wget -O St_Lawrence_river_discharge_${idate}.csv https://dd.meteo.gc.ca/hydrometric/csv/QC/hourly/QC_02OA016_hourly_hydrometric.csv
