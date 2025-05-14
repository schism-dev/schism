#!/bin/bash

#From Qian Wang (Feb 2025), for downloading HYCOM including ice info
WGET='/usr/bin/wget'
YEAR='1994'
MONTH='05'
DAY='30'
StartSeq='0'
EndSeq='2'

NCSS='http://ncss.hycom.org/thredds/ncss'
MODEL='GLBv0.08'
EXPT='expt_53.X'

VARS="var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v"
VARS1="var=surf_el"
VARS2="var=salinity&var=water_temp"
VARS3="var=water_u&var=water_v"
VARS4="var=sic&var=sih&var=siu&var=siv&var=sss&var=sst&var=ssu&var=ssv"
NORTH='north=90'
SOUTH='south=30'
##EAST='east=179.92'
EAST='east=180'
WEST='west=-180'

for PlusDay in `seq $StartSeq $EndSeq`; do
 
  MyTime=`date -d "$YEAR-$MONTH-$DAY 3 +$PlusDay days" +%Y-%m-%dT%H:%M:%SZ`
  Time="time=$MyTime"

  yyear="`echo $MyTime | cut -d '-' -f 1`"

  #OutFile="`echo $MyTime | cut -d 'T' -f 1`.nc"

  MyTime1=`date -d "$YEAR-$MONTH-$DAY +$PlusDay days" +%s`
  MyTime2=`date -d "1993-12-31" +%s`
  TimeDiff0=`expr $MyTime1 - $MyTime2`
  TimeDiff=`expr $TimeDiff0 / 86400`
  OutFile1="SSH_$TimeDiff.nc"
 
  URL="$NCSS/$MODEL/$EXPT/data/$yyear?$VARS1&$NORTH&$WEST&$EAST&$SOUTH&$Time&addLatLon=true&accept=netcdf4"
 
  if [ -s $OutFile1 ]; then
        echo "[warning] File $OutFile1 exists (skipping)"
  else
        wget -O $OutFile1 "$URL"
  fi

    OutFile2="TS_$TimeDiff.nc"
 
  URL="$NCSS/$MODEL/$EXPT/data/$yyear?$VARS2&$NORTH&$WEST&$EAST&$SOUTH&$Time&addLatLon=true&accept=netcdf4"
 
  if [ -s $OutFile2 ]; then
        echo "[warning] File $OutFile2 exists (skipping)"
  else
        wget -O $OutFile2 "$URL"
  fi

    OutFile3="UV_$TimeDiff.nc"
 
  URL="$NCSS/$MODEL/$EXPT/data/$yyear?$VARS3&$NORTH&$WEST&$EAST&$SOUTH&$Time&addLatLon=true&accept=netcdf4"
 
  if [ -s $OutFile3 ]; then
        echo "[warning] File $OutFile3 exists (skipping)"
  else
        wget -O $OutFile3 "$URL"
  fi
      OutFile4="ICE_$TimeDiff.nc"
 
  URL="$NCSS/$MODEL/$EXPT/data_ice/$yyear?$VARS4&$NORTH&$WEST&$EAST&$SOUTH&$Time&addLatLon=true&accept=netcdf4"
 
  if [ -s $OutFile4 ]; then
        echo "[warning] File $OutFile4 exists (skipping)"
  else
        wget -O $OutFile4 "$URL"
  fi
done

#http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1/2005?var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v&north=24&west=70&east=46&south=2&disableProjSubset=on&horizStride=1&time=2005-12-31T00%3A00%3A00Z&vertCoord=&accept=netcdf
#https://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data_ice/1994?var=sic&var=sih&var=siu&var=siv&var=sss&var=sst&var=ssu&var=ssv&north=90.0000&west=-180.0000&east=180&south=30&disableProjSubset=on&horizStride=1&time=1995-01-01T10%3A59%3A31.200Z&accept=netcdf4
