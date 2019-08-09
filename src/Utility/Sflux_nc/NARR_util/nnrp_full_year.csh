#! /bin/csh


######## define the script name to be called, etc ######################
set script = '/home/workspace/ccalmr37/mazulauf/new_data/scripts/get_nc_nnrp.01.csh'

set util_bin = /home/workspace/ccalmr/forecasts/bin/atmos_nc/bin/
set cjd_util       = $util_bin'cjd'


######## set the date variables ########################################
if ( $1 == "" ) then
  set year = (`/bin/date -u +%Y`)
  echo you must specify the year!
  exit
else
  set year = $1
endif

set cjd_start = (`$cjd_util $year 01 01`)
set cjd_stop  = (`$cjd_util $year 12 31`)

# for testing. . .
# set cjd_stop = `echo \($cjd_start + 1\) | bc` >& /dev/null


# top of loop over Julian days
set cjd = $cjd_start
while ($cjd <= $cjd_stop)

# determine calendar date
  set year  = (`$cjd_util -r $cjd | awk '{print $1}'`)
  set month = (`$cjd_util -r $cjd | awk '{print $2}'`)
  set day   = (`$cjd_util -r $cjd | awk '{print $3}'`)
  set date_string = $year'_'$month'_'$day

  echo
  echo Julian date and date_string : $cjd $date_string

# run script
  $script $year $month $day


# bottom of loop over Julian days
  set cjd = `echo \($cjd + 1\) | bc` >& /dev/null
end

