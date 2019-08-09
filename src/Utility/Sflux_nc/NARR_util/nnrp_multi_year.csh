#! /bin/csh


######## define the script name to be called, etc ######################
set script = '/home/workspace/ccalmr37/mazulauf/new_data/scripts/nnrp_full_year.csh'


######## set the date variables ########################################
if ( $1 == "" ) then
  echo you must specify the start_year!
  exit
else
  set start_year = $1
endif

if ( $2 == "" ) then
  echo you must specify the stop_year!
  exit
else
  set stop_year = $2
endif


# top of loop over years
set year = $start_year
while ($year <= $stop_year)

  echo
  echo Year : $year

# run script
  $script $year
# echo $script $year

# bottom of loop over years
  set year = `echo \($year + 1\) | bc` >& /dev/null
end

