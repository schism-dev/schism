#! /bin/csh

#####Usage: grab_data.csh <start_year> <stop_year>; may need to update the url below

if ( $1 == "" ) then
  echo enter start_year!
  exit
else
  set start_year = $1
endif

if ( $2 == "" ) then
  echo enter stop_year!
  exit
else
  set stop_year = $2
endif

set year = $start_year

while ($year <= $stop_year)

  set start_month = 2
  set stop_month = 5

  set month = $start_month
  while ($month <= $stop_month)

    if ($month >= 10) then
      set date = $year$month
    else
      set date = $year'0'$month
    endif

####    set url = 'http://dss.ucar.edu/download/Yinglong/'$date'.tar'
    set url = 'http://rda.ucar.edu/download/chifan/PaulT/'$date'.tar'
    echo $date $url
    wget $url

# bottom of loop of months for file retrieval
  @ next_month = $month + 1
  set month = $next_month
end



# bottom of loop of years for file retrieval
  @ next_year = $year + 1
  set year = $next_year
end

exit
