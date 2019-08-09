#! /bin/csh


######## define the model name, etc ####################################
set model = 'nam'
set grid_name = nam_grid_218


######## define specifics for rotation of wind fields ##################
set proj = 2
set scale = 1.0
set true_lat_1 = 25.0
set true_lat_2 = 25.0
set true_lon = -95.0


######## define specifics for domain reduction step ####################
set lon_min = -180
set lon_max = -100
set lat_min = -90
set lat_max = 90


######## set the date variables ########################################
if ( $1 == "" ) then
  set year = (`/bin/date -u +%Y`)
else
  set year = $1
endif

if ( $2 == "" ) then
  set month = (`/bin/date -u +%m`)
else
  set month = $2
endif

if ( $3 == "" ) then
  set day = (`/bin/date -u +%d`)
else
  set day = $3
endif

set date_string = $year'_'$month'_'$day
echo 'Date:  ' $date_string


######## define parent directory scripts, data, etc ####################
if ( $4 == "" ) then
  set script_parent = /home/workspace/ccalmr/forecasts/bin/atmos_nc/
else
  set script_parent = $4
endif

set grid_source = $script_parent'data/'$grid_name

echo script_parent:
echo     $script_parent

echo grid_source:
echo     $grid_source
echo " "


######## define parent directory for data destination, etc #############
if ( $5 == "" ) then
  set dst_parent = '/home/workspace/ccalmr37/sflux_data/'$model'/'
else
  set dst_parent = $5
endif

set work_dir   = $dst_parent'work/'

echo dst_parent:
echo     $dst_parent

echo work_dir:
echo     $work_dir
echo " "


######## define locations and names for utilities ######################
set util_bin       = $script_parent'bin/'
set wgrib          = $util_bin'wgrib'
set asc2netcdf     = $util_bin'asc2netcdf.01'
set nc_rotate_wind = $util_bin'nc_rotate_wind.01'
set nc_restrict    = $util_bin'nc_restrict.01'
set nc_fix_avg     = $util_bin'nc_fix_avg.01'
set nc_make_prate  = $util_bin'nc_make_prate.01'


######## define the URL of the data directory ##########################
set url_dir = 'ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.'$year$month$day'/'


######## create directories, etc #######################################
set dst_dir   = $dst_parent$year'_'$month'/'
mkdir -p $dst_dir
mkdir -p $work_dir
cd $work_dir


######## remove any preexisting files that we don't want/need ##########
\rm -f $work_dir*.txt $work_dir*.nc $work_dir*.deck


######## create links to data, etc #####################################
ln -sf $grid_source'.lon.txt' ./lon.txt
ln -sf $grid_source'.lat.txt' ./lat.txt


######## define the netcdf file names (total and local) ################
set air_file_total = $model'_air.total.'$date_string'.nc'
set air_file_local = $model'_air.local.'$date_string'.nc'
set rad_file_total = $model'_rad.total.'$date_string'.nc'
set rad_file_local = $model'_rad.local.'$date_string'.nc'
set prc_file_total = $model'_prc.total.'$date_string'.nc'
set prc_file_local = $model'_prc.local.'$date_string'.nc'
set apcp_file_total = $model'_apcp.total.'$date_string'.nc'


######## create total files (write lon and lat) #######################

# write lon to air_file_total
cat > asc2netcdf.deck << EOF
&inputs_list
in_data             = 'lon.txt',
out_file            = '$air_file_total',
glob_attr_name(1)   = 'Conventions',
glob_attr_type(1)   = 'char',
glob_attr_value(1)  = 'CF-1.0',
data_name           = 'lon',
time_var            = .false.,
year                = $year,
month               = $month,
day                 = $day,
data_attr_name(1)   = 'long_name',
data_attr_type(1)   = 'char',
data_attr_value(1)  = 'Longitude',
data_attr_name(2)   = 'standard_name',
data_attr_type(2)   = 'char',
data_attr_value(2)  = 'longitude',
data_attr_name(3)   = 'units',
data_attr_type(3)   = 'char',
data_attr_value(3)  = 'degrees_east',
/
EOF
$asc2netcdf asc2netcdf.deck

# write lat to air_file_total
cat > asc2netcdf.deck << EOF
&inputs_list
in_data             = 'lat.txt',
out_file            = '$air_file_total',
glob_attr_name(1)   = 'Conventions',
glob_attr_type(1)   = 'char',
glob_attr_value(1)  = 'CF-1.0',
data_name           = 'lat',
time_var            = .false.,
year                = $year,
month               = $month,
day                 = $day,
data_attr_name(1)   = 'long_name',
data_attr_type(1)   = 'char',
data_attr_value(1)  = 'Latitude',
data_attr_name(2)   = 'standard_name',
data_attr_type(2)   = 'char',
data_attr_value(2)  = 'latitude',
data_attr_name(3)   = 'units',
data_attr_type(3)   = 'char',
data_attr_value(3)  = 'degrees_north',
/
EOF
$asc2netcdf asc2netcdf.deck

# copy to other netcdf files
\cp -vaf $air_file_total $rad_file_total
\cp -vaf $air_file_total $prc_file_total
\cp -vaf $air_file_total $apcp_file_total


######## define starting and stopping hours (and step) #################
set hour_start = 0
set hour_stop = 84
set hour_step = 3
set hour_output_int = 12
set rad_start = 3
set prate_start = 3


######## loop over all of the times we wish to retrieve ################
set hour_num = $hour_start
set last_hour = $hour_start
while ($hour_num <= $hour_stop)

  if ($hour_num >= 100) then
    set hour = $hour_num
    set hour_3 = $hour_num
  else if ($hour_num >= 10) then
    set hour = $hour_num
    set hour_3 = 0$hour_num
  else
    set hour = 0$hour_num
    set hour_3 = 00$hour_num
  endif


######## define the file names for this hour ###########################
  set grib_prefix = 'nam.t00z.awip12'
  set grib_suffix = '.tm00'
  set grib_file = $grib_prefix$hour$grib_suffix
  \rm -f $grib_file

  echo ' ' 
  echo $date_string $hour


######## file retrieval block with automated retries ###################
# maximum number of attempts
  set max_attempt = 10

# number of seconds to sleep between attempts
  set sleep_time = 10

# initialize variables that determine whether to iterate or not
  set success = 0
  set attempt = 0

# repeat until success, or maximum attempts
  while ( !($success) && ($attempt < $max_attempt) )

# increment the attempt counter
    @ attempt++

# now grab the file for this time
    date
    echo 'retrieving ' $url_dir$grib_file -- attempt $attempt
    wget -N -nv --tries=10 $url_dir$grib_file

# were we successful?
    set success = (`filetest -e $grib_file`)

# if not successful, then sleep for a bit
    if ( !($success) ) then
      echo failure!  sleeping for $sleep_time
      sleep $sleep_time
    endif

# whitespace
    echo

  end

# if after all that, we failed, then exit the script
  if ( !($success) ) then
    echo file retrieval failure!
    exit
  endif

# if we made it to this point, then we should have been successful
  echo successful file retrieval
######## file retrieval block with automated retries ###################


######## extract data from grib_file (to ascii)  #######################
# extract the 10m u-velocity
  $wgrib -v $grib_file                                    \
      | grep UGRD                                         \
      | grep '10 m above gnd'                             \
      | $wgrib -i -text $grib_file                        \
         -o uwind.txt

# extract the 10m v-velocity
  $wgrib -v $grib_file                                    \
      | grep VGRD                                         \
      | grep '10 m above gnd'                             \
      | $wgrib -i -text $grib_file                        \
         -o vwind.txt

# extract the MSL surface pressure
  $wgrib -v $grib_file                                    \
      | grep PRMSL                                        \
      | $wgrib -i -text $grib_file                        \
         -o prmsl.txt

# extract the 2m temperature
  $wgrib -v $grib_file                                    \
      | grep TMP                                          \
      | grep '2 m above gnd'                              \
      | $wgrib -i -text $grib_file                        \
         -o stmp.txt

# extract the 2m humidity
  $wgrib -v $grib_file                                    \
      | grep SPFH                                         \
      | grep '2 m above gnd'                              \
      | $wgrib -i -text $grib_file                        \
         -o spfh.txt

# extract the downwelling LW flux at the surface
  $wgrib -v $grib_file                                    \
      | grep DLWRF                                        \
      | grep sfc                                          \
      | grep fcst                                         \
      | $wgrib -i -text $grib_file                        \
         -o dlwrf.txt

# extract the downwelling SW flux at the surface
  $wgrib -v $grib_file                                    \
      | grep DSWRF                                        \
      | grep sfc                                          \
      | grep fcst                                         \
      | $wgrib -i -text $grib_file                        \
         -o dswrf.txt

# extract the total precipitation accumulation at the surface
  set time_range = $last_hour"-"$hour_num"hr"
  $wgrib -v $grib_file                                    \
      | grep APCP                                         \
      | grep sfc                                          \
      | grep $time_range                                  \
      | $wgrib -i -text $grib_file                        \
         -o apcp.txt

# remove the grib file once we're done with it
  \rm -f $grib_file


######## add time-varying data to netcdf files #########################

# u
  cat > asc2netcdf.deck << EOF
  &inputs_list
  in_data             = 'uwind.txt',
  out_file            = '$air_file_total',
  data_name           = 'uwind',
  time_var            = .true.,
  year                = $year,
  month               = $month,
  day                 = $day,
  hour                = $hour_num,
  data_attr_name(1)   = 'long_name',
  data_attr_type(1)   = 'char',
  data_attr_value(1)  = 'Surface Eastward Air Velocity (10m AGL)',
  data_attr_name(2)   = 'standard_name',
  data_attr_type(2)   = 'char',
  data_attr_value(2)  = 'eastward_wind',
  data_attr_name(3)   = 'units',
  data_attr_type(3)   = 'char',
  data_attr_value(3)  = 'm/s',
  /
EOF
  $asc2netcdf asc2netcdf.deck

# v
  cat > asc2netcdf.deck << EOF
  &inputs_list
  in_data             = 'vwind.txt',
  out_file            = '$air_file_total',
  data_name           = 'vwind',
  time_var            = .true.,
  year                = $year,
  month               = $month,
  day                 = $day,
  hour                = $hour_num,
  data_attr_name(1)   = 'long_name',
  data_attr_type(1)   = 'char',
  data_attr_value(1)  = 'Surface Northward Air Velocity (10m AGL)',
  data_attr_name(2)   = 'standard_name',
  data_attr_type(2)   = 'char',
  data_attr_value(2)  = 'northward_wind',
  data_attr_name(3)   = 'units',
  data_attr_type(3)   = 'char',
  data_attr_value(3)  = 'm/s',
  /
EOF
  $asc2netcdf asc2netcdf.deck

# prmsl
  cat > asc2netcdf.deck << EOF
  &inputs_list
  in_data             = 'prmsl.txt',
  out_file            = '$air_file_total',
  data_name           = 'prmsl',
  time_var            = .true.,
  year                = $year,
  month               = $month,
  day                 = $day,
  hour                = $hour_num,
  data_attr_name(1)   = 'long_name',
  data_attr_type(1)   = 'char',
  data_attr_value(1)  = 'Pressure reduced to MSL',
  data_attr_name(2)   = 'standard_name',
  data_attr_type(2)   = 'char',
  data_attr_value(2)  = 'air_pressure_at_sea_level',
  data_attr_name(3)   = 'units',
  data_attr_type(3)   = 'char',
  data_attr_value(3)  = 'Pa',
  /
EOF
  $asc2netcdf asc2netcdf.deck

# stmp
  cat > asc2netcdf.deck << EOF
  &inputs_list
  in_data             = 'stmp.txt',
  out_file            = '$air_file_total',
  data_name           = 'stmp',
  time_var            = .true.,
  year                = $year,
  month               = $month,
  day                 = $day,
  hour                = $hour_num,
  data_attr_name(1)   = 'long_name',
  data_attr_type(1)   = 'char',
  data_attr_value(1)  = 'Surface Air Temperature (2m AGL)',
  data_attr_name(2)   = 'standard_name',
  data_attr_type(2)   = 'char',
  data_attr_value(2)  = 'air_temperature',
  data_attr_name(3)   = 'units',
  data_attr_type(3)   = 'char',
  data_attr_value(3)  = 'K',
  /
EOF
  $asc2netcdf asc2netcdf.deck

# spfh
  cat > asc2netcdf.deck << EOF
  &inputs_list
  in_data             = 'spfh.txt',
  out_file            = '$air_file_total',
  data_name           = 'spfh',
  time_var            = .true.,
  year                = $year,
  month               = $month,
  day                 = $day,
  hour                = $hour_num,
  data_attr_name(1)   = 'long_name',
  data_attr_type(1)   = 'char',
  data_attr_value(1)  = 'Surface Specific Humidity (2m AGL)',
  data_attr_name(2)   = 'standard_name',
  data_attr_type(2)   = 'char',
  data_attr_value(2)  = 'specific_humidity',
  data_attr_name(3)   = 'units',
  data_attr_type(3)   = 'char',
  data_attr_value(3)  = '1',
  /
EOF
  $asc2netcdf asc2netcdf.deck

# apcp
  cat > asc2netcdf.deck << EOF
  &inputs_list
  in_data             = 'apcp.txt',
  out_file            = '$apcp_file_total',
  data_name           = 'apcp',
  time_var            = .true.,
  year                = $year,
  month               = $month,
  day                 = $day,
  hour                = $hour_num,
  data_attr_name(1)   = 'long_name',
  data_attr_type(1)   = 'char',
  data_attr_value(1)  = 'Total Precipitation',
  data_attr_name(2)   = 'standard_name',
  data_attr_type(2)   = 'char',
  data_attr_value(2)  = 'precipitation_amount',
  data_attr_name(3)   = 'units',
  data_attr_type(3)   = 'char',
  data_attr_value(3)  = 'kg/m^2',
  /
EOF
  $asc2netcdf asc2netcdf.deck

  if ($hour_num >= $prate_start) then

# calculate the time that the average is valid for
    set time_avg = `echo \($hour_num - 0.5 \* $hour_step\) | bc` >& /dev/null
    echo $hour_num $time_avg

# convert precipitation from accumulated value to average rate
    $nc_make_prate  $apcp_file_total $hour_num $hour_step > prate_avg.txt

# prate
    cat > asc2netcdf.deck << EOF
    &inputs_list
    in_data             = 'prate_avg.txt',
    out_file            = '$prc_file_total',
    data_name           = 'prate',
    time_var            = .true.,
    year                = $year,
    month               = $month,
    day                 = $day,
    hour                = $time_avg,
    data_attr_name(1)   = 'long_name',
    data_attr_type(1)   = 'char',
    data_attr_value(1)  = 'Surface Precipitation Rate',
    data_attr_name(2)   = 'standard_name',
    data_attr_type(2)   = 'char',
    data_attr_value(2)  = 'precipitation_flux',
    data_attr_name(3)   = 'units',
    data_attr_type(3)   = 'char',
    data_attr_value(3)  = 'kg/m^2/s',
    /
EOF
    $asc2netcdf asc2netcdf.deck

  endif


  if ($hour_num >= $rad_start) then

# dlwrf
    cat > asc2netcdf.deck << EOF
    &inputs_list
    in_data             = 'dlwrf.txt',
    out_file            = '$rad_file_total',
    data_name           = 'dlwrf',
    time_var            = .true.,
    year                = $year,
    month               = $month,
    day                 = $day,
    hour                = $hour_num,
    data_attr_name(1)   = 'long_name',
    data_attr_type(1)   = 'char',
    data_attr_value(1)  = 'Downward Long Wave Radiation Flux',
    data_attr_name(2)   = 'standard_name',
    data_attr_type(2)   = 'char',
    data_attr_value(2)  = 'surface_downwelling_longwave_flux_in_air',
    data_attr_name(3)   = 'units',
    data_attr_type(3)   = 'char',
    data_attr_value(3)  = 'W/m^2',
    /
EOF
    $asc2netcdf asc2netcdf.deck

# dswrf
    cat > asc2netcdf.deck << EOF
    &inputs_list
    in_data             = 'dswrf.txt',
    out_file            = '$rad_file_total',
    data_name           = 'dswrf',
    time_var            = .true.,
    year                = $year,
    month               = $month,
    day                 = $day,
    hour                = $hour_num,
    data_attr_name(1)   = 'long_name',
    data_attr_type(1)   = 'char',
    data_attr_value(1)  = 'Downward Short Wave Radiation Flux',
    data_attr_name(2)   = 'standard_name',
    data_attr_type(2)   = 'char',
    data_attr_value(2)  = 'surface_downwelling_shortwave_flux_in_air',
    data_attr_name(3)   = 'units',
    data_attr_type(3)   = 'char',
    data_attr_value(3)  = 'W/m^2',
    /
EOF
    $asc2netcdf asc2netcdf.deck

  endif


######## check to see whether we're at an interval to save the output ##
######## or at the end of time loop - or whatever ######################

  set interval_modulus = `echo \($hour_num \% $hour_output_int\) | bc` >& /dev/null

  if ($interval_modulus == 0) then
    if ($hour_num == 0) then
      set complete_processing = false
    else
      set complete_processing = true
    endif
  else if ($hour_num == $hour_stop) then
    set complete_processing = true
  else
    set complete_processing = false
  endif

  
######## complete processing of data if at one of those times ##########
  if ($complete_processing == true) then

    echo
    echo completing processing and saving output at hour_num = $hour_num
    echo


######## copy air data to tmp file (so we don't modify original) ###############
    set tmp_air_file = tmp_air_file.nc
    \cp -f $air_file_total $tmp_air_file


######## remove any preexisting versions of the local files ############
    \rm -f $air_file_local $rad_file_local $prc_file_local


######## rotate winds to geographic coordinates (use tmp_file) #########
    cat > nc_rotate_wind.deck << EOF
    &inputs_list
    in_file             = '$tmp_air_file',
    proj                = $proj,
    scale               = $scale,
    true_lat(1)         = $true_lat_1,
    true_lat(2)         = $true_lat_2,
    true_lon            = $true_lon,
    /
EOF
    $nc_rotate_wind nc_rotate_wind.deck


######## reduce size of domain #########################################

# air
    cat > nc_restrict.deck << EOF
    &inputs_list
    in_file             = '$tmp_air_file',
    out_file            = '$air_file_local',
    lon_min             = $lon_min,
    lon_max             = $lon_max,
    lat_min             = $lat_min,
    lat_max             = $lat_max,
    data_name(1)        = 'uwind',
    data_name(2)        = 'vwind',
    data_name(3)        = 'prmsl',
    data_name(4)        = 'stmp',
    data_name(5)        = 'spfh',
    /
EOF
    $nc_restrict nc_restrict.deck

# rad
    cat > nc_restrict.deck << EOF
    &inputs_list
    in_file             = '$rad_file_total',
    out_file            = '$rad_file_local',
    lon_min             = $lon_min,
    lon_max             = $lon_max,
    lat_min             = $lat_min,
    lat_max             = $lat_max,
    data_name(1)        = 'dlwrf',
    data_name(2)        = 'dswrf',
    /
EOF
    $nc_restrict nc_restrict.deck

# prc
    cat > nc_restrict.deck << EOF
    &inputs_list
    in_file             = '$prc_file_total',
    out_file            = '$prc_file_local',
    lon_min             = $lon_min,
    lon_max             = $lon_max,
    lat_min             = $lat_min,
    lat_max             = $lat_max,
    data_name(1)        = 'prate',
    /
EOF
    $nc_restrict nc_restrict.deck


######## remove the temp air file #########################################
    \rm -f $tmp_air_file


######## copy the local (reduced size) netcdf files to archive #########
    \cp -f $air_file_local $dst_dir
    \cp -f $rad_file_local $dst_dir
    \cp -f $prc_file_local $dst_dir


  endif

# bottom of loop of hours for file retrieval (and data extraction)
  set last_hour = $hour_num
  @ next_hour = $hour_num + $hour_step
  set hour_num = $next_hour
end
######## bottom of time loop for file retrieval and data extraction ####


######## files were already moved to archive, clean up now #############
\rm -f $work_dir*.txt $work_dir*.nc $work_dir*.deck

