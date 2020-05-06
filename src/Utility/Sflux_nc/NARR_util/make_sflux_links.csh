#! /bin/csh

#### make_sflux_links.csh 3 <start year> <mm> <dd> <end year> <mm> <dd>
#### On a new system, change util_bin and data_grandparent below
######## define locations and names for utilities ######################

set util_bin = ~yinglong/git/schism/src/Utility/Sflux_nc/NARR_util/src/
set cjd_util       = $util_bin'cjd/cjd'


######## set command line options ######################################

# presently set up so that the minimum number of arguments is 7
# (see below)

if ( $7 == "" ) then
  echo
  echo This utility should be executed from within the run directory.
  echo
  echo It creates the sflux subdirectory, as well as the links to data
  echo and the sflux_inputs.txt file.
  echo
  echo you must supply the following command line arguments:
  echo
  echo \(1\) data_source_code
  echo \(2\) year_start
  echo \(3\) month_start
  echo \(4\) day_start
  echo \(5\) year_stop
  echo \(6\) month_stop
  echo \(7\) day_stop
  echo \(8\) utc_start  \(optional argument\)
  echo \(9\) start_hour \(optional argument\)
  echo
  echo Example:
  echo   % make_sflux_links.csh 1 2007 1 1 2007 1 7
  echo
  echo Numeric values for data_source_code:
  echo 1 - NAM
  echo 2 - GFS
  echo 3 - NARR
  echo 4 - NNRP
  exit
else

set data_source_code = $1
set year_start       = $2
set month_start      = $3
set day_start        = $4
set year_stop        = $5
set month_stop       = $6
set day_stop         = $7


# set optional arguments

if ( $8 == "" ) then
  set utc_start  = '8.0'
else
  set utc_start  = $8
endif

if ( $9 == "" ) then
  set start_hour = '0.0'
else
  set start_hour  = $9
endif

# output options to screen
echo
echo data_source_code = $data_source_code
echo
echo year_start       = $year_start
echo month_start      = $month_start
echo day_start        = $day_start
echo
echo year_stop        = $year_stop
echo month_stop       = $month_stop
echo day_stop         = $day_stop
echo
echo utc_start        = $utc_start
echo start_hour       = $start_hour


######## calculate starting and stopping Julian dates ##################

echo
echo Note: including additional day at end \(PST/UTC difference\)

set cjd_start = (`$cjd_util $year_start $month_start $day_start`)
set cjd_stop  = (`$cjd_util $year_stop  $month_stop  $day_stop `)

set cjd_stop = `echo \($cjd_stop + 1\) | bc` >& /dev/null

echo
echo cjd_start = $cjd_start
echo cjd_stop  = $cjd_stop


######## specify the models to be used as data sources, locations, etc #

set src_1 = ""
set src_2 = ""
set domain_1 = "."
set domain_2 = "."
set data_grandparent = '~yinglong/vims20/'

if ($data_source_code == "1") then
  set src_1 = "nam"
  set domain_1 = ".local."

else if ($data_source_code == "2") then
  set src_1 = "gfs"

  set cjd_switch    = (`$cjd_util 2005 10 20`)
  set cjd_switch_p1 = (`$cjd_util 2005 10 21`)

  if ($cjd_start < $cjd_switch) then
    if ($cjd_stop <= $cjd_switch) then
      set domain_1 = ".local."
    else if ($cjd_stop == $cjd_switch_p1) then
      set domain_1 = ".local."
      echo
      echo Removing additional day, due to GFS local/global switch!
      set cjd_stop = $cjd_switch
    else
      echo
      echo Run extends across GFS local/global switch!
      exit
    endif
  endif

else if ($data_source_code == "3") then
  set src_1 = "narr"

else if ($data_source_code == "4") then
  set src_1 = "nnrp"

else
  echo
  echo undefined data_source_code. . .
  exit
endif

echo
echo Model Sources
echo src_1 = $src_1
echo src_2 = $src_2


######## create the destination directory (replace if it exists) #######

# create name of destination directory
set local_dir = `pwd`'/'
set dest_dir  = $local_dir'sflux/'
set old_dest_dir  = $local_dir'sflux.old/'

# check to see if destination directory exists
set dest_dir_exist = (`filetest -e $dest_dir`)

# fail if directory already exists
if ( ($dest_dir_exist) ) then
  echo
  echo destination directory already exists!
  echo renaming it!
  \rm -rf $old_dest_dir
  \mv $dest_dir $old_dest_dir
endif

# create destination directory
mkdir $dest_dir

echo
echo created destination directory:
echo $dest_dir


######## create the sflux_inputs file ##################################

set inputs_file = 'sflux_inputs.txt'

echo
echo Creating sflux_inputs file at:
echo "  " $dest_dir$inputs_file

cat > $dest_dir$inputs_file << EOF
&sflux_inputs
/
EOF


######## loop over file types and dates, creating links ################

# loop over possible file types
foreach file_type (air rad prc)
  echo
  echo
  echo -----------------------------------------------------------------
  echo Creating links for file_type: $file_type

# initialize counter for files of this type
  set type_counter = 0

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

# determine where data will come from for this day (and go there)
    set data_parent = $data_grandparent$src_1'/'
    set data_dir    = $data_parent$year'_'$month'/'
    cd $data_dir

# loop over files that satisfy required name (this skips missing files)
    foreach in_file (`ls $src_1'_'$file_type$domain_1$date_string'.nc'`)
      set type_counter = `echo \($type_counter + 1\) | bc` >& /dev/null

      if ($type_counter >= 1000) then
        set type_counter_label = $type_counter
      else if ($type_counter >= 100) then
        set type_counter_label = 0$type_counter
      else if ($type_counter >= 10) then
        set type_counter_label = 00$type_counter
      else
        set type_counter_label = 000$type_counter
      endif

      set link_name = sflux_$file_type"_1."$type_counter_label.nc

      echo
      echo creating link from:
      echo "  " $data_dir$in_file
      echo to:
      echo "  " $dest_dir$link_name

      \ln -sf $data_dir$in_file $dest_dir$link_name

# end of foreach loop that selects files
    end

# bottom of loop over Julian days
    set cjd = `echo \($cjd + 1\) | bc` >& /dev/null
  end

# end of loop over file types
end

######## exit the shell ################################################

echo
exit

