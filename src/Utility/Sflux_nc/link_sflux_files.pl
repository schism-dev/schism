#! /usr/bin/perl -w

#Link sflux_air* etc to nc files generated from other programs for SCHISM runs. 
#Only works within a year; only a prototype.
@month_day=(31,28,31,30,31,30,31,31,30,31,30,31);
$year=2011;
$start_mon=7;
$start_day=1;
$days=2; #no. of SCHISM run days +1

if($year % 4 ==0) {$month_day[1]=29;}

$day=$start_day-1; $mon=$start_mon;
for($i=1;$i<=$days;$i++)
{
  $day++;
  if($day>$month_day[$mon-1])
  {
    $day=$day-$month_day[$mon-1];
    $mon++;
    if($mon>12) {die "Exceeded 1 year\n";}
  }
  $day_char1=sprintf("%2.2d",$day);
  $day_char2=sprintf("%3.3d",$i);
  $mon_char1=sprintf("%2.2d",$mon);
  print "Doing $year/$mon/$day\n";
  #system "ln -sf dwd2011$mon_char1$day_char1\_4selfe.nc sflux_air_1.$day_char2\.nc";
  system "ln -sf air_$year$mon_char1$day_char1\.nc sflux_air_1.$day_char2\.nc";
  system "ln -sf rad_$year$mon_char1$day_char1\.nc sflux_rad_1.$day_char2\.nc";
} #for

