#! /usr/bin/perl -w

# Simple script to link old sflux names to new format (May 2025). Move sflux/ sflux.0 first and then do this inside sflux/
# Change # of files below

for($i=1;$i<=32; $i++)
{
  $j=sprintf('%04d',$i);
  $file="sflux_air_1."."$j".".nc";
  $file2="sflux_air_1."."$i".".nc";
  print "doing $file\n";
  system "ln -sf ../sflux.0/$file $file2";

  $file="sflux_rad_1."."$j".".nc";
  $file2="sflux_rad_1."."$i".".nc";
  print "doing $file\n";
  system "ln -sf ../sflux.0/$file $file2";

  $file="sflux_prc_1."."$j".".nc";
  $file2="sflux_prc_1."."$i".".nc";
  print "doing $file\n";
  system "ln -sf ../sflux.0/$file $file2";

  #.txt
  system "ln -sf ../sflux.0/sflux*.txt";
}
