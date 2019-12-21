#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print("You may need to recompile 'gen_hot_3Dth_from_hycom.f90',\n");
print("see an example compiling cmd in the source code.\n");
print("Make sure you have netcdf libraries\n");
#{
#  local( $| ) = ( 1 );
#  print "Press <Enter> or <Return> to continue: \n";
#  my $resp = <STDIN>;
#}

#dirs
$hycom_dir="../HYCOM_IRENE_PERIOD/";

#$thisdir=cwd();

#UTM grid
system("ln -sf ../hgrid.* .");
system("ln -sf ../vgrid.in .");
system("ln -sf $hycom_dir/*.nc .");


system("./gen_gr3.pl");

#generate hotstart.nc and *D.th.nc
system("./gen_hot_3Dth_from_hycom");

unlink glob "../*D.th.nc";
unlink("./hotstart.nc");

system("cp *D.th.nc ../");

print("Done.\n")
