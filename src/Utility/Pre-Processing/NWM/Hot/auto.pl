#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print("You may need to recompile 'gen_source2.f90',\n");
print("see an example compiling cmd in the source code.\n");
print("Make sure you have netcdf libraries\n");
#{
#  local( $| ) = ( 1 );
#  print "Press <Enter> or <Return> to continue: \n";
#  my $resp = <STDIN>;
#}

#dirs
$script_dir="../Grid_manipulation/";
$hycom_dir="../HYCOM_IRENE_PERIOD/";

$thisdir=cwd();
chdir("..");
chdir($thisdir);

#UTM grid
system("ln -sf ../hgrid.* .");
system("ln -sf ../vgrid.in .");
system("ln -sf $hycom_dir/*.nc .");

chdir("Elev_IC");
system("./auto.pl");


chdir($thisdir);

#set elev.ic
system("ln -sf Elev_IC/elev.ic DB_elev_ic.gr3");

#set salt.ic
print("Interpolating salinity i.c. in Delaware Bay\n");
system("ln -sf DB_surf_S_ic.subset.gr3 bg.gr3");
system("rm fg.gr3"); 
#make fg.gr3, to be interpolated on
system("./gen_gr3.pl");
#set initial sal to 0 in DB before interpolation
system("$script_dir/auto_edit_region 0 DB.reg fg.gr3 0");
move("out.gr3","fg.gr3");
#make include.gr3 for salt interp
system("$script_dir/auto_edit_region 0 DB.reg hgrid.gr3 1 0");
move("out.gr3","include.gr3");
#do interp between 2 *.gr3
system("$script_dir/interpolate_unstructured");
system("mv fg.new DB_surf_S_ic.gr3");
print("Done interpolating salinity i.c. in Delaware Bay\n");

#set CB=1 and DB=2 in estuary.gr3
system("$script_dir/auto_edit_region 0 CB.reg hgrid.gr3 1 0");
move("out.gr3","estuary.gr3");
system("$script_dir/auto_edit_region 0 DB.reg estuary.gr3 2");
move("out.gr3","estuary.gr3");

#generate hotstart.nc and *D.th.nc
system("./gen_hot_3Dth_from_hycom");

unlink("../hotstart.nc");
unlink("./*D.th.nc");

copy("hotstart.nc","../hotstart.nc");

print("Done.\n")
