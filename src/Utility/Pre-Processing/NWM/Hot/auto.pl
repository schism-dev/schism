#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print("You may need to recompile 'gen_source2.f90'\n");
print("e.g: ifort -O2 -mcmodel=medium -assume byterecl -o gen_hot_3Dth_from_hycom gen_hot_3Dth_from_hycom.f90 ../UtilLib/compute_zcor.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff\n\n\n");
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
#set CB=1 and DB=2 in estuary.gr3
system("$script_dir/auto_edit_region 0 CB.reg hgrid.gr3 1 0");
move("out.gr3","estuary.gr3");
system("$script_dir/auto_edit_region 0 DB.reg estuary.gr3 2");
move("out.gr3","estuary.gr3");

#generate hotstart.nc and *D.th.nc
system("./gen_hot_3Dth_from_hycom");

unlink("../hotstart.nc");
unlink("../*D.th.nc");

copy("hotstart.nc","../hotstart.nc");
copy("*D.th.nc","../");

print("Done.\n")
