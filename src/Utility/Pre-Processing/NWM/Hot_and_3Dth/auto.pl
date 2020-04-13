#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print(">>>>>>>>>>>>Make sure you have netcdf libraries>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>Make you set the hycom data dir correctly:>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>Make sure the following excutables work on your machine:>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>    ./gen_hot_3Dth_from_hycom>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>    ./modify_hot>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>    ./Elev_IC/gen_elev>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>Recompile them if necessary;>>>>>>>>>>>>>\n");
print(">>>>>>>>>>>>see an example compiling cmd in each source code.>>>>>>>>>>>>>\n");
print("\n");

#{
#  local( $| ) = ( 1 );
#  print "Press <Enter> or <Return> to continue: \n";
#  my $resp = <STDIN>;
#}

#dirs
$script_dir="../Grid_manipulation/";
$hycom_dir="../HYCOM_FLORENCE_PERIOD/";

$thisdir=cwd();
chdir("..");
chdir($thisdir);

#UTM grid
system("ln -sf ../hgrid.* .");
system("ln -sf ../vgrid.in .");
system("ln -sf $hycom_dir/*.nc .");


#--------------------hycom---------------------
print("-------------------------------\n");
print(">>>>>>>>>>>>HYCOM>>>>>>>>>>>>>\n");
#set CB=1 and DB=1 in estuary.gr3
system("$script_dir/auto_edit_region 0 CB.reg hgrid.gr3 1 0");
move("out.gr3","estuary.gr3");
system("$script_dir/auto_edit_region 0 DB.reg estuary.gr3 1");
move("out.gr3","estuary.gr3");

system("./gen_hot_3Dth_from_hycom.exe");



#--------------------modify---------------------
print("-------------------------------\n");
print(">>>>>>>>>>>>Additional data>>>>>>>>>>>>>\n");

chdir("Elev_IC");
system("./auto.pl");


chdir($thisdir);

#set elev.ic
system("ln -sf Elev_IC/elev.ic DB_elev_ic.gr3");

#set salt.ic
print(">>>>>>>>>>>>Preparing i.c. for S >>>>>>>>>>>>>\n");
system("ln -sf DelawareBay_Data/DB_surf_S_ic.subset.gr3 bg.gr3");
system("rm fg.gr3"); 
#make fg.gr3, to be interpolated on
system("./gen_gr3.pl");
#
#set initial sal to 0 in coastal zones before interpolation
system("$script_dir/auto_edit_region 0 ic_S_0.reg fg.gr3 0");
move("out.gr3","fg.gr3");
#set initial sal to 0 in Missi. R.
system("$script_dir/auto_edit_region 0 ic_S_0_Missi.reg fg.gr3 0");
move("out.gr3","fg.gr3");
print(">>>>>>>>>>>>Done Setting 0 salinity i.c. in coastal zones>>>>>>>>>>>>>\n");

#make include.gr3 for salt interp inside DB
system("$script_dir/auto_edit_region 0 DB.reg hgrid.gr3 1 0");
move("out.gr3","include.gr3");
#do interp between 2 *.gr3
system("$script_dir/interpolate_unstructured");
system("mv fg.new DB_surf_S_ic.gr3");
print(">>>>>>>>>>>>Done interpolating salinity i.c. in Delaware Bay>>>>>>>>>>>>>\n");

#set CB=1 in estuary.gr3
system("$script_dir/auto_edit_region 0 CB.reg hgrid.gr3 1 0");
move("out.gr3","estuary.gr3");
#set "Other coastal zones"=2 in estuary.gr3
system("$script_dir/auto_edit_region 0 DB.reg estuary.gr3 2");
move("out.gr3","estuary.gr3");
system("$script_dir/auto_edit_region 0 ic_S_0.reg estuary.gr3 2");
move("out.gr3","estuary.gr3");
system("$script_dir/auto_edit_region 0 ic_S_0_Missi.reg estuary.gr3 2");
move("out.gr3","estuary.gr3");

system("./modify_hot");


system("rm ../hotstart.nc");
system("mv  hotstart.nc ../hotstart.nc");
system("rm ../*.th.nc");
system("mv  *.th.nc ../");
print(">>>>>>>>>>>>Done.>>>>>>>>>>>>>\n")
