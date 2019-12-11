#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#dirs
$thisdir=cwd();
chdir("..");
$rundir = cwd();
chdir($thisdir);
$script_dir="./";

#UTM grid
system("ln -sf ../hgrid.gr3 .");
system("ln -sf ../hgrid.ll .");
system("ln -sf ../vgrid.in .");

chdir("Nudge_gr3");
system("ln -sf $rundir/hgrid.ll hgrid.gr3");
system("./gen_nudge2 < gen_nudge2.in");
unlink("$rundir/SAL_nudge.gr3");
copy("nudge.gr3","$rundir/SAL_nudge.gr3");
unlink("$rundir/TEM_nudge.gr3");
copy("nudge.gr3","$rundir/TEM_nudge.gr3");


chdir($thisdir);

#set include zone, which is a bit larger than the actual nudging zone
system("$script_dir/auto_edit_region 0 include.reg hgrid.gr3 1 0");
move("out.gr3","include.gr3");

#generate nudging files
print("Make sure you have netcdf libraries\n");
{
  local( $| ) = ( 1 );
  print "Press <Enter> or <Return> to continue: \n";
  my $resp = <STDIN>;
}
system("./gen_nudge_from_hycom");

unlink("../SAL_nu.nc");
unlink("../TEM_nu.nc");

copy("SAL_nu.nc","../SAL_nu.nc");
copy("TEM_nu.nc","../TEM_nu.nc");

print("Done.\n")
