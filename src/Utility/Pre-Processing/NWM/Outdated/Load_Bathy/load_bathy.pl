#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#dirs-change these to your own dirs
$script_dir="./";
$thisdir = cwd();
$etopo_dir="../DEM/";
$usgs_dir="../DEM_USGS/";

system("ln -sf ../auto_edit_region .");




print "make sure you have hgrid.utm.26918 and bnd (boundary info) under your run dir (../)\n";
print "make sure you have USGS DEM in $usgs_dir, download link: https://drive.google.com/open?id=1plYccAcY0tk-Ze1hzc_jxfqi6zOkYNGB\n";
print "make sure the module 'proj' (PROJ.4 - Cartographic Projections Library) is loaded\n";
{
  local( $| ) = ( 1 );
  print "Press <Enter> or <Return> to continue: \n";
  my $resp = <STDIN>;
}

system("ln -sf ../auto_edit_region .");
system("cp ../hgrid.utm.26918 .");

print "./proj_wrap.pl epsg:26918 1 1 hgrid.utm.26918 hgrid.ll 1 0\n";
system "./proj_wrap.pl epsg:26918 1 1 hgrid.utm.26918 hgrid.ll 1 0";

print "cp hgrid.ll $etopo_dir/hgrid.old\n";
system "cp hgrid.ll $etopo_dir/hgrid.old";

print "chdir $etopo_dir\n";
chdir $etopo_dir;

print "./load_asc.pl\n";
system "./load_asc.pl";


print "cp hgrid.old hgrid.ll\n";
system "cp hgrid.old hgrid.ll";

print "readlink -f hgrid.ll\n";
system "readlink -f hgrid.ll";
print "head hgrid.ll\n";
system "head hgrid.ll";

print "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.gr3 1 0\n";
system "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.gr3 1 0";

print "cp hgrid.gr3 $usgs_dir/hgrid.old\n";
system "cp hgrid.gr3 $usgs_dir/hgrid.old";

print "chdir $usgs_dir\n";
chdir $usgs_dir;

print "./load_asc.pl\n";
system "./load_asc.pl";

print "cp hgrid.old hgrid.gr3\n";
system "cp hgrid.old hgrid.gr3";

print "readlink -f hgrid.gr3\n";
system "readlink -f hgrid.gr3";
print "head hgrid.gr3\n";
system "head hgrid.gr3";

print "./proj_wrap.pl epsg:26918 1 1 hgrid.gr3 hgrid.ll 1 0\n";
system "./proj_wrap.pl epsg:26918 1 1 hgrid.gr3 hgrid.ll 1 0";
print "head hgrid.ll\n";
system "head hgrid.ll";

#copy bathymetry loaded grid to thisdir
print "cp hgrid.ll $thisdir/\n";
system "cp hgrid.ll $thisdir/";
#
chdir $thisdir;
#convert to utm
print "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.utm.26918 1 0 \n";
system "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.utm.26918 1 0 ";

print "\n";
print "Changed to $thisdir \n";
print "hgrid.ll and hgrid.utm.26918 updated\n";
print "\n";

#manual:
#{
#  local( $| ) = ( 1 );
#  print "Manually load CB bathymetry using (CB_utm.reg) and set minimum depth to 5 m; then set minimum depths for DB-CB canal (5 m; skipped), ocean (5 m), SE corner (5 m; skipped), Riegelsville (-41 m), then\n";
#  my $resp = <STDIN>;
#}
system("$script_dir/auto_edit_region 1 CB_utm.reg hgrid.utm.26918 5");
move("out.gr3","hgrid.utm.26918");
system("$script_dir/auto_edit_region 1 DBCB_utm.reg hgrid.utm.26918 5");
move("out.gr3","hgrid.utm.26918");
system("$script_dir/auto_edit_region 1 Ocean_utm.reg hgrid.utm.26918 5");
move("out.gr3","hgrid.utm.26918");
system("$script_dir/auto_edit_region 1 Rieg.reg hgrid.utm.26918 -41");
move("out.gr3","hgrid.utm.26918");

print "\n";
print "\n";

print "./proj_wrap.pl epsg:26918 1 1 hgrid.utm.26918 hgrid.ll 1 0\n";
system "./proj_wrap.pl epsg:26918 1 1 hgrid.utm.26918 hgrid.ll 1 0";
print "ln -sf hgrid.ll hgrid.gr3\n";
system "ln -sf hgrid.ll hgrid.gr3";
print "./cpp < cpp.in\n";
system "./cpp < cpp.in";
print "mv out_hgrid.ll hgrid.cpp\n";
system "mv out_hgrid.ll hgrid.cpp";

system("cat bnd >> hgrid.gr3");
system("cat bnd >> hgrid.cpp");
system("cat bnd >> hgrid.utm.26918");
system("rm ../hgrid.*");
system("cp hgrid.* ../");


print "\n";
print "\n";

