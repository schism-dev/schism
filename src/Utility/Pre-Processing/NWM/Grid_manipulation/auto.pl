#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#dirs

print "make sure the proj/ module is loaded (module load proj)\n\n";
system "ln -sf ../hgrid.ll .";

print "converting lon/lat to UTM\n";
system "~/schism_trunk/src/Utility/Pre-Processing/proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.utm.26918 1 0";

print "converting lon/lat to cpp\n";
system "./cpp < cpp.in";
print "mv out_hgrid.ll hgrid.cpp\n";
system "mv out_hgrid.ll hgrid.cpp";

system "rm ../hgrid.utm.26918 ../hgrid.cpp";
system "cp hgrid.utm.26918 hgrid.cpp ../";

print "\n";
