#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#dirs

print "!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
print "make sure the proj/ module is loaded (module load proj)\n";
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
print "\n";
print "\n";

system "ln -sf ../hgrid.ll .";

print "converting lon/lat to UTM\n";
system "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.utm.26918 1 0";

print "converting lon/lat to cpp\n";
system "./cpp < cpp.in";
print "mv out_hgrid.ll hgrid.cpp\n";
system "mv out_hgrid.ll hgrid.cpp";

#put bnd from hgrid.gr3 into hgrid.cpp and hgrid.utm.*
system "cp hgrid.ll bnd";
#system "sed -i '0,/boundary/d' bnd";
system "sed -i '/bound/,\$\!d' bnd";
system "cat bnd >> hgrid.cpp";
system "cat bnd >> hgrid.utm.26918";

#put newly made hgrid.* under run dir and print the first few lines of each
system "rm ../hgrid.utm.26918 ../hgrid.cpp";
system "cp hgrid.utm.26918 hgrid.cpp ../";

print "First few lines of hgrid.*\n";
system "head ../hgrid.utm.26918 ../hgrid.cpp";

print "\n";
