#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#
chdir $thisdir;
#convert to utm
print "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.utm.26918 1 0 \n";
system "./proj_wrap.pl epsg:26918 2 1 hgrid.ll hgrid.utm.26918 1 0 ";

print "\n";
print "Changed to $thisdir \n";
print "hgrid.utm.26918 updated\n";
print "\n";

print "./cpp < cpp.in\n";
system "./cpp < cpp.in";
print "mv out_hgrid.ll hgrid.cpp\n";
system "mv out_hgrid.ll hgrid.cpp";

system("cat bnd >> hgrid.cpp");
system("cat bnd >> hgrid.utm.26918");

system("rm ../hgrid.cpp");
system("rm ../hgrid.utm.26918");

system("cp hgrid.cpp hgrid.utm.26918 ../");


print "\n";
print "\n";

