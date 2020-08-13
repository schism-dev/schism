#!/usr/bin/perl -w

# Load all .asc 
# Run in this dir

$code='./interpolate_depth_structured2';

system "cp -L hgrid.old hgrid.old.0"; #save a copy

open(IN,">in");
print IN "-1\n -0.000\n";
close(IN);
print "loading $a\n";
system("ln -sf etopo1.asc struc.grd");
system("$code < in ; mv hgrid.new hgrid.old;");
