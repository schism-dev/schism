#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#dirs
$script_dir="../../Grid_manipulation/";

#hgrid.gr3 is in lon/lat
system("ln -sf ../../hgrid.ll .");

#set inside region = 1 
system("$script_dir/auto_edit_prop sflux2source.reg hgrid.ll 1 0");
move("out.prop","sflux2source.prop");

unlink("../sflux2source.prop");
copy("sflux2source.prop","../sflux2source.prop");


