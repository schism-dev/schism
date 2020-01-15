#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print("You may need to recompile 'gen_slope_filter.f90'\n");
print("e.g: ifort -O2 -CB -o gen_slope_filter gen_slope_filter.f90\n\n\n");



#dirs
#$rundir = cwd();
$script_dir="../Grid_manipulation/";

#UTM grid
system("ln -sf ../hgrid.utm.26918 hgrid.gr3");

#set shapiro outside Delaware Bay and Chesapeake Bay
system("./gen_slope_filter < gen_slope_filter.in");
move("slope_filter.gr3","shapiro.gr3");

##set shapiro inside Delaware Bay and Chesapeake Bay
#system("$script_dir/auto_edit_region 1 0.1.reg shapiro.gr3 0.1");
#move("out.gr3","shapiro.gr3");
#system("$script_dir/auto_edit_region 1 0.2.reg shapiro.gr3 0.2");
#move("out.gr3","shapiro.gr3");
#
##set max shapiro = 0.5 inside Chesapeake Bay
#system("$script_dir/auto_edit_region 1 CB_0.5.reg shapiro.gr3 0.5");
#move("out.gr3","shapiro.gr3");

unlink("../shapiro.gr3");
copy("shapiro.gr3","../shapiro.gr3");

print("Done.\n");
