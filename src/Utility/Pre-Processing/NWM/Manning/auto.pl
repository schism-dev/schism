#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print("You may need to recompile 'gen_source2.f90'\n");
print("e.g: ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_source2 gen_source2.f90\n\n\n");


#dirs
#$rundir = cwd();
$script_dir="../Grid_manipulation/";

system("ln -sf ../hgrid.gr3 .");
system("ln -sf ../vgrid.in .");

#set Delaware Bay =1 in include.gr3
#system("$script_dir/auto_edit_region 0 include.reg hgrid.gr3 1 0");
#move("out.gr3","include.gr3");
#
system("./gen_source2 < gen.in");

unlink("../manning.gr3");
copy("bottom_friction.gr3","../manning.gr3");
