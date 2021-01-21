#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print("You are responsible for providing 'bottom_fric.in', see a sample in this dir\n");
print("You may need to recompile 'gen_bottom_fric.f90'\n");
print("e.g: ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_bottom_fric gen_bottom_fric.f90\n\n\n");


#dirs
#$rundir = cwd();
$script_dir="../Grid_manipulation/";

#hgrid.gr3 is in lon/lat
system("ln -sf ../hgrid.gr3 .");
system("ln -sf ../vgrid.in .");

#set Delaware Bay =1 in include.gr3
#system("$script_dir/auto_edit_region 0 include.reg hgrid.gr3 1 0");
#move("out.gr3","include.gr3");
#
unless (-e "bottom_fric.in") {
  print ("'bottom_fric.in' not found in the current dir, see examples in ./Sample_bottom_fric/\n");
  exit;
}
system("./gen_bottom_fric < bottom_fric.in");


#set special values in GoME regions
system("$script_dir/auto_edit_region 1 GoME_1.reg bottom_friction.gr3 0.05");
move("out.gr3","bottom_friction.gr3");
system("$script_dir/auto_edit_region 1 GoME_2.reg bottom_friction.gr3 0.05");
move("out.gr3","bottom_friction.gr3");


unlink("../manning.gr3");
copy("bottom_friction.gr3","../manning.gr3");


