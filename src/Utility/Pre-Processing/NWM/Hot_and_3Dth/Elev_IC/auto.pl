#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;


print("You may need to recompile 'gen_elev.f90'\n");
print("e.g: ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_elev gen_elev.f90\n\n\n");

#dirs
#$rundir = cwd();
$script_dir="../../Grid_manipulation/";

#UTM grid
system("ln -sf ../hgrid.utm.26918 hgrid.gr3");
system("ln -sf ../vgrid.in .");

#set z=1 in Delaware Bay
system("$script_dir/auto_edit_region 0 include.reg hgrid.gr3 1 0");
move("out.gr3","include.gr3");
system("./gen_elev");
unlink("source_sink.in","vsource.th","msource.th","rough.gr3");

#unlink("../elev.ic");
#copy("elev.ic","../elev.ic");

print("Done.\n")
