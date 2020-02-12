#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

#dirs
#$rundir = cwd();
$script_dir="../Grid_manipulation/";

#UTM grid
system("ln -sf ../hgrid.* .");

#set tvd.prop
system("$script_dir/auto_edit_prop tvd.rgn hgrid.ll 0 1");
move("out.prop","tvd.prop");
unlink("../tvd.prop");
copy("tvd.prop","../tvd.prop");

#set fluxflag.prop
system("$script_dir/auto_edit_prop fluxflag/0.reg hgrid.ll 0 -1");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/1-.reg hgrid.ll 1 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/0.reg hgrid.ll 0 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/1+.reg hgrid.ll 1 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/2-.reg hgrid.ll 2 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/2+.reg hgrid.ll 2 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/3-.reg hgrid.ll 3 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/3+.reg hgrid.ll 3 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/4-.reg hgrid.ll 4 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/4+.reg hgrid.ll 4 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/5-.reg hgrid.ll 5 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/5+.reg hgrid.ll 5 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/6-.reg hgrid.ll 6 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/6+.reg hgrid.ll 6 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/7-.reg hgrid.ll 7 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/7+.reg hgrid.ll 7 -9999");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/8.reg hgrid.ll 8 -9999");
move("out.prop","fluxflag.prop");

unlink("../fluxflag.prop");
copy("fluxflag.prop","../fluxflag.prop");

print("Done.\n")
