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

system("$script_dir/auto_edit_prop fluxflag/0-.reg hgrid.ll 0 -1");
move("out.prop","default.prop");
system("$script_dir/auto_edit_prop fluxflag/0+.reg hgrid.ll 1 -9999");
move("out.prop","default.prop");

for (my $i=1; $i <= 2; $i++) {
    $region="$i-.reg";
    system("$script_dir/auto_edit_prop fluxflag/$region hgrid.ll $i -9999");
    move("out.prop","default.prop");
    $region="$i+.reg";
    $i1=$i+1;
    system("$script_dir/auto_edit_prop fluxflag/$region hgrid.ll $i1 -9999");
    move("out.prop","default.prop");
}

move("default.prop","fluxflag.prop");

copy("fluxflag.prop","../fluxflag.prop");

print("Done.\n")


