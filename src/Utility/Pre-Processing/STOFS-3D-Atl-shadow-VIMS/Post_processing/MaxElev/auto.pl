#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

# ----------------inputs-------------------------
$start_stack = 53;
$end_stack = 60;
#
# $start_rec = 1;
# $end_rec = 24;

$rundir='../../../../Runs/RUN13k/';
# ----------------end inputs---------------------

# $file_suffix="stack$start_stack.$start_rec-$end_stack.$end_rec";
$file_suffix="stack$start_stack-$end_stack";

system("cp max_elev.in max_elev_$file_suffix.in");
system("sed -i 's/000 000 !start and end stacks/$start_stack $end_stack !start and end stacks/g' max_elev_$file_suffix.in");
# system("sed -i 's/00 00 !start rec/$start_rec $end_rec !start rec/g' max_elev_$file_suffix.in");

$thisdir=cwd();


system("ln -sf $rundir/hgrid* .");
system("ln -sf $rundir/vgrid* .");
system("cp hgrid* vgrid.in read_output10_allnodes_viz.exe max_elev_$file_suffix.in $rundir/outputs/");
chdir("$rundir/outputs/");
system("./read_output10_allnodes_viz.exe < max_elev_$file_suffix.in");


chdir($thisdir);
system("cp $rundir/outputs/elevation_max* .");

system("~/bin/linear_comb_2grids < linear.in ");
system("mv elevation_max.gr3 elev_max_$file_suffix.gr3");
system("mv maxinun.gr3 maxinun_$file_suffix.gr3");
system("~/bin/grd2sms.pl maxinun_$file_suffix.gr3 maxinun_$file_suffix.2dm");
system("readlink -f maxinun_$file_suffix.2dm");



