#!/usr/bin/perl
#
# Generate images in batch mode for outputs by MPI SELFE dir. structure; by Paul T. and modified by Joseph Z.
# If gifsicle is not available, only individual .gif frames are generated.
# Inputs: (1) hgrid.gr3 (in $rundir); 
#         (2) binary files in outputs/ (*elev.61 etc);  
#         (3) xmvis par file (for appearance; get a snapshot you like and then Files-->Save parameters)
#         May edit "skip" below; also adjust "-gifp" to avoid distortion of images.
# Outputs: anim_var.gif (animated gif) or individual frames.

if (@ARGV != 8) {
    print "Usage: $0 [run dir (where hgrid.gr3 is)] [vis6 par file] 
    [variable name (e.g. elev.61] [vertical level # (must=1 for 2D vars)] [start stack #] [end stack #]
    [# of records per file] [stride used in skipping records]\n";
    print "Example: $0 . frame.par salt.63 41 1 6 24 1\n";
    exit(1);
}

$rundir = $ARGV[0];
$par= $ARGV[1];
$var = $ARGV[2];
$lev = $ARGV[3];
$start = $ARGV[4];
$end= $ARGV[5];
$nsteps= $ARGV[6];
$stride= $ARGV[7];

system "rm -f var*.gif";

for ($i=$start;$i<=$end;$i++) {
#    $nsteps =24 ;#`/usr/local/ace/bin/getnsteps $rundir/outputs/$i\_$var`;
#    print "# of step in each file = $nsteps\n";
# batch file
    open(BATP, ">batch.par");
    print BATP "batch prefix \"var_\"\n";
    print BATP "batch run 1,$nsteps\n";
    close(BATP);
# read par file
    open(OUTP, ">read.par");
    print OUTP "with grid 0\n";
    print OUTP "read grid 0 \"$rundir/hgrid.gr3\"\n";
    print OUTP "grid boundary display on\n";
    print OUTP "read ELCIRC 0 \"$rundir/outputs/$i\_$var\" start 0 stop $nsteps skip $stride level $lev\n";
    #Repeat to load in more variables if necessary
    #print OUTP "read ELCIRC 1 \"$rundir/outputs/$i\_$var2\" start 0 stop 360 skip 1 level $lev\n";
    close(OUTP);
    $comm="/usr/local/ace/bin/xmvis6 -batch batch.par -par read.par -par $par -gif -gifp 0,1200,0,600,0.45,5,10,10,10\n";
    print "doing stack $i...\n";
    system "$comm";
} #for

if(-e "/usr/local/bin/gifsicle" || -e "/usr/local/gifsicle-1.7.1/bin/gifsicle") {
  $comm="gifsicle -O2 -d40 var*.gif > anim_var.gif";
  print "$comm\n";
  system "$comm";
  #Clean up
  system "rm -f var*.gif";
} #if
