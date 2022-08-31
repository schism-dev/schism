#!/usr/bin/perl
#
# Convert points from epsg to lat/lon and vice versa (NAD83 only)
# on any system that has proj/invproj
# NOTE: This code does not do datum conversions like NAD27 to NAD83
#
# pturner July 2006
# Modified by Joseph Zhang

# -I +init=nad27:3601 +units=m -f %.5f -W6

if(@ARGV != 7) {
  print "Usage: $0 [EPSG or SPCS code] [1: proj to ll; 2: ll to proj] 
           [1: input format is .gr3; 2: input xyz] [input .gr3/xyz file] 
           [output] [for xyz input, separator (1: space; 2: comma)] 
           [for xyz input, line number (0: without; 1: with)]\n";
  print "Example: (a) $0 epsg:26910 1 2 tmp.xyz ll.nxyz 1 0\n 
           (b) $0 nad83:1802 (ME state plane coord.) 1 2 tmp.xyz ll.nxyz 1 0\n
           For xyz input, output with space separator and with line numbers\n
           For .gr3 input, outputs is .gr3 format\n";
  exit(1);
}
else {
  print "$0 @ARGV\n";
}

$epsg = $ARGV[0];
$inv=$ARGV[1]; #1: proj to ll; 2: ll to proj
$form=$ARGV[2]; #1: input .gr3; 2: input xyz
$gridin = $ARGV[3]; #input .gr3/xyz file name
$llgrid = $ARGV[4]; #output name
$deli= $ARGV[5]; #for xyz input, separator (1: space; 2: comma)
$noden= $ARGV[6]; #for xyz input, line number (0: without; 1: with)

open(INP, "<$gridin"); @grid = <INP>; close(INP);
open(OUTP, ">$llgrid");
print "done reading input file...\n";

if($form==1) { #.gr3 input
  print OUTP $grid[0];
  print OUTP $grid[1];
  ($nmel, $nmnp) = split(' ', $grid[1]);
} 
else {$nmnp = @grid;}
 
for($i=0;$i<$nmnp;$i++) {
  if($form==1) { #.gr3 
    ($n,$x,$y,$d[$i])=split(' ',$grid[$i+2]);
  }
  else { #.xyz input
    if($noden==0) { #no line numbers
      if($deli==1) {($x, $y, $d[$i]) = split(' ', $grid[$i]);}
      else {($x, $y, $d[$i]) = split(',', $grid[$i]);}
    }
    else { #has line numbers
      if($deli==1) {($n, $x, $y, $d[$i]) = split(' ', $grid[$i]);}
      else {($n, $x, $y, $d[$i]) = split(',', $grid[$i]);}
    }
  } #$form

  $xyin[$i] = "$x $y\n";
} #for
print "done preparing x,y for proj...\n";

system "rm -f tmp tmp2";
open(TMP,">tmp"); print TMP @xyin; close(TMP);
print "done writing x,y file for proj...\n";

# Convert
if($inv==1) { #to ll
  $cmd="invproj +init=$epsg +datum=NAD83 +units=m -f %.8f tmp >& tmp2";
}
else { 
  $cmd="proj +init=$epsg +datum=NAD83 +units=m -f %.8f tmp >& tmp2";
}
print "$cmd\n";
system $cmd;
open(TMP,"<tmp2"); @xyout=<TMP>; close(TMP);
system "rm -f tmp tmp2";
print "done reading outputs from proj...\n";

# Output pts
for($i=0;$i<$nmnp;$i++) {
#  $buf = sprintf("%d %.8f %.8f %.3f\n", $i+1, $xyout[$i], $lat, $d);
  $i1=$i+1;
  @xyout2=split(' ',$xyout[$i]);
  print OUTP "$i1 @xyout2 $d[$i]\n";
} #for

# Output table for .gr3
if($form==1) {print OUTP @grid[$nmnp+2..$nmnp+$nmel+1];}

close(OUTP);
