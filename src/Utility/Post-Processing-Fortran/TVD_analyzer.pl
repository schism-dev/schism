#! /usr/bin/perl -w

#Read in outputs/nonfatal* and output violating elements 
#for a given time step threshold (Courant condition in Upwind or TVD).
#Run this in the dir where hgrid.gr3 is.
#Inputs: hgrid.gr3; outputs/nonfal*; screen
#Outputs: 
#        (1) TVD_upwind.txt (centroids of all violating elements; add header to become a bp file)

if (@ARGV != 4)
{ print "USAGE: $0 <# of CPUs> <dt allowed (sec)> <start time (days)> <end time (days)>\n"; exit(1); }

$cpu=$ARGV[0]; 
$dtb=$ARGV[1]; #time step threshold for reporting
$start=$ARGV[2];
$end=$ARGV[3];

open(IN,"<hgrid.gr3"); @all=<IN>; close(IN);
($ne,$np)=split(" ",$all[1]);
for($j=1;$j<=$ne;$j++) {$flag[$j]=0;}
for($i=1;$i<=$np;$i++)
{
  ($_,$x[$i],$y[$i],$dp)=split(" ",$all[$i+1]);
} #for i (node)
#print "Last node: $x[$np] $y[$np]\n";

for($i=1;$i<=$ne;$i++)
{
  ($_,$_,$nm[$i][1],$nm[$i][2],$nm[$i][3])=split(" ",$all[$i+1+$np]);
} #for j (element)
#print "Last elem: $nm[$ne][1] $nm[$ne][2] $nm[$ne][3]\n";

for($i=0;$i<$cpu;$i++)
{
  $id=sprintf("%04d",$i);
  open(IN,"<outputs/nonfatal_$id"); @all=<IN>; close(IN);
  #print "doing outputs/nonfatal_$id\n";
  $line=0;
  foreach $a (@all)
  {
    if($a =~ "TVD-upwind dtb info:")
    {
      #b[0..2]: TVD-upwind dtb info
      #b[3]: it - global time iteration #
      #b[4]: it_sub - subcylcing iteration # for transport
      #b[5]: global element # that has most strict Courant #
      #b[6]: vertical index for the prism that has most strict Courant #
      #b[7]: tracer index # (e.g., 1 for T, 2 for S)
      #b[8]: time step for transport
      #b[9]: current time in sec
      @b=split(" ",$a);
      #print "$a\n$b[3]\n";
      if($b[8]<$dtb && $line+1<=@all+1 && $b[9]/86400>=$start && $b[9]/86400<=$end) 
      {
        #@c=split(" ",$all[$line+1]); 
        $flag[$b[5]]=1;
      } #$b
    } #if $a
    $line++;
  } #foreach
} #for i

open(OUT,">TVD_upwind.txt");
print OUT "@ARGV[0..3]\n";
$bp=0; # # of pts
for($i=1;$i<=$ne;$i++) 
{
  if($flag[$i]==1) 
  {
    $bp++;
    $xc=0; $yc=0;
    for($j=1;$j<=3;$j++)
    {
      $nd=$nm[$i][$j];
      $xc=$xc+$x[$nd]/3;
      $yc=$yc+$y[$nd]/3;
    } #j
    print OUT "$bp $xc $yc $i\n";
  } #if $flag
} #for $i
close(OUT);
