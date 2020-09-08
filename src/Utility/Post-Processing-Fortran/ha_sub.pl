#!/usr/bin/perl -w

use Cwd;

# Called by ha_schism.pl to do extraction and harmonic analysis

if (@ARGV != 7) {
    print "Usage: $0 [full path to read_output8_allnodes_simple] [full path to tidal_analysis] 
          [full path to tidal_const.dat] [task ID] [start stack] [end stack] [# of pts]\n";
    print "Example: $0 /sciclone/home10/yinglong/bin/read_output8_allnodes_simple /sciclone/home10/yinglong/bin/tidal_analysis 
           /sciclone/home10/yinglong/tidal_const.dat filter_flag myout 2 4 43\n";
    exit(1);
}

$extract_code=$ARGV[0];
$ha_code=$ARGV[1];
$tidal_const=$ARGV[2];
$task_id=$ARGV[3]; #1-based
$start=$ARGV[4];
$end=$ARGV[5];
$npts=$ARGV[6];

$curr_dir=getcwd;
#Extract time series
$t=sprintf('%03d',$task_id);
system "rm -f outputs/in_$t";
system "cp -L $curr_dir/filter_flag_$t outputs/";
open(IN,">outputs/in_$t");
$cmd="1\n 2\n 1\n elev\n $start $end\n myout_$t\n filter_flag_$t\n";
print IN "$cmd";
close(IN);
system "rm -f outputs/myout_$t";
system "cd outputs; $extract_code < in_$t >& fout_$t; cd ../";
#Output is myout_$t

#HA
system "rm -f M2_$t K1_$t";
open(AP1,">M2_$t");
open(AP2,">K1_$t");
for($i=1; $i<=$npts; $i++)
{
  $j=$i+1;
  $cmd="awk \'{print \$1,\$$j}\' myout_$t > node_$t.out";
  system "cd outputs; $cmd; rm -f ap_$t; cd ../";
  system "cd outputs; $ha_code node_$t.out $tidal_const ap_$t; cd ../";
  open(AP,"<outputs/ap_$t") or die "Cannot open ap_$t\n";
  @ap=<AP>; close(AP);
  foreach $b (@ap)
  {
    ($name,$_,$_)=split(" ",$b);
    chomp $b;
    if($b =~ "M2" && length $name==2) {print AP1 "$i $b\n";} 
    elsif($b =~ "K1" && length $name==2) {print AP2 "$i $b\n";}
  }
}#for
close(AP1);
close(AP2);

#Output a marker for main perl
print "Done ha_sub.pl\n";
