#! /usr/bin/perl -w

# This script can be run during or after a parallel SCHISM run (e.g., using mpiexec or qsub);
# launch it shortly after the SCHISM run starts to make sure outputs/local_to_global_* are there (for
# the same reason, purge old outputs/* before launching this script to get right # of CPUs etc).
# It is useful when you don't know the strucutre of the cluster.
# It runs as a daemon and will combine outputs along the way.
# Run on the run dir (i.e., one level above outputs/) on any system (but
# make sure the combine code below is compatible).
# Assume that the code will output *(end_stack+1).nc (otherwise the script will hang).

#use Cwd;

if(@ARGV < 2) 
{ 
  print "$0 <start # of stacks> <end # of stacks> <optional wet/dry (0: last wet values)>\n";
  exit(1);
}
print "$0 @ARGV\n";
$start_stack=$ARGV[0]; $end_stack=$ARGV[1]; 
$iwetdry=0;
if(@ARGV==3) {$iwetdry=$ARGV[2];}

if(!-e "outputs") {die "No outputs dir!";}
##Get nproc and ntracers
##open(FL,"<outputs/local_to_global_0000"); @all=<FL>; close(FL);
##($_,$_,$_,$nproc,@ntr[0..9])=split(" ",$all[0]);

$code="~yinglong/bin/combine_output11";

for($next_stack=$start_stack+1; $next_stack<=$end_stack+1; $next_stack++)
{
  while(!-e "outputs/schout_0000_$next_stack\.nc") {sleep 120;} #sleep 2 min.
#  sleep 180; #wait a little longer to make sure outputs are written
  $current_stack=$next_stack-1;
  system "cd outputs; $code -b $current_stack -e $current_stack -w $iwetdry > combine.out";
  print "done combining stack $current_stack...\n";
} #for

#if(-e "outputs/$end_stack\_hvel.64") {system "cd outputs/; echo schout_[0-9]???_*.nc | xargs rm -f";}
