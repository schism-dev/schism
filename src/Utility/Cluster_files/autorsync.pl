#! /usr/bin/perl -w

# This script can be run before a parallel SCHISM run (make sure purge outputs/ first), 
# to automatically rsync outputs back to a remote system (e.g. Sciclone), assuming 
# passwdless ssh is set up. It will also create separate dirs for each output stack to reduce file count. 
# It runs as a daemon. 
# Run on the run dir (i.e., one level above outputs/) on any system.  
# Assume that the code will output *(end_stack+1).nc (otherwise the script will hang).
#
#use Cwd;
#$pwd=cwd();

if(@ARGV != 3) 
{ 
  print "$0 <start # of stacks> <end # of stacks> <remote dir ending with outputs>\n";
  print "e.g.: $0 1 40 /sciclone/home10/yinglong/schism10/RUN02a/outputs\n";
  exit(1);
}
print "$0 @ARGV\n";
$start_stack=$ARGV[0]; $end_stack=$ARGV[1]; 
$remote_dir=$ARGV[2];

#Remote account
$chinook='yinglong@chinook.sciclone.wm.edu:';

if(!-e "outputs") {die "No outputs dir!";}

while(!-e "./outputs/local_to_global_0000") {sleep 20;} #sec
#sleep 10;

$count=0;
for($next_stack=$start_stack+1; $next_stack<=$end_stack+1; $next_stack++)
{
  $count=$count+1;
  while(!-e "outputs/schout_0000_$next_stack\.nc") {
    system "touch core.dummy";
    system "rm -rf core.*";
    sleep 20;
  } 
  system "touch core.dummy";
  system "rm -rf core.*";
  if($count==1) {
   system "rsync -az ./outputs/local* $chinook$remote_dir";
  }

  sleep 180; #wait a little longer to make sure all rank outputs are finished
  $current_stack=$next_stack-1;
  $cmd="rm -rf outputs.$current_stack; mkdir outputs.$current_stack; mv -f ./outputs/schout_????_$current_stack\.nc outputs.$current_stack";
  system "$cmd";
  $cmd="rsync -az ./outputs.$current_stack/schout_????_$current_stack\.nc $chinook$remote_dir/";
  print "$cmd\n";
  system "$cmd";
#For autocomb to start
#  $cmd="rsync -az ./outputs/schout_0000_$next_stack\.nc $chinook$remote_dir/";
#  print "$cmd\n";
#  system "$cmd";
  #system "rsync -av ./outputs/hot* $chinook$remote_dir/ ";
#  system "rsync -av ./outputs/mirror.out $chinook$remote_dir/mirror.out.$current_stack ";
#  system "rsync -av ./outputs/subcycling.out $chinook$remote_dir/subcycling.out.$current_stack ";
  print "done rsync'ing stack $current_stack...\n";
} #for

#Final sync
system "rsync -az ./outputs/  $chinook$remote_dir/ ";
#system "touch core.dummy";
system "rm -rf core.*";
#chdir("./outputs/");
#system "sbatch run_comb";
# system "mv hotstart_i* $pwd ";
#if(-e "outputs/$end_stack\_hvel.64") {system "cd outputs/; echo schout_[0-9]???_*.nc | xargs rm -f";}
