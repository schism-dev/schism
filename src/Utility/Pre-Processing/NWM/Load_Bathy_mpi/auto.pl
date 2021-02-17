#! /usr/bin/perl -w

#load some functions
use POSIX ();
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;
use Cwd 'abs_path';

$pwd = cwd();

#-------------------inputs------------------------------
my $DEM_dir = "~/work/stampede2/NWM/DEM/DEM/";
# sample on james: $DEM_dir = "/ches/data10/whuang07/Case1/DEMs/DEM_pre/DEM/";
# sample on sciclone: $DEM_dir = "/sciclone/data10/whuang07/NWM/DEM/DEM/";
# sample on stampede2: $DEM_dir = "~/work/stampede2/NWM/DEM/DEM/";
#-------------------end inputs--------------------------

if ($DEM_dir eq "") {
    print(">>> Source DEM folder needs to be specified in the 'input' section of this script\n");
    print(">>> Contact the VIMS group if you're not using W&M/VIMS HPCs or Stampede2.\n");
    print(">>> Aborting ...\n");
    exit;
}

# copy hgrid
system('rm hgrid.old*');
system('cp ./hgrid.ll ./hgrid.old');

while (!(-e 'hgrid.old')) {
    print("hgrid.old not found\n");
    print(">>> Aborting ...\n");
    exit;
    #print("fix the issue then ");
    #WaitForKey();
}

#find out how many cores per node on this cluster
$host=`hostname`;
print("running on $host\n");

if ($host =~ /cyclops/) {
    $cores_per_node = 24;
    $nodes_bloat = 2.75;  #use more nodes to secure enough memory
    $nodes_min = POSIX::ceil($ndems/$cores_per_node);
    print "Minimal number of nodes for Batch: $nodes_min\n";
    $nodes = POSIX::ceil($ndems/$cores_per_node*$nodes_bloat);
    print "Number of nodes actually requested for Batch: $nodes\n";
}elsif ($host =~ /(james)/){
    $cores_per_node = `pbsnodes | grep -A 3 jm01 | grep np`;
    $nodes_bloat = 1.0;  #use more nodes to secure enough memory
    $nodes = 10;
    print "Number of nodes actually requested for Batch: $nodes\n";
}elsif ($host =~ /chesapeake/){
    $cores_per_node = `pbsnodes | grep -A 3 pt01 | grep np`;
    $nodes_bloat = 1.1;  #use more nodes to secure enough memory
    $nodes_min = POSIX::ceil($ndems/$cores_per_node);
    print "Minimal number of nodes for Batch: $nodes_min\n";
    $nodes = POSIX::ceil($ndems/$cores_per_node*$nodes_bloat);
    print "Number of nodes actually requested for Batch: $nodes\n";
}elsif($host =~ /stampede2/) {
    #$cores_per_socket=`scontrol show node  | grep -B 9  skx-dev | head -n 10 | grep CoresPerSocket`
    #$n_sockets=`scontrol show node  | grep -B 9  skx-dev | head -n 10 | grep Sockets`
    #the above goes too far, but will probably be useful in some cases.
    $cores_per_node=48;
    #print("Automation on Stampede2 not yet tested\n");
    #print("Aborting ...\n");
    #exit;
    $nodes_bloat = 1.0;  #use more nodes to secure enough memory
    $nodes = 4;
    print "Number of nodes actually requested for Batch: $nodes\n";
}else{
    print("The script has not been tested on $host\n");
    print("Aborting ...\n");
    exit;
}

$cores_per_node=~s/[^0-9]//g;
print "cores_per_node=$cores_per_node\n";

# !!!not copying symlink*pl,
# which is a temporary fix before the dem_* naming is fixed
#prepare DEM folder
#system("rm -rf ./DEM/");
#system("mkdir DEM");
chdir("DEM");
system("rm -rf ./DEM/*.asc");

# system("ln -sf $DEM_dir/*.asc .");
# link absolute path
(@asc_files)=glob "$DEM_dir/*.asc";
foreach (@asc_files) {
    $file=$_;
    $fullPN=abs_path($file);
    system "ln -sf $fullPN .";
}
system "rm dem_[0-9][0-9][0-9][0-9].asc";
# system("cp $DEM_dir/*.pl .");
chdir($pwd);

#make a copy of the hgrid.old
system("cp hgrid.old hgrid.old.0");
system("rm hgrid.new");


#load coastal DEMs in several batches:
print "---------Doing DEMs Batch--------------\n";
print "head and tail of hgrid.old\n";
system("head -n 2 hgrid.old");
system("tail -n 1 hgrid.old");

#link batch #i DEMs
chdir("DEM");
print("Calling ./symlink_dems.pl\n");
system("./symlink_dems.pl > /dev/null");
chdir($pwd);
system("rm dem*.asc");
system("rm *.out");
system("ln -sf DEM/dem_????.asc .");

#count number of dems
(@dems)=glob "dem*.asc";
$ndems = @dems;
print "Number of total DEMs tiles: $ndems\n";

#write dems.in
open(my $fh, '>', 'dems.in');
print $fh "$ndems !# of DEMs\n";
print $fh "$nodes !# of compute nodes for load balancing\n";
close $fh;

#------submit the jobs-----------------
system("rm run_submit");
if($host =~ /chesapeake/){
  open PW, "> ./run_submit";
  print PW "#!/usr/bin/tcsh\n";
  print PW "module list\n";
  print PW "cd $pwd/\n";
  print PW "mvp2run -v -C 0.05 ./interpolate_depth_structured2_mpi > err2.out\n";
  system "qsub -N Load_Bathy run_submit -l nodes=$nodes:c18:ppn=$cores_per_node  -l walltime=00:50:00";
}elsif ($host =~ /(james)/){
  open PW, "> ./run_submit";
  print PW "#!/usr/bin/tcsh\n";
  print PW "module purge\n";
  print PW "module load use.own schism2\n";
  print PW "module list\n";
  print PW "cd $pwd/\n";
  print PW "mvp2run -v -C 0.05 ./interpolate_depth_structured2_mpi > err2.out\n";
  system "qsub -N Load_Bathy run_submit -l nodes=$nodes:james:ppn=$cores_per_node  -l walltime=00:50:00";
}elsif ($host =~ /cyclops/) {
  open PW, "> ./run_submit";
  print PW "#!/usr/bin/tcsh\n";
  print PW "#SBATCH --nodes=$nodes --ntasks-per-node=$cores_per_node\n";
  print PW "#SBATCH --constraint=cyclops\n";
  print PW "#SBATCH -t 00:50:00\n";
  print PW "module list\n";
  print PW "cd $pwd/\n";
  print PW "srun ./interpolate_depth_structured2_mpi > err2.out\n";
  system "sbatch run_submit";
}elsif ($host =~ m/stampede2/) {
  $total_cores = $nodes*$cores_per_node;
  open PW, "> ./run_submit";
  print PW "#!/bin/bash\n";
  print PW "#SBATCH -J CORIE # Job name\n";
  print PW "#SBATCH -o err2.o%j       # Name of stdout output file\n";
  print PW "#SBATCH -e err2.e%j       # Name of stderr error file\n";
  print PW "#SBATCH -p skx-dev\n";
  print PW "#SBATCH -N $nodes \n";
  print PW "#SBATCH -n $total_cores\n";
  print PW "#SBATCH -t 1:00:00 \n";
  print PW "module list\n";
  print PW "pwd\n";
  print PW "date\n";
  print PW "ibrun ./interpolate_depth_structured2_mpi >& err2.out\n";
  system "sbatch run_submit";
}


#wait for hgrid.new
while(!-e "./hgrid.new") {sleep 30;}
print "hgrid.new generated.\n";

#wait another 180s to make sure hgrid.new is fully written
$nLines0 = `wc -l < hgrid.old.0 | bc`;

$nLines = `wc -l < hgrid.new | bc`;
while ($nLines0!=$nLines){
   sleep 10;
   $nLines = `wc -l < hgrid.new | bc`;
}
print "Head and tail of hgrid.new.\n";
system("head -n 2 hgrid.new");
system("tail -n 1 hgrid.new");

print "copying hgrid.new to hgrid.old.$i\n";
system("cp hgrid.new hgrid.old.1");
$nLines = `wc -l < hgrid.old.1 | bc`;
while ($nLines0!=$nLines){
  sleep 10;
  $nLines = `wc -l < hgrid.old.1 | bc`;
}
print "done copying hgrid.new to hgrid.old.1\n";

print "copying hgrid.new to hgrid.old\n";
system("cp hgrid.new hgrid.old");
$nLines = `wc -l < hgrid.old | bc`;
while ($nLines0!=$nLines){
  sleep 10;
  $nLines = `wc -l < hgrid.old | bc`;
}
print "done copying hgrid.new to hgrid.old\n";
system("rm hgrid.new");

system("cp hgrid.old hgrid.new");
system("rm hgrid.old");
system("rm *.out");
system("rm *.asc");
system("rm DEM/*.asc");
# system("rm -rf ./DEM/");

# set minmum depth to be -10m, minimum depth around Caribean Sea to be 5m
#system("./auto_edit_region 1 dummy.reg hgrid.new -10");
#system("cp out.gr3 hgrid.new; rm out.gr3"); # use cp and rm instead of mv, because mv has some weird problem on james
# system("./auto_edit_region 1 min_5m_ll.reg hgrid.new 5");
# system("cp out.gr3 hgrid.new; rm out.gr3");

@regions=("min_5m_ll.reg","BergenPoint.reg","SabinePass.reg");
@min_vals=(5,2,7);
for my $i (0..$#regions) {
    $region=$regions[$i];
    $min_val=$min_vals[$i];
    system("./auto_edit_region 1 $region hgrid.new $min_val");
    system("cp out.gr3 hgrid.new; rm out.gr3");
}

#sub WaitForKey() {
#    print "press any key to continue\n";
#    chomp($key = <STDIN>);
#}

