#! /usr/bin/perl -w

#load some functions
use POSIX ();
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;
use Cwd 'abs_path';

$pwd = cwd();

#-------------------inputs------------------------------
my $DEM_dir = "";
# sample: $DEM_dir = '/ches/data10/whuang07/Case1/DEMs/DEM_pre/DEM/';
# sample: $DEM_dir = '/sciclone/data10/whuang07/NWM/DEM/DEM/';
#-------------------end inputs--------------------------

if ($DEM_dir eq "") {
    print(">>> Source DEM folder needs to be specified in the 'input' section of this script\n");
    print(">>> Contact the VIMS group if you're not using W&M/VIMS HPCs or Stampede2.\n");
    print(">>> Aborting ...\n");
    exit;
}

# copy hgrid
system('rm hgrid.old*');
system('cp ../hgrid.ll ./hgrid.old');

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
}elsif ($host =~ /(james)/){
    $cores_per_node = `pbsnodes | grep -A 3 jm01 | grep np`;
    $nodes_bloat = 1.3;  #use more nodes to secure enough memory
}elsif ($host =~ /chesapeake/){
    $cores_per_node = `pbsnodes | grep -A 3 pt01 | grep np`;
    $nodes_bloat = 1.1;  #use more nodes to secure enough memory
}elsif($host =~ /stampede2/) {
    #$cores_per_socket=`scontrol show node  | grep -B 9  skx-dev | head -n 10 | grep CoresPerSocket`
    #$n_sockets=`scontrol show node  | grep -B 9  skx-dev | head -n 10 | grep Sockets`
    #the above goes too far, but will probably be useful in some cases.
    $cores_per_node=48;
    print("Automation on Stampede2 not yet tested\n");
    print("Aborting ...\n");
    exit;
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

#load etopo and cut off at 5 m:
open(IN,">etopo.in");
print IN "-1\n -0.000\n";
close(IN);
print "loading etopo ...\n";
system("ln -sf DEM/etopo1.asc struc.grd");
system("./interpolate_depth_structured2 < etopo.in; mv hgrid.new hgrid.old");
print "cut off at 5 m ...\n";
# 5 m cut-off in/outside the dummy region, i.e., depth=max(5, depth) for all nodes
system("./auto_edit_region 1 dummy.reg hgrid.old 5 5");
system("cp out.gr3 hgrid.old.etopo_cutoff_5m");
system("mv out.gr3 hgrid.old");

#load coastal DEMs in several batches:
for ($i=1;$i<=3;$i++){
    print "---------Doing DEMs Batch #$i--------------\n";
    print "head and tail of hgrid.old\n";
    system("head -n 2 hgrid.old");
    system("tail -n 1 hgrid.old");

#link batch #i DEMs
    chdir("DEM");
    print("Calling ./symlink_dems_part$i.pl\n");
    system("./symlink_dems_part$i.pl > /dev/null");
    chdir($pwd);
    system("rm dem*.asc");
    system("rm *.out");
    system("ln -sf DEM/dem_????.asc .");

#count number of dems
    (@dems)=glob "dem*.asc";
    $ndems = @dems;
    print "Number of DEMs for Batch $i: $ndems\n";

#calculate how many nodes ($nodes) are required.
    $nodes_min = POSIX::ceil($ndems/$cores_per_node);
    print "Minimal number of nodes for Batch $i: $nodes_min\n";
    $nodes = POSIX::ceil($ndems/$cores_per_node*$nodes_bloat);
    print "Number of nodes actually requested for Batch $i: $nodes\n";

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
    }

#wait for hgrid.new
    while(!-e "./hgrid.new") {sleep 30;}
    print "hgrid.new generated.\n";

#wait another 180s to make sure hgrid.new is fully written
    sleep 180;
    print "Head and tail of hgrid.new.\n";
    system("head -n 2 hgrid.new");
    system("tail -n 1 hgrid.new");

    system("cp hgrid.new hgrid.old.$i");
    sleep 30;
    system("mv hgrid.new hgrid.old");
    sleep 30;
}
system("mv hgrid.old hgrid.new");
system("rm *.out");
system("rm *.asc");
system("rm DEM/*.asc");
# system("rm -rf ./DEM/");

 
#sub WaitForKey() {
#    print "press any key to continue\n";
#    chomp($key = <STDIN>);
#}

