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


#------------set tvd.prop--------------------
@reg_list = glob("tvd/*.r*");
$list_size = @reg_list;
if ($list_size == 0){
    print("No regions are provided in Prop/tvd/.\n");
    print("Inside the provided regions, high-order transport schemes are NOT used.\n");
    print("Do you want to set 1 for all elements in tvd.prop (use higher-orde transport schemes for the whole domain)?\n");
    print("y or n\n");
    chomp ($_=<STDIN>);
    if( /^[Yy](?:es)?$/ ) { # set tvd.prop=1 uniformly
        if ( ! open inFile, "< hgrid.gr3" ){
            die "can not open hgird.gr3: $!";
        }
        if ( ! open outFile, "> tvd.prop" ){
            die "can not open tvd.prop: $!";
        } 
        #blank line 
        $line=<inFile>;
        #get ne,np
        $line=<inFile>;
        $line =~ s/^\s+|\s+$//g;
        ($ne,$np) = split(/ +/, $line);
        print("ne, np: $ne, $np\n");
        # write tvd.prop
        for (my $i=1; $i <= $ne; $i++) {
            print outFile "$i 1\n";
        }
    }else{
        print("Put *.reg files in Prop/tvd/, then try again\n");
        print("See samples in Prop/Sample_tvd/\n");
        print("Aborting ...\n");
        exit
    }
}else{ #regions provided
    $dummy_reg = $reg_list[0];
  # initialize to 1 (i.e., use tvd)
    system("$script_dir/auto_edit_prop $dummy_reg hgrid.ll 1 1"); 
    move("out.prop","default.prop");
  # set local regions to 0
    foreach (@reg_list) {
        $this_reg = $_;
        print("$script_dir/auto_edit_prop $this_reg hgrid.ll 0 -9999\n");
        system("$script_dir/auto_edit_prop $this_reg hgrid.ll 0 -9999");
        move("out.prop","default.prop");
    }
    move("default.prop","tvd.prop");
    unlink("../tvd.prop");
    copy("tvd.prop","../tvd.prop");
}

#-----------------------fluxflag-------------------
@reg_list = glob("fluxflag/*+.reg");
$list_size = @reg_list;
if ($list_size == 0){
  print("Fluxflag regions not provided. See samples in Prop/Sample_fluxflag_*/\n");
  print("If this is intentional, you should turn off 'iflux' in param.nml instead.\n");
  print("Aborting ...\n");
  exit
}
$dummy_reg = $reg_list[0];
#initialize background as -1
system("$script_dir/auto_edit_prop $dummy_reg hgrid.ll -1 -1");
move("out.prop","default.prop");

#dump the table of transects into fluxflag.prop.table
if ( ! open outFile, "> fluxflag.prop.table" ){
    die "can not open fluxflag.prop.table: $!";
} 
my $n = 0;
# loop through each transect
foreach (@reg_list) {
    $m=$n+1;
    # only upstream regs (+) are in the array
    $upstream_reg = $_;
    # replace + with - to get the downstream reg's name
    $downstream_reg = $upstream_reg;
    $downstream_reg =~ s/\+/\-/g;

    print outFile "$upstream_reg: $m; $downstream_reg: $n\n";

    # downstream
    print("$script_dir/auto_edit_prop $downstream_reg hgrid.ll $n -9999\n");
    system("$script_dir/auto_edit_prop $downstream_reg hgrid.ll $n -9999");
    move("out.prop","default.prop");

    # upstream
    print("$script_dir/auto_edit_prop $upstream_reg hgrid.ll $m -9999\n");
    system("$script_dir/auto_edit_prop $upstream_reg hgrid.ll $m -9999");
    move("out.prop","default.prop");

    $n = $n+1;
}
move("default.prop","fluxflag.prop");

unlink("../fluxflag.prop");
copy("fluxflag.prop","../fluxflag.prop");


print("Done.\n")
