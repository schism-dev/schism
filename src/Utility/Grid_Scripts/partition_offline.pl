#!/usr/bin/perl -w

# Prepare partition.prop for offline mode (no ParMETIS). Drives a FORTRAN and a C code (gpmetis)
# Run this script on any Sciclone/Ches system
# Inputs: hgrid.gr3 (with b.c. part); vgrid.in
# Outputs: partition.prop.<corecount> (copy or link it to partition.prop)
if ( @ARGV !=1) {
    print "Usage: $0 [core_count (compute only excluding scribes] \n";
    print "e.g.: $0 3334\n";
    exit(1);
}

$cores=$ARGV[0];

$prep="~yinglong/bin/metis_prep"; #FORTRAN script to generate graphinfo
$gpmetis="~yinglong/bin/gpmetis"; #C code of METIS

system "rm -f graphinfo; $prep";
system "rm -f graphinfo.part.$cores";
$cmd="$gpmetis graphinfo $cores -ufactor=1.01 -seed=15";
print "$cmd\n";
system $cmd;

system "awk \'{print NR,\$0}' graphinfo.part.$cores > partition.prop.$cores";

