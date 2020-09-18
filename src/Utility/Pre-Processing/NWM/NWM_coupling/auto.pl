#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

$thisdir=cwd();

print("You probably need an interactive session if you're on a cluster\n\n");

print("You may need to recompile 'coupling_nwm.f90' upon first use.\n");
print("e.g: ifort -O2 -CB -g -traceback -o coupling_nwm  ~/git/schism/src/Utility/UtilLib/julian_date.f90 ~/git/schism/src/Utility/UtilLib/schism_geometry.f90 ~/git/schism/src/Utility/UtilLib/pt_in_poly_test.f90 coupling_nwm.f90 -I\$NETCDF/include -I\$NETCDF_FORTRAN/ include -L\$NETCDF_FORTRAN/lib -L\$NETCDF/lib -L\$NETCDF/lib -lnetcdf -lnetcdff \n\n");

print("You may need to recompile 'combine_sink_source.f90' upon first use.\n");
print("e.g: ifort -CB -O2 -qopenmp -o combine_sink_source combine_sink_source.F90\n\n");

#grid
system("ln -sf ../hgrid.ll hgrid.ll");

#generate source/sink; note there're some cmd line inputs required
#receive cmd line input arguments for coupling_nwm
print "Days needed:\n";
chomp(my $nday = <STDIN>);
print "Starting date (dd mm yyyy):\n";
chomp(my $start_d = <STDIN>);
#write to couple.in
open(OUT,">couple.in");
print OUT "1e-3\n $nday\n $start_d\n";
close(OUT);
#link NWM outputs
print "NWM outputs dir:\n";
chomp(my $nwm_dir = <STDIN>);
system("ln -sf $nwm_dir NWM_DATA");
#run the fortran script coupling_nwm
system("./coupling_nwm < couple.in");

#combine nearby source/sink
system("ln -sf ../hgrid.cpp .");
system("./combine_sink_source < comb.in");
chdir '..';
system("ln -sf $thisdir/source_sink.in .");
system("ln -sf $thisdir/msource.th .");
system("ln -sf $thisdir/vsource.th.1 ./vsource.th");
system("ln -sf $thisdir/vsink.th.1 ./vsink.th");

print("Done.\n")
