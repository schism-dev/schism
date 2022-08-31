#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;

print(">>>>>>>>>>You may need to recompile 'gen_vqs.f90',>>>>>>>>>>>>>>\n");
print(">>>>>>>>>>see an example compiling cmd in the source code.>>>>>>>>>>>>>>\n");
#{
#  local( $| ) = ( 1 );
#  print "Press <Enter> or <Return> to continue: \n";
#  my $resp = <STDIN>;
#}

#$thisdir=cwd();

#UTM grid
system("ln -sf ../hgrid.utm.* hgrid.gr3");

#generate vgrid
system("./gen_vqs < in");
print(">>>>>>>>>>outputing on a transect>>>>>>>>>>>>>>\n");

print(">>>>>>>>>>removing old vgrid.in from run dir>>>>>>>>>>>>>>\n");
unlink("../vgrid.in");

print(">>>>>>>>>>copying vgrid.in to the run dir>>>>>>>>>>>>>>\n");
copy("vgrid.in","../vgrid.in");

print(">>>>>>>>>>Done.>>>>>>>>>>>>>>\n")
