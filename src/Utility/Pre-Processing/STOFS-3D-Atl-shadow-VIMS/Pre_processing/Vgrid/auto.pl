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

#generate vgrid
print(">>>>>>>>>>generate vgrid>>>>>>>>>>>>>>\n");
system("./gen_vqs < in");

copy("vgrid.in","vgrid.in.old");

print(">>>>>>>>>>change format>>>>>>>>>>>>>>\n");
system("./change_vgrid < in");
copy("vgrid.in.new","vgrid.in");

print(">>>>>>>>>>Done.>>>>>>>>>>>>>>\n")
