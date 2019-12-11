#!/usr/bin/perl -w

# Load all .asc 
# Run in this dir
# Inputs: hgrid.old (UTM 18N, NAVD88)
# Output: revised hgrid.old (NGVD29); original input is renamed as hgrid.old.0

$code='./interpolate_depth_structured2';

system "cp -L hgrid.old hgrid.old.0"; #save a copy
unlink "out";

open(IN,">in");
print IN "-1\n -0.000\n";
close(IN);
for ($i=0; $i<=19; $i++)
{
  #if ($i==1 || $i==5) {
  $a='nj_de_tbdem_resamp3m_split'.$i.'.tif.asc';
  if (-e $a) {
    print "loading $a\n";
    system("ln -sf $a struc.grd");
    system("$code < in >> out; mv hgrid.new hgrid.old;");
  }
  #}
} #foreach
