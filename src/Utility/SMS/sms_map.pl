#!/usr/bin/perl
#Manipulate x,y coordinates in .map files. Example here for CPP proj to lon/lat

if (@ARGV !=4) {
	die "usage: $0 infile (.map) outfile(.map) CPP_lon CPP_lat\n";
}
$file = $ARGV[0];
$outfile = $ARGV[1];
$lam0=$ARGV[2];
$phi0=$ARGV[3];

open(IN,"<$file"); @all=<IN>; close(IN);
open(OUT,">$outfile");

$linenum=0;
while ($linenum<@all) {
  @char=split(" ",$all[$linenum]);
  if($char[0] =~ "XY") {
    @char2[0..1]=projection($lam0,$phi0,$char[1],$char[2]);
    print OUT "XY @char2 0.0\n";
  } 
  elsif($char[0] =~ "ARCVERTICES") {
    $np=$char[1];
    print OUT $all[$linenum];
    for($i=1;$i<=$np;$i++) {
      @char=split(" ",$all[$linenum+$i]);
      @char2[0..1]=projection($lam0,$phi0,$char[0],$char[1]);
      print OUT "@char2 0\n";
    } #for
    $linenum=$linenum+$np;
  } 
  else {
    print OUT $all[$linenum];
  }#if
  
  $linenum++;
} #while
# =~ "ARCVERTICES" #followed by np \n x y 0 ....
close(OUT); 

#Modify this routine as you wish for reproj
sub projection {
  ($lam00,$phi00,$x,$y)=@_;
  $rearth=6378206.4;
  $pi=3.1415926;

  $lam=$lam00/180*$pi;
  $phi=$phi00/180*$pi;
#  $lam0=-124.35/180*$pi; $phi0=43.36/180*$pi;
  $xl=$lam+$x/$rearth/cos($phi);
  $yl=$y/$rearth;
  $xout=$xl/$pi*180;
  $yout=$yl/$pi*180;
#  $lam0=-124.498333333/180*$pi; $phi0=42.73833333333/180*$pi;
#$xout=$rearth*($xl-$lam0)*cos($phi0);
#$yout=$yl*$rearth;
  return ($xout,$yout);
} #sub
