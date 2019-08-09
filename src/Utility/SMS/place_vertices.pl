#!/usr/bin/perl
#   Copyright 2014 College of William and Mary
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


#Read in a map file (SMS 11.xx) and manipulate vertices along each arc
if (@ARGV !=2) {
	die "usage: $0 infile(.map) outfile\n";
}
$file = $ARGV[0];
$outfile = $ARGV[1];

open(IN,$file);
@lines = <IN>;
close(IN);
open(OUT,">$outfile");

$nline=@lines;
for($i=0;$i<$nline;$i++) 
{
  @a=split(" ",$lines[$i]);
  if($a[0] ~= "NODE")
  {
  }
  elsif($a[0] ~= "ARC")
} #for


print OUT "MESH2D\n";
chomp $lines[1];
$lines[1]=~s/^\s+//;
$lines[1]=~s/\s+/ /g;
($e,$n)=split(" ",$lines[1]);
$starte = $n+2;
print "$lines[1]\n$e $n $starte\n";
for ($i = $starte; $i<$n+$e+2; $i++){
chomp $lines[$i];
$lines[$i]=~s/^\s+//;
$lines[$i]=~s/\s+/ /g;
	($elemn,$elem34,$e1,$e2,$e3,$e4)=split(" ",$lines[$i]);
	if ($elem34 == 4) {print OUT "E$elem34"."Q $elemn $e1 $e2 $e3 $e4 1\n";}
	elsif ($elem34 == 3) {print OUT "E$elem34"."T $elemn $e1 $e2 $e3 1\n";}
}
for ($i = 2; $i<$starte; $i++){
chomp $lines[$i];
$lines[$i]=~s/^\s+//;
$lines[$i]=~s/\s+/ /g;
	print OUT "ND $lines[$i]\n";
}
