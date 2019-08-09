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





#
# $Revision: 1.4 $
#
# Convert SMS grid+boundary to .gr3 format.
#
# Author - cseaton 01-2005.
#
# $Log: 2dm2gr3.pl,v $
# Revision 1.4  2005/02/06 16:19:56  pturner
# Adjustments to 2dm2gr3.pl to account for case in river names.
#
# Revision 1.3  2005/01/16 20:21:33  pturner
# Working.
#
# Revision 1.2  2005/01/15 17:04:32  pturner
# Changes by cseaton to add support for rivers. Perltidy.
#
# Revision 1.1.1.1  2005/01/15 16:54:33  pturner
# Forecast NEFS support codes
#
#

if ( @ARGV < 3 ) {
    print "Usage: $0 [SMS version (0: before 9.; 1: 10.xx; 2: 11.xx)] 
    [SMS grid file in meter with boundary information] [.gr3 file in meters] [river1] ... [riverN]\n";
    print "River name is checked against boundary types.\n";
    exit(1);
}

$vers=$ARGV[0];
$file    = $ARGV[1];
$outfile = $ARGV[2];
@rivernames = @ARGV[ 3 .. $#ARGV ];

#Determine # of line to skip at the beginning
if($vers==0) #before 9.xx
  {$i=1;}   
elsif($vers==1) #v10.xx
  {$i=2;}   
elsif($vers==2) #v11.xx
  {$i=3;}   
else
  {die "Unknown SMS version $vers\n";}

print "$file $outfile\n";

open( IN,  "$file" )     || die "Can't open $file\n";
open( OUT, ">$outfile" ) || die "Can't open $outfile\n";

@lines = <IN>;
close(IN);
$elemflag = 'yes';
$elem     = '';
print "nlines = " . @lines . "\n";

while ( $elemflag eq 'yes' ) {
    $tmp = $lines[ $i++ ];

    #print $tmp;
    chomp $tmp;
    $tmp =~ s/\s+/ /g;
    $tmp =~ s/^\s//;
    if ( $tmp =~ /^E(\d)/ ) {
        $en       = $1;
        @tmp      = split ( " ", $tmp );
        $elem_num = $tmp[1];
        $elem .= "$elem_num $en";
        for ( $j = 0 ; $j < $en ; $j++ ) {
            $elem .= " $tmp[$j+2]";
        }
        $elem .= "\n";
    } else {
        $elemflag = 'no';
        --$i;
    }
}
$nodeflag = 'yes';
$node     = '';
while ( $nodeflag eq 'yes' ) {
    $tmp = $lines[ $i++ ];

    #print $tmp;
    chomp $tmp;
    $tmp =~ s/\s+/ /g;
    $tmp =~ s/^\s//;
    if ( $tmp =~ /^ND/ ) {
        @tmp = split ( " ", $tmp );
        $node_num = $tmp[1];
        $node .= "$node_num $tmp[2] $tmp[3] $tmp[4]\n";
    } else {
        $nodeflag = 'no';
        --$i;
    }
}
$boundaryflag        = 'yes';
@bound               = ();
@boundlen            = ();
@boundname           = ();
@boundtype           = ();
@boundname_ordered   = ();
$boundind            = 0;
$boundind_end        = 'no';
$bound[$boundind]    = '';
$boundlen[$boundind] = 0;

while ( $boundaryflag eq 'yes' ) {
    $tmp = $lines[ $i++ ];

    #print $tmp;
    chomp $tmp;
    $tmp =~ s/\s+/ /g;
    $tmp =~ s/^\s//;
    if ( $tmp =~ /^NS/ ) {
        @tmp = split ( " ", $tmp );
        if ( $tmp[$#tmp] < 0 ) {
            $boundind_end =
              'yes';    # mark that this will be the end of the boundary
            $tmp[$#tmp] *= -1;    #change sign
        }
        for ( $tmpind = 1 ; $tmpind < @tmp ; $tmpind++ ) {
            print "$tmpind $tmp[$tmpind]\n";
            $bound[$boundind] .= "$tmp[$tmpind]\n";
        }
        $boundlen[$boundind] += $#tmp;
        if ( $boundind_end eq 'yes' ) {
            $bound[ ++$boundind ] = '';
            $boundlen[$boundind]  = 0;
            $boundind_end         = 'no';
        }
    } else {
        $boundaryflag = 'no';
        --$i;
    }
}

while ( $lines[ ++$i ] !~ /BD 1/ && $i < @lines ) { print $lines[ $i - 1 ] }
while ( $lines[$i] !~ /BEG2DMBC/ && $i < @lines ) {
    $tmp = $lines[ $i++ ];
    if ( $tmp =~ /BD 1/ ) {
        chomp $tmp;
        $tmp =~ s/\s+/ /g;
        $tmp =~ s/^\s//;
        @tmp = split ( " ", $tmp );
        $tmp[2] =~ s/"//g;
        $boundname[ $tmp[3] - 1 ] = lc($tmp[2]);
        print "boundary type $tmp[3] assigned name of  $boundname[$tmp[3]-1]\n";
    }
}
while ( $lines[ $i++ ] !~ /BEG2DMBC/ && $i < @lines ) { print $lines[ $i - 1 ] }
$bound_index = 0;
while ( $lines[$i] =~ /BCS/ ) {
    $tmp = $lines[ $i++ ];
    print $tmp;
    chomp $tmp;
    $tmp =~ s/\s+/ /g;
    $tmp =~ s/^\s//;
    @tmp = split ( " ", $tmp );
    $boundname_ordered[ $bound_index ] = $boundname[ $tmp[2] - 1 ];
    print "boundary "
      . ( $bound_index )
      . " assigned a name of $boundname[$tmp[2]-1]\n";
    print "boundary type $tmp[2] is $boundname[$tmp[2]-1]\n";

    if ( $boundname[ $tmp[2] - 1 ] eq 'ocean' ) {
        $boundtype[ $bound_index ] = 'ocean';
        print "designating boundary $bound_index  as ocean\n";
    } elsif ( 'island' eq $boundname[ $tmp[2] - 1 ] ) {
        $boundtype[ $bound_index ] = 'island';
        print "designating boundary $bound_index as island\n";
    } elsif ( 'land' eq $boundname[ $tmp[2] - 1 ] ) {
        $boundtype[ $bound_index ] = 'land';
        print "designating boundary $bound_index as land\n";
    } else {
        $riverflag = 'no';
        foreach $river (@rivernames) {
            if ( $river eq $boundname[ $tmp[2] - 1 ] ) {
                $boundtype[ $bound_index ] = 'river';
                $riverflag = 'yes';
                print "designating boundary $bound_index as river\n";
            }
        }
        if ( $riverflag eq 'no' ) {
            $boundtype[ $bound_index ] = 'land';
            print "designating boundary $bound_index as dry river\n";
        }
    }
    ++$bound_index;

}
@boundname       = @boundname_ordered;
$ocean_bound     = '';
$river_bound     = '';
$open_bound_nb   = 0;
$open_bound_nn   = 0;
$land_bound      = '';
$island_bound    = '';
$closed_bound_nb = 0;
$closed_bound_nn = 0;
$bound_n         = 1;
$closed_bound    = '';
for ( $ind = 0 ; $ind < @boundtype ; $ind++ ) {
    print "$ind $boundtype[$ind] $boundname[$ind] ";
    if ( $boundtype[$ind] =~ /river/ ) {    # open boundary
        print "is a river\n";
        $open_bound_nb++;
        $river_bound .=
"$boundlen[$ind] = Number of nodes for open boundary $open_bound_nb $boundname[$ind]\n$bound[$ind]";
        $open_bound_nn += $boundlen[$ind];
    } elsif ( $boundtype[$ind] =~ /ocean/ ) {    # open boundary
        print "is the ocean\n";
        $open_bound_nb++;
        $ocean_bound .=
"$boundlen[$ind] = Number of nodes for ocean boundary $open_bound_nb\n$bound[$ind]";
        $open_bound_nn += $boundlen[$ind];
    } elsif ( $boundtype[$ind] =~ /island/ ) {
        print "is an island\n";
        $closed_bound_nb++;
        $island_bound .=
"$boundlen[$ind] = Number of nodes for island boundary $closed_bound_nb\n$bound[$ind]";
        $closed_bound_nn += $boundlen[$ind];
    } elsif ( $boundtype[$ind] =~ /land/ ) {     # land boundary
        print "is land\n";
        $closed_bound_nb++;
        $land_bound .=
"$boundlen[$ind] = Number of nodes for land boundary $closed_bound_nb\n$bound[$ind]";
        $closed_bound_nn += $boundlen[$ind];
    } else {
        die "unknown boundary type $boundtype[$ind] in $ind boundary";
    }
}
if ( $open_bound_nb + $closed_bound_nb > 0 ) {
    $bound =
      "$open_bound_nb = number of open boundaries\n"
      . "$open_bound_nn = total number of open boundary nodes\n"
      . $ocean_bound
      . $river_bound
      . "$closed_bound_nb = number of closed boundaries\n"
      . "$closed_bound_nn = total number of closed boundary nodes\n"
      . $land_bound
      . $island_bound;
}
print OUT "$file\n$elem_num $node_num\n$node$elem$bound";
close(OUT);
