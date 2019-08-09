#!/usr/bin/perl
#
# Convert hgrid.gr3 or standard grid format to landxml
#

if (@ARGV != 2) {
    die "Usage: $0 [input hgrid.gr3] [output TIN in LandXML .xml]";
}

my ($np, $ne);

my $hgrid = $ARGV[0];
my $xmlgrid = $ARGV[1];

open(INP, "<$hgrid") or die "Unable to open file $hgrid for reading";
open(OUTP, ">$xmlgrid") or die "Unable to open file $xmlgrid for writing";

my $header = <<EOF;
<?xml version="1.0"?>
<LandXML xmlns="http://www.landxml.org/schema/LandXML-1.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.landxml.org/schema/LandXML-1.0 http://www.landxml.org/schema/LandXML-1.0/LandXML-1.0.xsd" version="1.0" date="2011-08-18" time="14:16:29" readOnly="false" language="English">
	<Project name="hgrid.gr3 conversion" desc=""/>
	<Units>
		<Metric linearUnit="meter" areaUnit="squareMetert" volumeUnit="cubicMeter" temperatureUnit="celsius" pressureUnit="milliBars" angularUnit="decimal degrees" directionUnit="decimal degrees"/>
	</Units>
	<Application name="grid2landxml.pl" manufacturer="CMOP" version="1.0" desc="pre-release build" manufacturerURL="www.stccmop.org" timeStamp="2011-08-18T09:16:29">
		<Author createdBy="grid2landxml.pl" createdByEmail="selfe@stccmop.org" company="CMOP" companyURL="www.stccmop.org" timeStamp="2011-07-29T11:16:29"/>
	</Application>
	<Surfaces>
		<Surface name="hgrid" desc="Computational grid">
			<Definition surfType="TIN" elevMax="0" elevMin="-8000">
				<Pnts>
EOF

print OUTP $header;

my $buf = <INP>;
$buf = <INP>;
($ne, $np) = split(' ', $buf);

for (my $i=0;$i<$np;$i++) {
    $buf = <INP>;
# <P id="1">495828.300000 220076.500000 3.0000000e+01</P>
    ($n, $x, $y, $d) = split(' ', $buf);
 
    my $nn = $i+1;
    print OUTP "<P id=\"$nn\">$y $x $d</P>\n";
}
print OUTP "</Pnts>\n";
print OUTP "<FACES>\n";
for (my $i=0;$i<$ne;$i++) {
    $buf = <INP>;
    ($n, $jnk, $n1, $n2, $n3) = split(' ', $buf);
    print OUTP "<F>$n1 $n2 $n3</F>\n";
}

my $footer = <<EOF;
</FACES>
                        </Definition>
                </Surface>
        </Surfaces>
</LandXML>
EOF
print OUTP $footer;
