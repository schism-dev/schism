#! /usr/bin/perl -w

$hgrid="hgrid.gr3";

@gr3_files=("fg.gr3");
@val=      ("33");

$iFile=0;

foreach (@gr3_files) {
    $gr3=$_;

    if ( ! open inFile, "< $hgrid" ){
        die "can not open $hgrid: $!";
    }
     
    unlink $gr3;
    if ( ! open outFile, "> $gr3" ){
        die "can not open $gr3: $!";
    } 

    #blank line 
    $line=<inFile>;
    print outFile "$val[$iFile]\n";
    #ne,np
    $line=<inFile>;
    print outFile $line;
    $line =~ s/^\s+|\s+$//g;
    ($ne,$np1) = split(/ +/, $line);

    for (my $i=0; $i < $np1; $i++) {
        $line=<inFile>;
        $line =~ s/^\s+|\s+$//g;
        my ($id,$x,$y,$z) = split(/ +/, $line);
        print outFile "$id $x $y $val[$iFile]\n";
    }
    for (my $i=0; $i < $ne; $i++) {
        $line=<inFile>;
        print outFile $line;
    }

    $iFile=$iFile+1;
    close inFile; close outFile;
}

#prop
#if ( ! open outFile, "> tvd.prop" ){
#    die "can not open tvd.prop";
#} 
#for (my $i=1; $i <= $ne; $i++) {
#    print outFile "$i 1\n";
#}
