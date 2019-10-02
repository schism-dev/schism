#! /usr/bin/perl -w

# Add a '0' for sflux*.nc for new versions.
# Run inside sflux

opendir(DIR, ".") or die "can't opendir: $!";
while (defined($file = readdir(DIR))) {
  if($file =~ "sflux_" && $file =~ ".nc") {
    ($a,$b,$c)=split("[.]",$file);
    #print "$a.0$b.$c\n";
    if(length($b)==3) {
      print "doing $file\n";
      system "ln -sf $file $a.0$b.$c";
    }
  }
}
closedir(DIR);
