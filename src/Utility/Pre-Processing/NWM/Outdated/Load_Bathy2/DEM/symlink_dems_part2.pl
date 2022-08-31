#! /usr/bin/perl -w

# Create symlinks for DEMs (in order of loading)
@order=("ge_ll","sc_ll","gulf_1","gulf_3","cb_ll","db_ll","new_england");

system "rm -rf dem_????.asc";
@files=glob("*.asc");

$nfiles=-1;
for($i=0;$i<@order;$i++) {
  foreach $file (@files) {
    if($file =~ "$order[$i]") {
      $nfiles=$nfiles+1;
      $list[$nfiles]="$file";
    }
  } #foreach
} #for $i

for($i=0;$i<@list;$i++) {
  print "Final list: $list[$i]\n";
  $j=sprintf('%04d',$i);
#  print "doing dem_$j\n";
  system "ln -sf $list[$i] dem_$j.asc";
}#for $i

