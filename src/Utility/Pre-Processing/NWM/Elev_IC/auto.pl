#! /usr/bin/perl -w

#load some functions
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd;


print("You may need to recompile 'gen_elev.f90'\n");
print("e.g: ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_elev gen_elev.f90\n\n\n");

system("ln -sf ../hgrid.ll hgrid.gr3");

system("./gen_elev");

unlink("../elev.ic");
copy("elev.ic","../elev.ic");

print("Done.\n")
