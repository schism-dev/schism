#!/usr/bin/env python
import string
import os.path as pth


def cull_depends(infile,outfile, targets):
    inp = open(infile,'r')
    inlines = inp.readlines()
    inp.close()
    out = open(outfile,'w')
    assert inlines[0].startswith("# DO NOT DELETE")
    out.write(inlines[0])
    for line in inlines[1:]:
        target,depends = line.split(":")
	#print "depends: %s" % [pth.splitext(x)[0] for x in depends.strip().split(" ")]
        moddeps = [x for x in depends.strip().split(" ") if pth.splitext(x)[0] in targets]
	#print "targets: %s" % targets
        #print "target: %s" % target
	#print moddeps
        if len(moddeps) > 0:
            out.write("%s: %s\n" % (target,string.join(moddeps," ")))
        
if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    outfile = sys.argv[2]
    targets = [pth.splitext(x)[0] for x in sys.argv[3].split()]
    #print "targets: %s" %targets
    print("Culling dependencies")
    cull_depends(infile,outfile,targets)
