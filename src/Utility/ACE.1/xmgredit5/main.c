/*
 * ACE/gredit - 2d finite element grid generation
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *                      All Rights Reserved.
 *
 */

/*
 * main.c - entry point for gredit
 *
 */

#ifndef lint
static char RCSid[] = "$Id: main.c,v 1.11 2011/09/14 17:45:01 pturner Exp $";
char xmgredit_version[] = "ACE/gredit (xmgredit) version 5.1 (Wed Sep 14 10:20:23 PDT 2011)";

#endif

#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "globals.h"

static void usage(char *progn);

extern char compiled_version[];
char version[128];

int glflag = 1;
int dosolidbox = 1;

double mapscale = 1.0;

extern int save_isoline;
extern int draw_circum;

/*
 * for image display TODO move to globals.h or .h file
 */
int readimage = 0;
char image_filename[1024];
char bcfile[2048];
double imagex;
double imagey;

main(int argc, char **argv)
{
    extern int revflag;
    int i = 0, c, readback = 0, did_bath = 0, wflag = 0, lflag = 0,
     tflag = 0, cflag = 0, cbc = 0;
    int dofork = 1, doedge = 0, dopom = 0;
    double x1, y1, x2, y2;
    char *s, *getenv(const char *);
    char coastfile[512];
    char gridfile[512];
    char backfile[512];
    char depthfile[512];
    depthfile[0] = 0;
    
    printf("\n%s\nSemi-automatic generation of 2D finite element grids\nCopyright 1990-2004 Oregon Health and Science University\nAll Rights Reserved\n\n", xmgredit_version);
    printf("\nSweep2 - by Steven Fortune.  Copyright (c) 1994 by AT&T Bell Laboratories\n\n");
#ifdef DO_TRIANGLE
    printf("\nTriangle - A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator\nCopyright 1996, Jonathan Richard Shewchuk\nSchool of Computer Science\nCarnegie Mellon University\n\n");
#endif
    printf("Loading...");
    fflush(stdout);

    strcpy(version, xmgredit_version);
    if (argc == 1) {
	usage(argv[0]);
    }
    set_program_defaults();
    initialize_screen(&argc, argv);
    if (argc >= 2) {
	for (i = 1; i < argc; i++) {
	    if (argv[i][0] == '-' || argv[i][0] == '+') {
/* Turn on debugging */
		if (argmatch(argv[i], "-debug", 6)) {
		    i++;
		    debuglevel = atoi(argv[i]);
		} else if (argmatch(argv[i], "-mindist", 2)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for minimum distance between build points\n");
			usage(argv[0]);
			exit(0);
		    } else {
/* When triangulating, eliminate points too close to each other */
			i++;
			mindist = atof(argv[i]);
			printf("Using minimum distance for triangulating of %lf\n", mindist);
		    }
                } else if (argmatch(argv[i], "-bc", 3)) {
/* Compute boundary and read bc file */
                    if (i == argc) {
                        fprintf(stderr, "Missing argument for BC file\n");
                        usage(argv[0]);
                        exit(1);
                    } else {
                        i++;
                        cbc = 1;
                        strcpy(bcfile, argv[i]);
                    }
		} else if (argmatch(argv[i], "-swap", 5)) {
/* For binary files create on systems with different byte order */
		    swapBytes = 1;
		} else if (argmatch(argv[i], "-noswap", 7)) {
/* For binary files create on systems with different byte order */
		    swapBytes = 0;
		} else if (argmatch(argv[i], "-nosolidbox", 11)) {
/* Hack to shut down filled boxes */
		    dosolidbox = 0;
		} else if (argmatch(argv[i], "-load", 2)) {
/* Load grid to background grid */
		    lflag = 1;
		} else if (argmatch(argv[i], "-nofork", 7)) {
/* Don't fork a copy of gredit, useful for debugging */
		    dofork = 0;
		} else if (argmatch(argv[i], "-background", 12)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for background grid file\n");
			usage(argv[0]);
			exit(0);
		    } else {
/* read a different grid for the background grid */
			i++;
			readback = 1;
			strcpy(backfile, argv[i]);
		    }
		} else if (argmatch(argv[i], "-nofixedscale", 13)) {
/* Allow for varying scales in X and Y */
		    fixed_scale = 0;
		} else if (argmatch(argv[i], "-bs", 3)) {
/* Set backing store */
		    xlibsetbackingstore(1);
		} else if (argmatch(argv[i], "-edgelist", 3)) {
		    doedge = 1;
                } else if (argmatch(argv[i], "-stations", 4)) {
                    i++;
                    if (i == argc) {
                        fprintf(stderr, "Missing argument for station file\n");
                        usage(argv[0]);
                        exit(1);
                    } else {
                        ReadStations(argv[i], &nsta, &sta);
                    }
		} else if (argmatch(argv[i], "-bathymetry", 2)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for triangulating buildpoints.\n");
			usage(argv[0]);
			exit(0);
		    } else {
/* triangulate build points */
			i++;
			if (readbuild(curbuild, argv[i])) {
			    elim_build_dupes(curbuild, mindist);
			    triang(curbuild);
			    load_grid(backgrid, curbuild);
			    did_bath = 1;
			}
		    }
		} else if (argmatch(argv[i], "-coastline", 2)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for coastal outline file.\n");
			usage(argv[0]);
			exit(0);
		    } else {
/* read a coastal outline file */
			cflag = 1;
			i++;
			strcpy(coastfile, argv[i]);
		    }
		} else if (argmatch(argv[i], "-dpe", 4)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for depths file.\n");
			usage(argv[0]);
			exit(0);
		    } else {
			i++;
			strcpy(depthfile, argv[i]);
		    }
		} else if (argmatch(argv[i], "-override", 2)) { 
		    if (i == argc) {
			fprintf(stderr, "Missing argument build points override\n");
			usage(argv[0]);
			exit(0);
		    } else {
			i++;
/* read a subset of build points */
/* The value of noverride is the number of points to read */
			noverride = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-q", 2)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for quad cutoff.\n");
			usage(argv[0]);
			exit(0);
		    } else {
			i++;
/* experiment with a parameter used to create quadrangles */
			qcutoff = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-saveisoline", 6)) {
/* Save a particular isoline value */
		    if (i == argc) {
			fprintf(stderr, "Missing argument for save isoline.\n");
			usage(argv[0]);
			exit(0);
		    } else {
			i++;
			save_isoline = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-mapscale", 9)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for mapscale.\n");
			usage(argv[0]);
			exit(0);
		    } else {
			i++;
/* experiment with the mapscale */
			mapscale = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-f", 2)) {
/* File is binary */
		    tflag = 1;
		} else if (argmatch(argv[i], "-binary", 7)) {
/* File is binary */
		    tflag = 1;
		} else if (argmatch(argv[i], "-w", 2)) {
/* Make a 2 element grid based on build points */
		    wflag = 1;
		} else if (argmatch(argv[i], "-version", 2)) {
/* Print the version */
		    printf("%s\n\n", compiled_version);
		    exit(0);
		} else if (argmatch(argv[i], "-generate", 2)) {
/* Generate a grid */
		    wflag = 1;
		} else if (argmatch(argv[i], "-World", 2)) {
		    if (i == argc) {
			fprintf(stderr, "Missing arguments for world scaling\n");
			usage(argv[0]);
			exit(0);
		    } else {
/* Generate a grid based on the corners given on the command line */
			i++;
			wflag = 2;
			x1 = atof(argv[i++]);
			y1 = atof(argv[i++]);
			x2 = atof(argv[i++]);
			y2 = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-d", 2)) {
/* For testing and debugging */
		    wflag = 2;
		    x1 = 0.0;
		    y1 = 0.0;
		    x2 = 10.0;
		    y2 = 10.0;
		    dofork = 0;
		} else if (argmatch(argv[i], "-rvideo", 2)) {
/* reverse video */
		    revflag = 1;
		} else if (argmatch(argv[i], "-nlist", 3)) {
		    extern int displaynl;
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing arguments for node list file\n");
			usage(argv[0]);
			exit(0);
		    } else {
/* display node list */
			if (!ReadPNodeList(argv[i])) {
			    displaynl = 1;
			}
		    }
		} else if (argmatch(argv[i], "-belel", 6)) {
		    if (i == argc) {
			fprintf(stderr, "Missing arguments for belel\n");
			usage(argv[0]);
			exit(0);
		    } else {
			i++;
/* Set the tolerance for point location in an element */
			belel_tol = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-image", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for image file\n");
			usage(argv[0]);
		    } else {
/* read an image */
			readimage = TRUE;
			strcpy(image_filename, argv[i]);
		    }
		} else if (argmatch(argv[i], "-imagexy", 8)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for imagexy\n");
			usage(argv[0]);
		    }
		    imagex = atof(argv[i]);
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for imagexy\n");
			usage(argv[0]);
		    }
		    imagey = atof(argv[i]);
		} else if (argmatch(argv[i], "-pom", 4)) {
/* Activate the pom interface */
		    dopom = 1;
		} else if (argmatch(argv[i], "-circum", 7)) {
/* Draw circumcenters, testing only */
		    draw_circum = 1;
		} else {
		    fprintf(stderr, "Command line argument not found\n");
		    usage(argv[0]);
		    exit(0);
		}

	    } else {
		/* a data file (we hope) */
		if (!(i >= argc)) {
		    strcpy(gridfile, argv[i]);
		}		/* end if */
	    }			/* end if '-' */
	}			/* end for */
    }				/* end if argc >= 2) */
    if (!wflag && !lflag) {
	if (tflag) {
	    if (!readgridbin(curgrid, gridfile)) {
		fprintf(stderr, "Error reading file %s\n", gridfile);
		exit(1);
	    }
	} else {
	    if (dopom) {
		if (!readpomgrid(curgrid, gridfile)) {
		    fprintf(stderr, "Error reading file %s\n", gridfile);
		    exit(1);
		}
	    } else {
		if (!readgrid(curgrid, gridfile)) {
		    fprintf(stderr, "Error reading file %s\n", gridfile);
		    exit(1);
		}
	    }
	}
    } else if (wflag == 1) {
	if (!gengrid(curgrid, curbuild)) {
	    fprintf(stderr, "Unable to create grid\n");
	    exit(1);
	}
    } else if (wflag == 2) {
	if (!gengrid2(curgrid, x1, y1, x2, y2)) {
	    fprintf(stderr, "Unable to create grid\n");
	    exit(1);
	}
    } else if (lflag && did_bath) {
	copy_grid(backgrid, curgrid);
    }
    if (readback) {
	if (!readgridbin(backgrid, backfile)) {
	    fprintf(stderr, "Error reading file %s\n", backfile);
	    exit(1);
	}
    }
    if (!readback && !did_bath) {
	copy_grid(curgrid, MAXGRIDS);
    }
    if (cflag) {
	readboundary(MAXGRIDS, coastfile, 0);
    }
/*
   getassoc_elements(0);
 */
    if (cbc) {
        FILE *fp;
        compute_boundary(0);
        file_bc(bcfile);
        fp = fopen("bcs.bnd", "w");
        if (fp != NULL) {
            write_adcirc_openb(0, fp);
            fclose(fp);
        }
    }

    hselectfont(4);

    if (dofork && (fork() != 0)) {
        exit(0);
    }
    printf("...\n");
    if (doedge) {
	getedgelist2(0);
    }
    if (depthfile[0]) {
	readdepths(0, depthfile);
    }
    do_main_loop();
/* NOTREACHED */
}

static void usage(char *progn)
{
    fprintf(stderr, "%s -m [min_dist] -b [bathymetry_file] -w -f      edit_grid_filename\n", progn);
    fprintf(stderr, "-m [min_dist]   ...  Minumum distance to allow between points when triangulating\n");
    fprintf(stderr, "-b [bathymetry] ...  File of bathymetry to triangulate\n");
    fprintf(stderr, "-w              ...  generate an editable grid from the limits\n");
    fprintf(stderr, "-generate       ...  generate an editable grid from the limits\n");
    fprintf(stderr, "                     found in the file of bathymetry\n");
    fprintf(stderr, "-W xmin ymin xmax ymax ...  Generate a grid based on the corners given on the command line\n");
    fprintf(stderr, "-f              ...  indicates the grid in 'edit_grid_filename' is\n");
    fprintf(stderr, "                     in binary format\n");
    fprintf(stderr, "-binary         ...  same as -f\n");
    fprintf(stderr, "-swap           ...  For binary files created on systems with different byte order\n");
    fprintf(stderr, "-noswap         ...  For binary files created on systems with different byte order\n");

    fprintf(stderr, "-belel [belel_tol] ... Set the tolerance for point location in an element\n");
    fprintf(stderr, "-pom            ... interpret edit_grid as a POM file\n");
    fprintf(stderr, "-version        ... Print the version\n");
    fprintf(stderr, "-l              ... Load grid to background grid\n");
    fprintf(stderr, "-background [file] ... Read 'file' into the background grid\n");
    fprintf(stderr, "-rvideo         ... reverse video\n");
    fprintf(stderr, "-nosolidbox     ... Hack to shut down filled boxes\n");
    fprintf(stderr, "-nofork         ... Don't fork a copy of gredit, useful for debugging\n");
    fprintf(stderr, "-nofixedscale   ... Allow for varying scales in X and Y \n");
    fprintf(stderr, "-bs             ... Set X11 backing store \n");
    fprintf(stderr, "-edgelist       ... write edges to 'grid.tmp' \n");
    fprintf(stderr, "-stations [stationfile] ... Read station file: stations\\nlabel\\nx y\\n... \n");
    fprintf(stderr, "-coastline [coastfile]  ... read a coastal outline file\n");
    fprintf(stderr, "-dpe [depthfile]   ... read a depthfile \n");
    fprintf(stderr, "-o [numpoints]     ... Override the build points and read only 'numpoints' \n");
    fprintf(stderr, "-q [quadcutoff]    ... experiment with a parameter used to create quadrangles\n");
    fprintf(stderr, "-saveisoline [i]   ... Save a particular isoline value\n");
    fprintf(stderr, "-mapscale [mapscale] ... experiment with the mapscale\n");
    fprintf(stderr, "-d                 ... set some debugging parameters\n");
    fprintf(stderr, "-nlist [nodePropFile] ... display node list\n");
    fprintf(stderr, "-image [image.xwd] ... read an xwd format image\n");
    fprintf(stderr, "-imagexy [x] [y]   ... locate the image\n");
    fprintf(stderr, "-circum            ... Draw circumcenters, testing only\n");

    fprintf(stderr, "\nSee http://www.stccmop.org/ace/gredit/ online documentation.\n");
    fprintf(stderr, "For comments, suggestions, and bug reports contact gredit@ccalmr.ogi.edu\n\n");
    exit(0);
}
