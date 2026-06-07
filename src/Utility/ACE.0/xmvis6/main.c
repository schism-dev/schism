/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2004 Oregon Health and Science University
 * All Rights Reserved
 *
 */

#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/time.h>
#include <pwd.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "defines.h"
#include "globals.h"

#ifndef lint
static char RCSid[] = "$Id: main.c,v 1.30 2011/09/08 19:19:21 pturner Exp $";
#endif

char acevis_version[] = "ACE/vis 6.01 (Sep 8, 2011)";

char *version = acevis_version;

int batchflag = 0;
char batchfile[1024];
char parmfile[1024];
int reading_parms = 0;
int setwin = 0;
int geoflag = 0;
char geoprintstr[1024];

double dcor = 0.0;

double mapscale = 1.0;

/*
 * for image display TODO move to globals.h or .h file
 */
int readimage = 0;
char image_filename[1024];
double imagex;
double imagey;
int readgridfile = 0;

extern int gifdtrans, gifdinterlace;	/* for the gd gif driver */
extern double belel_tol;

int docoords = 0; /* hack to get gif image coords for stations */

char logfilename[2048];
void usage(char *progn);
void writelogfile(char *msg);
void sighandler(int signo);

main(int argc, char **argv)
{
    int i;
    int dofork = 1;
    int readbound = 0;
    int pread = 0;
    int cur_graph = 0;
    int autoscale[200];		/* check for autoscaling override */
    int noautoscale[200];	/* check for no autoscaling override */
    char ACEvisrc_file[256], *s, *getenv(const char *);
    int gno;
    int do_pathl = 0;
    FILE *fp;

    putenv("TZ=UTC");

    sprintf(logfilename, "/usr/local/ace/log/xmvis6-%d.log", (int) getpid());
/* log this execution */
    writelogfile("Start xmvis6");

/* log current directory */
    if (getcwd(buf, 1024) != NULL) {
	writelogfile(buf);
    }

/* set the signal handler */
    signal(SIGSEGV, sighandler);
    signal(SIGINT, sighandler);
    signal(SIGTERM, sighandler);
    strcpy(statusstr, "Program start.");

    statusstr[0] = 0;
    for (i = 0; i < argc; i++) {
	strcat(statusstr, " ");
	strcat(statusstr, argv[i]);
    }
    writelogfile(statusstr);

    strcpy(save_images_fname, "xmv");

    hdevice = GR_PS_L;

    for (i = 0; i < maxgraph; i++) {
	autoscale[i] = 0;
	noautoscale[i] = 0;
    }

    /*
     * initialize the X toolkit unless in batch mode
     */
    for (i = 0; i < argc; i++) {
	if (argmatch(argv[i], "-batch", 4)) {
	    batchflag = TRUE;
	    inwin = 0;
	}
    }
    if (batchflag) {
	if (argc == 1) {
	    fprintf(stderr, "%s: No command line arguments, nothing to do!\n", argv[0]);
	    exit(1);
	}
    } else {
	fprintf(stderr, "\n%s Visualization of Flow and Transport\nCopyright 1990-2004 Oregon Health and Science University\nAll Rights Reserved\n\n", acevis_version);
	initialize_screen(&argc, argv);
    }


    /* initialize colormap data */
    initialize_cms_data();

    /* initialize plots, strings, graphs */
    set_program_defaults();
    /*
     * initialize drogue colors
     */
    init_drog_colors();

    /* initialize the rng */
    srand48(100L);

    /* initialize device, here tdevice is always 0 = Xlib */
    device = tdevice = 0;

    /* check for startup file in effective users home dir */
    if ((s = getenv("HOME")) != NULL) {
	strcpy(ACEvisrc_file, s);
	strcat(ACEvisrc_file, "/.ACEvisrc");
	if ((fp = fopen(ACEvisrc_file, "r")) != NULL) {
	    fclose(fp);
	    getparms(cur_graph, ACEvisrc_file);
	}
    }
    /*
     * check for changed printer spooling strings
     */
    if ((s = getenv("GR_PS_PRSTR")) != NULL) {
	strcpy(ps_prstr, s);
    }
    /*
     * check for changed hardcopy device
     */
    if ((s = getenv("GR_HDEV")) != NULL) {
	hdevice = atoi(s);
    }
    set_printer(hdevice, NULL);
    set_graph_active(cur_graph);

    if (argc > 0) {
	for (i = 1; i < argc; i++) {
	    if (argv[i][0] == '-' || argv[i][0] == '+') {

		if (argmatch(argv[i], "-debug", 6)) {
		    i++;
		    debuglevel = atoi(argv[i]);
		} else if (argmatch(argv[i], "-geo", 4)) {
		    geoflag = 1;
		} else if (argmatch(argv[i], "-batch", 3)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for batch file\n");
			usage(argv[0]);
			exit(1);
		    } else {
			strcpy(batchfile, argv[i]);
			if (batchfile[0]) {
			    batchflag = 1;
			}
			/* runbatch(batchfile); /* testing */
		    }
		} else if (argmatch(argv[i], "-window", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing arguments for setting window dimensions\n");
			usage(argv[0]);
			exit(1);
		    } else {
			winsetwidth = atoi(argv[i]);
			i++;
			winsetheight = atoi(argv[i]);
			setwin = 1;
		    }
		} else if (argmatch(argv[i], "-version", 2)) {
		    extern char compiled_version[];	/* defined in vers.c */

		    printf("%s\n\n", compiled_version);
		    exit(1);
		} else if (argmatch(argv[i], "-swap", 5)) {
		    swapBytes = 1;
		} else if (argmatch(argv[i], "-noswap", 7)) {
		    swapBytes = 0;
		} else if (argmatch(argv[i], "-bs", 3)) {
		    xlibsetbackingstore(1);
		    do_backingstore = 1;
		} else if (argmatch(argv[i], "-belel", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for belel tolerance (for finding elements).\n");
			usage(argv[0]);
			exit(1);
		    } else {
			belel_tol = -atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-mapscale", 9)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for mapscale.\n");
			usage(argv[0]);
			exit(1);
		    } else {
			mapscale = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-drogues", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for Drogues file\n");
			usage(argv[0]);
		    } else {
			if (readdrogues(0, argv[i], -1, 0, 0)) {
			    set_clock(0, drogues[0].start, drogues[0].stop, drogues[0].step, drogues[0].nsteps);
			    load_clock(DROGUES, 0);
			}
		    }
		} else if (argmatch(argv[i], "-pathlines", 10)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for Drogues file\n");
			usage(argv[0]);
		    } else {
			do_pathl = atoi(argv[i]);
			i++;
			if (readdrogues(do_pathl, argv[i], -1, 0, 0)) {
			    set_clock(0, drogues[0].start, drogues[0].stop, drogues[0].step, drogues[0].nsteps);
			    load_clock(DROGUES, 0);
			}
		    }

		} else if (argmatch(argv[i], "-flow", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for SELFE .64 file\n");
			usage(argv[0]);
		    } else {
                        if (ReadElcircSurfNew(0, argv[i], 0, 0, 95, 1, 0, 0.0, 0)) {
                            g[cg].flowt[0].display = NODES;
                            set_clock(0, flowt[0].start, flowt[0].stop, flowt[0].step, flowt[0].nsteps);
                            load_clock(ADCIRC, 0);
			}
		    }
		} else if (argmatch(argv[i], "-surf", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for SELFE .63 file\n");
			usage(argv[0]);
		    } else {
                        if (ReadElcircSurfNew(0, argv[i], 0, 0, 95, 1, 0, 0.0, 0)) {
                            g[cg].flowt[0].display_elev = ON;
                            set_clock(0, flowt[0].start, flowt[0].stop, flowt[0].step, flowt[0].nsteps);
                            load_clock(ADCIRC, 0);
                        }
		    }
		} else if (argmatch(argv[i], "-bot", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for SELFE .63 file\n");
			usage(argv[0]);
		    } else {
                        if (ReadElcircSurfNew(0, argv[i], 1, 0, 95, 1, 0, 0.0, 0)) {
                            g[cg].flowt[0].display_elev = ON;
                            set_clock(0, flowt[0].start, flowt[0].stop, flowt[0].step, flowt[0].nsteps);
                            load_clock(ADCIRC, 0);
                        }
		    }
		} else if (argmatch(argv[i], "-tea", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for TEA-NL file\n");
			usage(argv[0]);
		    } else {
			readflow(0, 0, argv[i]);
		    }
		} else if (argmatch(argv[i], "-teanl", 6)) {
		    int do_teanl = 0;

		    i++;
		    do_teanl = atoi(argv[i]);
		    if (do_teanl < 0) {
			fprintf(stderr, "TEA-NL flow number < 0, should be from 1 - %d\n", MAXTEANL);
			usage(argv[0]);
		    }
		    if (i == argc) {
			fprintf(stderr, "Missing argument for TEA-NL file\n");
			usage(argv[0]);
		    } else {
			i++;
			readflow(do_teanl, 0, argv[i]);
		    }
		} else if (argmatch(argv[i], "-tides", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for Tidal stations file\n");
			usage(argv[0]);
		    } else {
			ReadTideStations(argv[i]);
		    }
		} else if (argmatch(argv[i], "-time", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for time\n");
			usage(argv[0]);
			exit(1);
		    } else {
			timeclock.start = atof(argv[i]);
			i++;
			timeclock.step = atof(argv[i]);
			i++;
			timeclock.nsteps = atoi(argv[i]);
			set_clock(1, timeclock.start, timeclock.step * timeclock.nsteps, timeclock.step, timeclock.nsteps);
		    }
		} else if (argmatch(argv[i], "-coords", 7)) {
		    docoords = 1;
		} else if (argmatch(argv[i], "-gif", 4)) {
		    set_printer(15, NULL);
		} else if (argmatch(argv[i], "-gift", 5)) {
		    gifdtrans = 1;
		} else if (argmatch(argv[i], "-gifi", 5)) {
		    gifdinterlace = 1;
		} else if (argmatch(argv[i], "-gifp", 5)) {
		    int xmin;
		    int xmax;
		    int ymin;
		    int ymax;
		    double charsize;
		    int symsize;
		    int xticl;
		    int yticl;
		    int arrowlength;
		    i++;
		    strcpy(buf, argv[i]);
		    sscanf(buf, "%d,%d,%d,%d,%lf,%d,%d,%d,%d", &xmin, &xmax, &ymin, &ymax, &charsize, &symsize, &xticl, &yticl, &arrowlength);
		    gdsetdevice(xmin, xmax, ymin, ymax, charsize, symsize, xticl, yticl, arrowlength);
		} else if (argmatch(argv[i], "-printfile", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing file name for printing\n");
			usage(argv[0]);
		    } else {
			ptofile = TRUE;
			strcpy(printstr, argv[i]);
		    }
		} else if (argmatch(argv[i], "-vhist", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for time history file\n");
			usage(argv[0]);
		    } else {
			readflowh2d(0, 0, argv[i]);
			set_clock(0, flowh[0].start, flowh[0].stop, flowh[0].step, flowh[0].nsteps);
			load_clock(HISTORY, FLOW, 0);
			g[cg].flowh[0].display = ON;
		    }
		} else if (argmatch(argv[i], "-ela", 4)) {
		    int do_ela = 0;
		    i++;
		    do_ela = atoi(argv[i]);
		    if (do_ela < 0) {
			fprintf(stderr, "ELA concentration number < 0, should be from 1 - %d\n", MAXELA);
			usage(argv[0]);
		    }
		    if (i == argc) {
			fprintf(stderr, "Missing argument for ELA file\n");
			usage(argv[0]);
		    } else {
			i++;
			if (readelaconcbin(do_ela, 0, argv[i])) {
			    set_clock(0, elaconc[do_ela].start, elaconc[do_ela].stop, elaconc[do_ela].step, elaconc[do_ela].nsteps);
			    load_clock(ELA, 0);
			    autoscale_isolines(elaconc[do_ela].cmin, elaconc[do_ela].cmax, &g[cg].elaconc[do_ela].ip);
			} else {
			    fprintf(stderr, "Unable to read ELA concentration file %s\n", argv[i]);
			}
		    }
		} else if (argmatch(argv[i], "-bound", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for boundary file\n");
			usage(argv[0]);
			exit(1);
		    } else {
			if (readboundary2(0, argv[i])) {
			} else {
			    fprintf(stderr, "Can't read boundary file\n");
			    usage(argv[0]);
			}
		    }
		} else if (argmatch(argv[i], "-stations", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for station file\n");
			usage(argv[0]);
			exit(1);
		    } else {
			ReadStations(argv[i], &nsta, &sta);
		    }
		} else if (argmatch(argv[i], "-build", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for build points file\n");
			usage(argv[0]);
			exit(1);
		    } else {
			ReadBuildPoints(argv[i], &nsta, &sta);
		    }
		} else if (argmatch(argv[i], "-pexec", 4)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for command execution\n");
			usage(argv[0]);
		    } else {
			i++;
			execute_command(argv[i]);
		    }
		} else if (argmatch(argv[i], "-rvideo", 2)) {
		    revflag = 1;
		} else if (argmatch(argv[i], "+phase", 6)) {
		    flowf[curteanl].sign_of_phase = -1;
		} else if (argmatch(argv[i], "-phase", 6)) {
		    flowf[curteanl].sign_of_phase = -1;
		} else if (argmatch(argv[i], "-dcor", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for TEANL depth correction\n");
			usage(argv[0]);
		    } else {
			flowf[0].dcor = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-param", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameter file name\n");
			usage(argv[0]);
		    } else {
			strcpy(parmfile, argv[i]);
			if (!getparms(cur_graph, argv[i])) {
			    g[cur_graph].parmsread = FALSE;
			    pread = 0;
			    fprintf(stderr, "Unable to read parameter file %s\n", argv[i]);
			    exit(1);
			} else {
			    pread = 1;
			    g[cur_graph].parmsread = TRUE;
			    set_display_items();
			}
		    }
		} else if (argmatch(argv[i], "-graph", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameter for graph select\n");
			usage(argv[0]);
		    } else {
			sscanf(argv[i], "%d", &gno);
			if (gno >= 0 && gno < maxgraph) {
			    cg = cur_graph = gno;
			    set_graph_active(gno);
			} else {
			    fprintf(stderr, "Graph number must be between 0 and %d\n", maxgraph - 1);
			}
		    }
		} else if (argmatch(argv[i], "-graphtype", 7)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for graph type\n");
			usage(argv[0]);
		    } else {
			if (!strcmp("xy", argv[i])) {
			    g[cur_graph].type = XY;
			} else if (!strcmp("fixed", argv[i])) {
			    g[cur_graph].type = XYFIXED;
			} else if (!strcmp("logx", argv[i])) {
			    g[cur_graph].type = LOGX;
			} else if (!strcmp("logy", argv[i])) {
			    g[cur_graph].type = LOGY;
			} else if (!strcmp("logxy", argv[i])) {
			    g[cur_graph].type = LOGXY;
			} else {
			    fprintf(stderr, "%s: Improper argument for -graphtype should be one of 'xy', 'logx', 'logy', 'logxy'\n", argv[0]);
			    usage(argv[0]);
			}
		    }
		} else if (argmatch(argv[i], "-nofork", 7)) {
		    dofork = 0;
		} else if (argmatch(argv[i], "-grid", 5)) {
		    int do_grid = 0;
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for grid\n");
			usage(argv[0]);
		    } else {
			do_grid = atoi(argv[i]);
			i++;
			readgrid(do_grid, argv[i]);
			compute_boundary(do_grid);
		    }
		} else if (argmatch(argv[i], "-image", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for image file\n");
			usage(argv[0]);
		    } else {
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
		} else {
		    fprintf(stderr, "Unknown command line option %s\n", argv[i]);
		    usage(argv[0]);
		}
	    } else {
		if (!readgridfile && !readgrid(0, argv[i])) {
		    fprintf(stderr, "Warning, no grid file read from %s\n", argv[i]);
		} else {
		    if (readgridfile) {
			fprintf(stderr, "grid file %s already read\n", argv[i]);
		    }
		    readgridfile = 1;
		}
	    }
	}
    }
    if (readgridfile) {
	if (!g[cg].parmsread) {
	    autoscale_grid(cg, 0);
	    set_defaults(cg);
	}
    }

    cg = 0;			/* set current graph */

    dobackground = 1;
    hselectfont(4);
    if (!readbound) {
	if (object_isactive(GRID, 0)) {
	    compute_boundary(0);
	} else {
	    //fprintf(stderr, "Grid 0, is not active, boundary not computed, plot may need to be autoscaled\n");
	}
    }
#ifndef WIN32
/*
    if (dofork && (fork() != 0)) {
	exit(1);
    }
*/
#endif
    if (!batchflag) {
	do_main_loop(argc, argv);
    } else {			/* batch job */
	runbatch(batchfile);
    }
}

void usage(char *progn)
{
    fprintf(stderr, "Usage: %s -batch [batch file] -par [parameter file] [grid file]\n", progn);
    fprintf(stderr, "     -batch [batch file] \n");
    fprintf(stderr, "     -par [parameter file]\n");
    fprintf(stderr, "[grid file] is mandatory, other files and flags are optional.\n");
    fprintf(stderr, "For comments, suggestions, and bug reports contact ace_bugs@ccalmr.ogi.edu (subject ACE/vis)\n\n");
    exit(1);
}

void writelogfile(char *msg)
{
    int i;
    time_t t;
    uid_t uid;
    struct passwd *p;
    char *sp = NULL;
    FILE *fp;
    if (!debuglevel) {
	return;
    }
    fp = fopen(logfilename, "ab");
    if (fp == NULL) {
	return;
    }
    t = time(NULL);
    uid = getuid();
    p = getpwuid(uid);
    sp = ctime(&t);
    if (sp != NULL && strlen(sp) > 0) {
	for (i = 0; i < strlen(sp); i++) {
	    if (sp[i] == '\r' || sp[i] == '\n') {
		sp[i] = 0;
		break;
	    }
	}
	if (msg != NULL) {
	    fprintf(fp, "xmvis6 %s %s \"%s\"\n", p->pw_name, sp, msg);
	} else {
	    fprintf(fp, "xmvis6 %s %s\n", p->pw_name, sp);
	}
    }
    fclose(fp);
}

void sighandler(int signo)
{
    time_t t;
    uid_t uid;
    struct passwd *p;
    FILE *fp;
    if (!debuglevel) {
	exit(1);
    }
    fp = fopen(logfilename, "ab");
    if (fp == NULL) {
	exit(1);
    }
    t = time(NULL);
    uid = getuid();
    p = getpwuid(uid);
    switch (signo) {
    case SIGSEGV:
	fprintf(fp, "xmvis6 SEGMENTATION fault: %s at %s", p->pw_name, ctime(&t));
	fprintf(fp, "\nCrash string: %s\n", statusstr);
	fprintf(stderr, "\nSorry, xmvis6 SEGMENTATION fault, send report to pturner@ccalmr.ogi.edu,\ninclude as much information as possible about the circumstances that lead to this crash.\n");
	fprintf(stderr, "\nStatus string: %s\n", statusstr);
	fclose(fp);
	exit(1);
	break;
    case SIGTERM:
	fprintf(fp, "xmvis6 SEGTERM: %s at %s", p->pw_name, ctime(&t));
	fprintf(fp, "\nStatus string: %s\n", statusstr);
	fclose(fp);
	exit(1);
	break;
    case SIGINT:
	fprintf(fp, "xmvis6 SEGINT: %s at %s", p->pw_name, ctime(&t));
	fprintf(fp, "\nStatus string: %s\n", statusstr);
	fclose(fp);
	exit(1);
	break;
    }
}
