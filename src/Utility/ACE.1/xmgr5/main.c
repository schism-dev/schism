/* $Id: main.c,v 1.7 2008/10/08 22:53:19 pturner Exp $

   (C) COPYRIGHT 1991-1995 Paul J Turner
   All Rights Reserved

   Copyright 1995 by Paul J Turner, Portland, OR

   All Rights Reserved

   Permission to use, copy, modify, and distribute this software and its
   documentation for any purpose and without fee is hereby granted,
   provided that the above copyright notice appear in all copies and that
   both that copyright notice and this permission notice appear in
   supporting documentation.

   PAUL J TURNER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
   ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
   PAUL J TURNER BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
   ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
   WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
   ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
   SOFTWARE.

   * ACE/gr - Graphics for Exploratory Data Analysis
   *
   * Paul J. Turner
   * pturner@teleport.com
   *
 */

/* for globals.h */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>

#define MAIN

#include "globals.h"
#include "patchlevel.h"

char version[] = "ACE/gr v5";

static void usage(char *progname);

extern int epsflag;		/* defined in ps.c - force eps to be written */

#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)
int readcdf = 0;
char netcdf_name[512], xvar_name[128], yvar_name[128];
#endif

extern char dbuser[], dbpasswd[], dbase[], dbhost[];

void main(int argc, char *argv[])
{
    FILE *fp;
    int i, j;
    int gno;
    int cur_graph = cg;		/* default graph is graph 0 */
    int loadlegend = 0;		/* legend on and load file names */
    int gcols = 1, grows = 1;
    int autoscale[200];		/* check for autoscaling override */
    int noautoscale[200];	/* check for no autoscaling override */
    int grbatch = 0;		/* if executed as 'grbatch' then TRUE */
    int remove = 0;		/* remove file after read */
    int noprint = 0;		/* if grbatch, then don't print if true */
    extern int gotbatch;	/* flag if a parameter file to execute */
    extern int realtime;	/* if data is to be plotted as it is read in */
    extern int inpipe;		/* if xmgr is to participate in a pipe */
    extern char batchfile[];	/* name of file to execute */
    extern int density_flag;	/* temp, for density plots, defined in
				 * plotone.c */
    extern int render_tri;	/* Render every three as a tri */
    extern int ptofile;

    extern int gifdtrans, gifdinterlace; /* for the gd gif driver */

    char xvgrrc_file[256], *s;

    getcwd(currentdir, 1024);	/* get the starting directory */

#ifdef PGSQL
    SetDBDefaults();
#endif

#ifndef GR_HELPFILE
    strcpy(help_file, "ACEgr.html");
#else
    strcpy(help_file, GR_HELPFILE);
#endif

#ifndef GR_HELPVIEWER
    strcpy(help_viewer, "netscape");
#else
    strcpy(help_viewer, GR_HELPVIEWER);
#endif

    for (i = 0; i < maxgraph; i++) {
	autoscale[i] = 0;
	noautoscale[i] = 0;
    }

    for (i = 0; i < argc; i++) {
        if (argmatch(argv[i], "-batch", 4)) {
	    grbatch = TRUE;
        }
    }
    if (grbatch) {
    } else {
/*
	fprintf(stderr, "(C) Copyright 1991-2002 Paul J Turner\n");
	fprintf(stderr, "All Rights Reserved\n");*/
	initialize_screen(&argc, argv);
    }

    /* initialize plots, strings, graphs */
    set_program_defaults();

    /* initialize monitors */
    initMonitor();

    /* initialize colormap data */
    initialize_cms_data();

    /* initialize the rng */
    srand48(100L);

    /* initialize device, here tdevice is always 0 = Xlib */
    device = tdevice = 0;

    /* check for startup file in local directory */
    if ((fp = fopen(".xmgrrc", "r")) != NULL) {
	fclose(fp);
	getparms(cur_graph, ".xmgrrc");
    } else {
	/* check for startup file in effective users home dir */
	if ((s = getenv("HOME")) != NULL) {
	    strcpy(xvgrrc_file, s);
	    strcat(xvgrrc_file, "/.xmgrrc");
	    if ((fp = fopen(xvgrrc_file, "r")) != NULL) {
		fclose(fp);
		getparms(cur_graph, xvgrrc_file);
	    }
	} else {
	    if (!grbatch) {
	    }
	}
    }
    plfile[0] = 0;		/* parameter file name */

    /*
     * ACE/gr home directory
     */
    if ((s = getenv("GR_HOME")) != NULL) {
	strcpy(acegrdir, s);
    } else {
	if (!grbatch) {
	}
    }

    /*
     * check for changed printer spooling strings
     */
    if ((s = getenv("GR_PS_PRSTR")) != NULL) {
	strcpy(ps_prstr, s);
    }
    if ((s = getenv("GR_MIF_PRSTR")) != NULL) {
	strcpy(mif_prstr, s);
    }
    /*
     * check for changed hardcopy device
     */
    if ((s = getenv("GR_HDEV")) != NULL) {
	hdevice = atoi(s);
    }
    set_printer(hdevice, NULL);
    /*
     * check for changed help file location
     */
    if ((s = getenv("GR_HELPFILE")) != NULL) {
	strcpy(help_file, s);
    }
    /*
     * check for changed help file viewer
     */
    if ((s = getenv("GR_HELPVIEWER")) != NULL) {
	strcpy(help_viewer, s);
    }
    set_graph_active(cur_graph);

    if (argc >= 2) {
	for (i = 1; i < argc; i++) {
	    if (argv[i][0] == '-') {
		if (argmatch(argv[i], "-version", 5)) {
		    fprintf(stderr, "%s PL %d\n", version, PATCHLEVEL);
		    exit(0);
		}
		if (argmatch(argv[i], "-debug", 6)) {

// YYDEBUG broken shit main.c:219: undefined reference to `yydebug'			
//#ifdef LOCAL
//#define YYDEBUG
//#endif

//#ifdef YYDEBUG
//		    extern int yydebug;
//#endif

		    i++;
		    debuglevel = atoi(argv[i]);
		    if (debuglevel == 4) {	/* turn on debugging in
						 * pars.y */

//#ifdef YYDEBUG
//			yydebug = 1;
//#endif

		    }
		} else if (argmatch(argv[i], "-autoscale", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for autoscale flag, should be one of 'x', 'y', 'xy'\n");
			usage(argv[0]);
		    } else {
			if (!strcmp("x", argv[i])) {
			    autoscale[cur_graph] = 1;
			} else if (!strcmp("y", argv[i])) {
			    autoscale[cur_graph] = 2;
			} else if (!strcmp("xy", argv[i])) {
			    autoscale[cur_graph] = 3;
			} else {
			    fprintf(stderr, "%s: Improper argument for -a flag should be one of 'x', 'y', 'xy'\n", argv[0]);
			    usage(argv[0]);
			}
		    }
		} else if (argmatch(argv[i], "-noauto", 7)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for no autoscale flag, should be one of 'x', 'y', 'xy'\n");
			usage(argv[0]);
		    } else {
			if (!strcmp("x", argv[i])) {
			    noautoscale[cur_graph] = 1;
			} else if (!strcmp("y", argv[i])) {
			    noautoscale[cur_graph] = 2;
			} else if (!strcmp("xy", argv[i])) {
			    noautoscale[cur_graph] = 3;
			} else {
			    fprintf(stderr, "%s: Improper argument for -noauto flag should be one of 'x', 'y', 'xy'\n", argv[0]);
			    usage(argv[0]);
			}
		    }
		} else if (argmatch(argv[i], "-batch", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for batch file\n");
			usage(argv[0]);
		    } else {
			gotbatch = TRUE;
			ptofile = TRUE;
			strcpy(batchfile, argv[i]);
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
		    imagex = atoi(argv[i]);
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for imagexy\n");
			usage(argv[0]);
		    }
		    imagey = atoi(argv[i]);
		} else if (argmatch(argv[i], "-digit", 6)) {
		    extern int digitflag;	/* TODO to be removed */
		    digitflag = TRUE;
		} else if (argmatch(argv[i], "-pipe", 5)) {
		    int flags = 0;
		    inpipe = TRUE;
/*
   fcntl(0, F_SETFL, flags | O_NONBLOCK);
 */
		} else if (argmatch(argv[i], "-noprint", 8)) {
		    noprint = TRUE;
		} else if (argmatch(argv[i], "-logwindow", 5)) {
		    logwindow = TRUE;
		} else if (argmatch(argv[i], "-nologwindow", 7)) {
		    logwindow = FALSE;
		} else if (argmatch(argv[i], "-hot", 4)) {
		} else if (argmatch(argv[i], "-usehelp", 8)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for use help\n");
			usage(argv[0]);
		    } else {
			use_help = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-npipe", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for named pipe\n");
			usage(argv[0]);
		    } else {
			strcpy(pipe_name, argv[i]);
			named_pipe = 1;
		    }
#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)
		} else if (argmatch(argv[i], "-netcdf", 7) || argmatch(argv[i], "-hdf", 4)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for netcdf file\n");
			usage(argv[0]);
		    } else {
			strcpy(netcdf_name, argv[i]);
			readcdf = 1;
		    }
		} else if (argmatch(argv[i], "-ADP", 4)) {
		    char varname[128];
		    int step;
		    i++;
		    strcpy(varname, argv[i]);
		    i++;
		    step = atoi(argv[i]);
		    netcdfreadADP(cg, -1, netcdf_name, varname, step);
		} else if (argmatch(argv[i], "-ADP2", 5)) {
		    char varname[128];
		    int ind, start, end;
		    i++;
		    strcpy(varname, argv[i]);
		    i++;
		    start = atoi(argv[i]);
		    i++;
		    end = atoi(argv[i]);
		    for (ind = start; ind <= end; ind++) {
			netcdfreadADP(cg, -1, netcdf_name, varname, ind);
		    }
		} else if (argmatch(argv[i], "-netcdfxy", 9) || argmatch(argv[i], "-hdfxy", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for netcdf X variable name\n");
			usage(argv[0]);
		    } else {
			strcpy(xvar_name, argv[i]);
		    }
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for netcdf Y variable name\n");
			usage(argv[0]);
		    } else {
			strcpy(yvar_name, argv[i]);
		    }
		    if (strcmp(xvar_name, "null")) {
			readnetcdf(cg, -1, netcdf_name, xvar_name, yvar_name, -1, -1, 1);
		    } else {
			readnetcdf(cg, -1, netcdf_name, NULL, yvar_name, -1, -1, 1);
		    }
#endif				/* HAVE_NETCDF */
		} else if (argmatch(argv[i], "-timer", 6)) {
		    extern int timer_delay;	/* TODO move to globals.h */

		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for time delay\n");
			usage(argv[0]);
		    } else {
			timer_delay = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-realtime", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for realtime plotting\n");
			usage(argv[0]);
		    } else {
			realtime = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-monitor", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for monitoring a growing file\n");
			usage(argv[0]);
		    } else {
			createQuickMonitor(argv[i]);
		    }
		} else if (argmatch(argv[i], "-maxsets", 8)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for max number of sets\n");
			usage(argv[0]);
		    } else {
			maxplot = atoi(argv[i]);
			realloc_plots(maxplot);
		    }
		} else if (argmatch(argv[i], "-graphsets", 10)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for max number of sets for he current graph\n");
			usage(argv[0]);
		    } else {
			realloc_graph_plots(cur_graph, atoi(argv[i]));
		    }
		} else if (argmatch(argv[i], "-maxgraph", 8)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for max number of graphs\n");
			usage(argv[0]);
		    } else {
			maxgraph = atoi(argv[i]);
			realloc_graphs();
		    }
		} else if (argmatch(argv[i], "-maxlines", 8)) {
		    int itmp;
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for max number of lines\n");
			usage(argv[0]);
		    } else {
			itmp = atoi(argv[i]);
			realloc_lines(itmp);
		    }
		} else if (argmatch(argv[i], "-maxboxes", 8)) {
		    int itmp;
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for max number of boxes\n");
			usage(argv[0]);
		    } else {
			itmp = atoi(argv[i]);
			realloc_boxes(itmp);
		    }
		} else if (argmatch(argv[i], "-maxstr", 7)) {
		    int itmp;
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for max number of strings\n");
			usage(argv[0]);
		    } else {
			itmp = atoi(argv[i]);
			realloc_strings(itmp);
		    }
/*
   } else if (argmatch(argv[i], "-defaultcmap", 12)) {
   use_defaultcmap = 0;
 */
		} else if (argmatch(argv[i], "-vertext", 8)) {
		    use_xvertext = 1;
		} else if (argmatch(argv[i], "-remove", 7)) {
		    remove = 1;
		} else if (argmatch(argv[i], "-fixed", 5)) {
		    i++;
		    canvasw = atoi(argv[i]);
		    i++;
		    canvash = atoi(argv[i]);
		    page_layout = FIXED;
		} else if (argmatch(argv[i], "-landscape", 10)) {
		    page_layout = LANDSCAPE;
		} else if (argmatch(argv[i], "-portrait", 9)) {
		    page_layout = PORTRAIT;
		} else if (argmatch(argv[i], "-free", 5)) {
		    page_layout = FREE;
		} else if (argmatch(argv[i], "-noask", 5)) {
		    noask = 1;
		} else if (argmatch(argv[i], "-mono", 5)) {
		    monomode = 1;
		} else if (argmatch(argv[i], "-bs", 3)) {
		    backingstore = 1;
		} else if (argmatch(argv[i], "-nobs", 5)) {
		    backingstore = 0;
		} else if (argmatch(argv[i], "-dc", 3)) {
		    allow_dc = 1;
		} else if (argmatch(argv[i], "-nodc", 5)) {
		    allow_dc = 0;
		} else if (argmatch(argv[i], "-redraw", 7)) {
		    allow_refresh = 1;
		} else if (argmatch(argv[i], "-noredraw", 9)) {
		    allow_refresh = 0;
		} else if (argmatch(argv[i], "-maxcolors", 10)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for number of colors\n");
			usage(argv[0]);
		    }
		    maxcolors = atoi(argv[i]);
		    if (maxcolors > 256) {
			fprintf(stderr, "Number of colors exceeds 256, set to 256\n");
			maxcolors = 256;
		    }
		    initialize_cms_data();
		} else if (argmatch(argv[i], "-GXxor", 6)) {
		    invert = 0;
		} else if (argmatch(argv[i], "-GXinvert", 6)) {
		    invert = 1;
		} else if (argmatch(argv[i], "-eps", 4)) {
		    epsflag = 1;
		} else if (argmatch(argv[i], "-gif", 4)) {
		    set_printer(5, NULL);
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
		    sscanf(buf, "%d,%d,%d,%d,%lf,%d,%d,%d,%d",
			 &xmin, &xmax, &ymin, &ymax, &charsize, &symsize,
			   &xticl, &yticl, &arrowlength);
		    gdsetdevice(xmin, xmax, ymin, ymax, charsize, symsize,
				xticl, yticl, arrowlength);
		} else if (argmatch(argv[i], "-device", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for hardcopy device select flag\n");
			usage(argv[0]);
		    }
		    set_printer(atoi(argv[i]), NULL);
		} else if (argmatch(argv[i], "-log", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for log plots flag\n");
			usage(argv[0]);
		    }
		    if (!strcmp("x", argv[i])) {
			g[cur_graph].type = LOGX;
		    } else if (!strcmp("y", argv[i])) {
			g[cur_graph].type = LOGY;
		    } else if (!strcmp("xy", argv[i])) {
			g[cur_graph].type = LOGXY;
		    } else {
			fprintf(stderr, "%s: Improper argument for -l flag should be one of 'x', 'y', 'xy'\n", argv[0]);
			g[cur_graph].type = XY;
		    }
		} else if (argmatch(argv[i], "-printfile", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing file name for printing\n");
			usage(argv[0]);
		    } else {
			set_printer(FILEP, argv[i]);
		    }
		} else if (argmatch(argv[i], "-hardcopy", 2)) {
		    grbatch = TRUE;
		} else if (argmatch(argv[i], "-pexec", 6)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for exec\n");
			usage(argv[0]);
		    } else {
			char pstring[256];
			int ilen;

			i++;
			strcpy(pstring, argv[i]);
			ilen = strlen(pstring);
			pstring[ilen] = '\n';
			pstring[ilen + 1] = 0;
			read_param(pstring);
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
		} else if (argmatch(argv[i], "-block", 6)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameter for block data\n");
			usage(argv[0]);
		    } else {
			if (getdata(cur_graph, argv[i], cursource, BLOCK)) {
			}
		    }
		} else if (argmatch(argv[i], "-bxy", 4)) {
		    char blocksetcols[32];
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameter for block data set creation\n");
			usage(argv[0]);
		    }
		    strcpy(blocksetcols, argv[i]);
		    create_set_fromblock(cg, curtype, blocksetcols);
		} else if (argmatch(argv[i], "-xy", 3)) {
		    curtype = XY;
		} else if (argmatch(argv[i], "-bin", 4)) {
		    curtype = BIN;
		} else if (argmatch(argv[i], "-xydx", 5)) {
		    curtype = XYDX;
		} else if (argmatch(argv[i], "-xydy", 5)) {
		    curtype = XYDY;
		} else if (argmatch(argv[i], "-xydxdx", 7)) {
		    curtype = XYDXDX;
		} else if (argmatch(argv[i], "-xydydy", 7)) {
		    curtype = XYDYDY;
		} else if (argmatch(argv[i], "-xydxdy", 7)) {
		    curtype = XYDXDY;
		} else if (argmatch(argv[i], "-xyuv", 5)) {
		    curtype = XYUV;
		} else if (argmatch(argv[i], "-xyz", 4)) {
		    curtype = XYZ;
		} else if (argmatch(argv[i], "-xyztri", 4)) {
		    curtype = XYZ;
		    render_tri = TRUE;
		} else if (argmatch(argv[i], "-xyd", 4)) {
		    curtype = XYZ;
		    density_flag = TRUE;
		} else if (argmatch(argv[i], "-xybox", 6)) {
		    curtype = XYBOX;
		} else if (argmatch(argv[i], "-rectgrid", 9)) {
		    curtype = RECTGRID;
		} else if (argmatch(argv[i], "-fegrid", 7)) {
		    curtype = FEGRID;
		} else if (argmatch(argv[i], "-boxplot", 8)) {
		    curtype = XYBOXPLOT;
		} else if (argmatch(argv[i], "-xyr", 4)) {
		    curtype = XYRT;
		} else if (argmatch(argv[i], "-ihl", 4)) {
		    curtype = IHL;
		} else if (argmatch(argv[i], "-bp", 3)) {
		    curtype = BP;
		} else if (argmatch(argv[i], "-nxy", 4)) {
		    curtype = NXY;
		} else if (argmatch(argv[i], "-ctd", 4)) {
		    curtype = CTD;
		} else if (argmatch(argv[i], "-hilo", 5)) {
		    curtype = XYHILO;
		} else if (argmatch(argv[i], "-type", 2)) {
		    /* other file types here */
		    i++;
		    if (argmatch(argv[i], "bin", 3)) {
			curtype = BIN;
		    } else if (argmatch(argv[i], "xydx", 4)) {
			curtype = XYDX;
		    } else if (argmatch(argv[i], "xydy", 4)) {
			curtype = XYDY;
		    } else if (argmatch(argv[i], "xydxdx", 6)) {
			curtype = XYDXDX;
		    } else if (argmatch(argv[i], "xydydy", 6)) {
			curtype = XYDYDY;
		    } else if (argmatch(argv[i], "xydxdy", 6)) {
			curtype = XYDXDY;
		    } else if (argmatch(argv[i], "xyz", 3)) {
			curtype = XYZ;
		    } else if (argmatch(argv[i], "xyr", 3)) {
			curtype = XYRT;
		    } else if (argmatch(argv[i], "hilo", 4)) {
			curtype = XYHILO;
		    } else if (argmatch(argv[i], "boxplot", 4)) {
			curtype = XYBOXPLOT;
		    } else {
			fprintf(stderr, "%s: Unknown file type '%s' assuming XY\n", argv[0], argv[i]);
			curtype = XY;
		    }
		} else if (argmatch(argv[i], "-missing", 7)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for missing data value.\n");
			usage(argv[0]);
		    } else {
			i++;
			hasmissing = 1;
			missing = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-nomissing", 7)) {
			hasmissing = 0;
		} else if (argmatch(argv[i], "-graphtype", 7)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for graph type\n");
			usage(argv[0]);
		    } else {
			i++;
			if (!strcmp("xy", argv[i])) {
			    g[cur_graph].type = XY;
			} else if (!strcmp("fixed", argv[i])) {
			    g[cur_graph].type = XYFIXED;
			} else if (!strcmp("polar", argv[i])) {
			    g[cur_graph].type = POLAR;
			} else if (!strcmp("bar", argv[i])) {
			    g[cur_graph].type = BAR;
			} else if (!strcmp("hbar", argv[i])) {
			    g[cur_graph].type = HBAR;
			} else if (!strcmp("stackedbar", argv[i])) {
			    g[cur_graph].type = STACKEDBAR;
			} else if (!strcmp("stackedhbar", argv[i])) {
			    g[cur_graph].type = STACKEDHBAR;
			} else if (!strcmp("logx", argv[i])) {
			    g[cur_graph].type = LOGX;
			} else if (!strcmp("logy", argv[i])) {
			    g[cur_graph].type = LOGY;
			} else if (!strcmp("logxy", argv[i])) {
			    g[cur_graph].type = LOGXY;
			} else if (!strcmp("boxplot", argv[i])) {
			    g[cur_graph].type = BOXPLOT;
			} else if (!strcmp("hboxplot", argv[i])) {
			    g[cur_graph].type = HBOXPLOT;
			} else {
			    fprintf(stderr, "%s: Improper argument for -graphtype should be one of 'xy', 'logx', 'logy', 'logxy', 'bar', 'stackedbar'\n", argv[0]);
			    usage(argv[0]);
			}
		    }
		} else if (argmatch(argv[i], "-arrange", 7)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for graph arrangement\n");
			usage(argv[0]);
		    } else {
			i++;
			grows = atoi(argv[i]);
			i++;
			gcols = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-cols", 5)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for graph column arrangement\n");
			usage(argv[0]);
		    } else {
			i++;
			gcols = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-rows", 5)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for graph row arrangement\n");
			usage(argv[0]);
		    } else {
			i++;
			grows = atoi(argv[i]);
		    }
		} else if (argmatch(argv[i], "-legend", 4)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for -legend\n");
			usage(argv[0]);
		    } else {
			i++;
			if (!strcmp("load", argv[i])) {
			    loadlegend = TRUE;
			    g[cur_graph].l.active = ON;
			} else {
			    fprintf(stderr, "Improper argument for -legend\n");
			    usage(argv[0]);
			}
		    }
		} else if (argmatch(argv[i], "-rvideo", 7)) {
		    revflag = 1;
		} else if (argmatch(argv[i], "-corie", 6)) {
		    SetCorieTime(1);
		} else if (argmatch(argv[i], "-starttime", 6)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for -starttime YYMMDDHHMMSS\n");
			usage(argv[0]);
		    } else {
			i++;
			SetStartTime(1, argv[i]);
		    }
#ifdef PGSQL
		} else if (argmatch(argv[i], "-db", 3)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for db statement.\n");
			usage(argv[0]);
		    } else {
                        char *s;
			i++;
			s = strstr(argv[i], ",");
			if (s != NULL) {
			    s[0] = ' ';
                        }
			sscanf(argv[i], "%s %s", dbhost, dbase);
			printf("Connecting to [%s] %s: %s\n", argv[i], dbhost, dbase);
		    }
		} else if (argmatch(argv[i], "-dbuser", 7)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for db user,passwd statement.\n");
			usage(argv[0]);
		    } else {
                        char *s;
			i++;
			s = strstr(argv[i], ",");
			if (s != NULL) {
			    s[0] = ' ';
                        }
			sscanf(argv[i], "%s %s", dbuser, dbpasswd);
			printf("Connecting to [%s] %s: %s (%s)\n", argv[i], dbhost, dbase, dbuser);
		    }
		} else if (argmatch(argv[i], "-sql", 4)) {
		    if (i == argc) {
			fprintf(stderr, "Missing argument for sql statement.\n");
			usage(argv[0]);
		    } else {
			i++;
			ReadDBSQL(-1, -1, argv[i]);
		    }
#endif
		} else if (argmatch(argv[i], "-param", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameter file name\n");
			usage(argv[0]);
		    } else {
			strcpy(plfile, argv[i]);
			if (!getparms(cur_graph, plfile)) {
			    g[cur_graph].parmsread = FALSE;
			    fprintf(stderr, "Unable to read parameter file %s\n", plfile);
			} else {
			    g[cur_graph].parmsread = TRUE;
			}
		    }
		} else if (argmatch(argv[i], "-results", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing results file name\n");
			usage(argv[0]);
		    } else {
			strcpy(resfile, argv[i]);
/*
   *  open resfile if -results option given
 */
			if (!fexists(resfile)) {
			    if ((resfp = fopen(resfile, "w")) == NULL) {
				fprintf(stderr, "Unable to open file %s", resfile);
				exit(1);
			    }
			}
		    }
		} else if (argmatch(argv[i], "-saveall", 8)) {
		    extern char sformat[];
		    char savefile[256];
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing save file name\n");
			usage(argv[0]);
		    } else {
			strcpy(savefile, argv[i]);
			do_writesets(maxgraph, -1, 1, savefile, sformat);
		    }
		} else if (argmatch(argv[i], "-wd", 3)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameters for working directory\n");
			usage(argv[0]);
		    } else {
			if (isdir(argv[i])) {
			    strcpy(buf, argv[i]);
			    if (chdir(buf) >= 0) {
				strcpy(workingdir, buf);
			    } else {
				fprintf(stderr, "Can't change to directory %s, fatal error", buf);
				exit(1);
			    }
			} else {
			    fprintf(stderr, "%s is not a directory, fatal error\n", argv[i]);
			    exit(1);
			}
		    }
		} else if (argmatch(argv[i], "-source", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing argument for data source parameter\n");
			usage(argv[0]);
		    }
		    if (argmatch(argv[i], "pipe", 4)) {
			cursource = PIPE;
		    } else if (argmatch(argv[i], "disk", 4)) {
			cursource = DISK;
		    } else if (argmatch(argv[i], "stdin", 5)) {
			cursource = 2;
		    }
		    /* we are in a pipe */
		    if (cursource == 2) {
			if (getdata(cur_graph, "STDIN", cursource, curtype)) {
			}
		    }
		} else if (argmatch(argv[i], "-viewport", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameters for viewport setting\n");
			usage(argv[0]);
		    } else {
			g[cur_graph].v.xv1 = atof(argv[i++]);
			g[cur_graph].v.yv1 = atof(argv[i++]);
			g[cur_graph].v.xv2 = atof(argv[i++]);
			g[cur_graph].v.yv2 = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-world", 2)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing parameters for world setting\n");
			usage(argv[0]);
		    } else {
			g[cur_graph].w.xg1 = atof(argv[i++]);
			g[cur_graph].w.yg1 = atof(argv[i++]);
			g[cur_graph].w.xg2 = atof(argv[i++]);
			g[cur_graph].w.yg2 = atof(argv[i]);
		    }
		} else if (argmatch(argv[i], "-seed", 5)) {
		    i++;
		    if (i == argc) {
			fprintf(stderr, "Missing seed for srand48()\n");
			usage(argv[0]);
		    } else {
			srand48(atol(argv[i]));		/* note atol() */
		    }
		} else if (argmatch(argv[i], "-help", 5)) {
		    printf("\n%s understands all standard X Toolkit command-line options\n\n", argv[0]);
		    exit(0);
		} else if (argmatch(argv[i], "-usage", 2)) {
		    usage(argv[0]);
		} else {
		    fprintf(stderr, "No option %s\n", argv[i]);
		    usage(argv[0]);
		}
	    } else {
		if (!(i == argc)) {
		    if (getdata(cur_graph, argv[i], cursource, curtype)) {
			if (remove) {
			    unlink(argv[i]);
			}
		    }
		}		/* end if */
		remove = 0;
	    }			/* end else */
	}			/* end for */
    }				/* end if */
    for (i = 0; i < maxgraph; i++) {
	if (isactive_graph(i) && (activeset(i))) {
	    if (g[i].parmsread == FALSE) {
		autoscale_xybox(i);
		if (noautoscale[i]) {
		    switch (noautoscale[i]) {
		    case 1:
			defaulty(i, -1);
			default_axis(i, g[i].auto_type, Y_AXIS);
			default_axis(i, g[i].auto_type, ZY_AXIS);
			break;
		    case 2:
			defaultx(i, -1);
			default_axis(i, g[i].auto_type, X_AXIS);
			default_axis(i, g[i].auto_type, ZX_AXIS);
			break;
		    case 3:
			break;
		    }
		} else {
		    defaultgraph(i);
		    default_axis(i, g[i].auto_type, X_AXIS);
		    default_axis(i, g[i].auto_type, ZX_AXIS);
		    default_axis(i, g[i].auto_type, Y_AXIS);
		    default_axis(i, g[i].auto_type, ZY_AXIS);
		}
	    }
/*
 * if auto scaling type set with -a option, then scale appropriate axis
 */
	    else {
		if (autoscale[i]) {
		    switch (autoscale[i]) {
		    case 1:
			defaultx(i, -1);
			default_axis(i, g[i].auto_type, X_AXIS);
			default_axis(i, g[i].auto_type, ZX_AXIS);
			break;
		    case 2:
			defaulty(i, -1);
			default_axis(i, g[i].auto_type, Y_AXIS);
			default_axis(i, g[i].auto_type, ZY_AXIS);
			break;
		    case 3:
			defaultgraph(i);
			default_axis(i, g[i].auto_type, X_AXIS);
			default_axis(i, g[i].auto_type, ZX_AXIS);
			default_axis(i, g[i].auto_type, Y_AXIS);
			default_axis(i, g[i].auto_type, ZY_AXIS);
			break;
		    }
		}
	    }
	}
    }				/* end for */
    cg = 0;			/* default is graph 0 */
/*
 * load legend
 */
    if (loadlegend) {
	for (i = 0; i < maxgraph; i++) {
	    if (isactive_graph(i) && (activeset(i))) {
		for (j = 0; j < maxplot; j++) {
		    if (isactive(i, j)) {
			set_plotstr_string(&g[i].l.str[j], g[i].p[j].comments);
		    }
		}
	    }
	}
    }
/*
 * straighten our cursource if a pipe was used
 */
    if (cursource == 2) {
	cursource = DISK;
    }
/*
 * arrange graphs if grows & gcols set
 */
    arrange_graphs(grows, gcols);
/*
 * initialize the Hershey fonts
 */
    hselectfont(grdefaults.font);
/*
 * initialize the world and viewport
 */
    defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, islogx(cg), islogy(cg));
    viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
/*
 * if -h on command line or executed as grbatch just plot the graph and quit
 */
    if (grbatch || hardcopyflag) {
	if (hdevice == 0) {
	    fprintf(stderr,
		    "%s: Device 0 (Xlib) can't be used for batch plotting\n", argv[0]);
	    exit(1);
	}
	if (inpipe) {
	    getdata(cg, "STDIN", 2, XY);
	    inpipe = 0;
	}
	if (gotbatch && batchfile[0]) {
	    runbatch(batchfile);
	}
	if (!noprint) {
	    noerase = 1;
	    do_hardcopy();
	}
	if (resfp) {
	    fclose(resfp);
	}
	exit(0);
    }
/*
 * go window things up - do_main_loop is in x[v,m]gr.c
 */
    do_main_loop();
}

static void usage(char *progname)
{
    fprintf(stderr, "Usage of %s command line arguments: \n", progname);
    fprintf(stderr, "-maxsets   [number_of_sets]           Set the number of data sets per graph (minimum is 30)\n");
    fprintf(stderr, "-maxgraph  [number_of_graphs]         Set the number of graphs for this session (minimum is 10)\n");
    fprintf(stderr, "-autoscale [x|y|xy]                   Override any parameter file settings\n");
    fprintf(stderr, "-noauto    [x|y|xy]                   Supress autoscaling for the specified axis\n");
    fprintf(stderr, "-arrange   [rows] [cols]              Arrange the graphs in a grid rows by cols\n");
    fprintf(stderr, "-cols      [cols]\n");
    fprintf(stderr, "-rows      [rows]\n");
    fprintf(stderr, "-batch     [batch_file]               Execute batch_file on start up\n");
    fprintf(stderr, "-noask                                Assume the answer is yes to all requests - if the operation would overwrite a file, ACE/gr will do so with out prompting\n");
    fprintf(stderr, "-pipe                                 Read data from stdin on startup\n");
    fprintf(stderr, "-logwindow                            Open the log window\n");
    fprintf(stderr, "-nologwindow                          No log window, overrides resource setting\n");
    fprintf(stderr, "-device    [hardcopy device number]\n");
    fprintf(stderr, "-hardcopy                             No interactive session, just print and quit\n");
    fprintf(stderr, "-eps                                  Set the PostScript driver to write EPS\n");
    fprintf(stderr, "-log       [x|y|xy]                   Set the graph type to logarithmic\n");
    fprintf(stderr, "-legend    [load]                     Turn the graph legend on\n");
    fprintf(stderr, "-printfile [file for hardcopy output]\n");
    fprintf(stderr, "-graph     [graph number]             Set the current graph number\n");
    fprintf(stderr, "-graphsets [number_of_sets]           Set the number of data sets for the current graph\n");
    fprintf(stderr, "-graphtype [xy|bar|stackedbar|hbar|stackedhbar]  Set the type of the current graph\n");
    fprintf(stderr, "-world     [xmin ymin xmax ymax]      Set the user coordinate system for the current graph\n");
    fprintf(stderr, "-view      [xmin ymin xmax ymax]      Set the viewport for the current graph\n");
    fprintf(stderr, "-results   [results_file]             write the results from regression to results_file\n");
    fprintf(stderr, "-source    [disk|pipe|stdin]          Source of next data file\n");
    fprintf(stderr, "-param     [parameter_file]           Load parameters from parameter_file to the current graph\n");
    fprintf(stderr, "-pexec     [parameter_string]         Interpret string as a parameter setting\n");
    fprintf(stderr, "-type      [xy|xydx|xydy|xydxdx|xydydy|hilo] Set the type of the next data file\n");
    fprintf(stderr, "-ihl       [ihl_file]                 Assume data file is in IHL format (local)\n");
    fprintf(stderr, "-bp       [bp_file]                   Assume data file is in Build Points format (local)\n");
    fprintf(stderr, "-xy       [xy_file]                   Assume data file is in X Y format - sets are separated by lines containing non-numeric data\n");
    fprintf(stderr, "-nxy       [nxy_file]                  Assume data file is in X Y1 Y2 Y3 ... format\n");
    fprintf(stderr, "-xydx      [xydx_file]                Assume data file is in X Y DX format\n");
    fprintf(stderr, "-xydy      [xydy_file]                Assume data file is in X Y DY format\n");
    fprintf(stderr, "-xydxdx    [xydxdx_file]              Assume data file is in X Y DX1 DX2 format\n");
    fprintf(stderr, "-xydydy    [xydydy_file]              Assume data file is in X Y DY1 DY2 format\n");
    fprintf(stderr, "-xydxdy    [xydxdy_file]              Assume data file is in X Y DX DY format\n");
    fprintf(stderr, "-xyz       [xyz_file]                 Assume data file is in X Y Z format\n");
    fprintf(stderr, "-xyd       [xyd_file]                 Assume data file is in X Y density format\n");
    fprintf(stderr, "-xyr       [xyr_file]                 Assume data file is in X Y RADIUS format\n");
    fprintf(stderr, "-block     [block_data]               Assume data file is block data\n");
    fprintf(stderr, "-bxy       [x:y:etc.]                 Form a set from the current block data set using the current set type from columns given in the argument\n");
    fprintf(stderr, "-hilo      [hilo_file]                Assume data is in X HI LO OPEN CLOSE format\n");
    fprintf(stderr, "-boxplot   [boxplot_file]             Assume data is in X MEDIAN Y1 Y2 Y3 Y4 format\n");
#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)
    fprintf(stderr, "-netcdf    [netcdf file]              Assume data file is bnetCDF format\n");
    fprintf(stderr, "-netcdfxy  [X var name] [Y var name]  If -netcdf was used previously, read from the netCDF file, 'X var name' and 'Y var name' and create a set. If 'X var name' equals \"null\" then load the index of Y to X\n");
#endif
    fprintf(stderr, "-rvideo                               Exchange the color indices for black and white\n");
    fprintf(stderr, "-mono                                 Run %s in monochrome mode (affects the display only)\n", progname);
    fprintf(stderr, "-mssing                               Missing data value, creates breaks in lines when drawing a set. \n", progname);
    fprintf(stderr, "-seed      [seed_value]               Integer seed for random number generator\n");
    fprintf(stderr, "-GXxor                                Use xor to draw rubberband lines and graph focus markers\n");
    fprintf(stderr, "-GXinvert                             Use invert to draw rubberband lines and graph focus markers\n");
    fprintf(stderr, "-bs                                   Do backing store\n");
    fprintf(stderr, "-nobs                                 Suppress backing store\n");
    fprintf(stderr, "-dc                                   Allow double click operations on the canvas\n");
    fprintf(stderr, "-nodc                                 Disallow double click operations on the canvas\n");
    fprintf(stderr, "-maxcolors  [max_colors]              Set the number of colors to allocate (minimum is 17)\n");
    fprintf(stderr, "-redraw                               Do a redraw for refreshing the canvas when the server doesn't do backing store\n");
    fprintf(stderr, "-noredraw                             Don't do a redraw for refreshing the canvas when the server doesn't do backing store\n");
    fprintf(stderr, "-debug     [debug_level]              Set debugging options\n");
    fprintf(stderr, "-image     [image_file]               Argument is the name of an X Window dump (.xwd format)\n");
    fprintf(stderr, "-imagexy   [X] [Y]                    Arguments are the position of the image in pixels, where (0,0) is the upper left corner of the display and y increases down the screen\n");
    fprintf(stderr, "-noprint                              In batch mode, do not print\n");
    fprintf(stderr, "-usage                                This message\n");
    exit(0);
}
