/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/*
 *
 * globals.h - global variables for gredit
 *
 */

/*
 * the following definitions are found  the graphics libraries
 */
extern char ps_prstr[], noprint[];

extern int use_colors;
extern int revflag;
extern int doclear;
extern int noerase;
extern int overlay;
extern int bgcolor;
extern int inwin;
extern int backingstore;
extern int win_h, win_w;
extern int save_images;
extern int save_images_type;
extern int save_images_count;
extern char save_images_fname[];

extern int mapisolbath[];
extern int mapisolconc[];

int readgrid(int gridno, char *fname);
int ReadElcirc(int flowno, char *s, int level, int start, int stop, int skip, int missing, double mval, int append);
int ReadNodeDataNew(int n, int startstep, int stopstep, int append);
int ReadXYDataNew(int n, int startstep, int stopstep, int append);
int ReadElcircDepth(int flowno, char *fname, char *zname, double depth, int start, int stop, int skip, int missing, double mval, int append);

#ifdef MAIN

#ifndef lint
static char RCSid_globals[] = "$Id: globals.h,v 1.4 2004/10/13 23:13:10 pturner Exp $";
#endif

int change_gno;
int change_type;
int hardcopyflag;
int cursortype;
int cursor_oldx;
int cursor_oldy;
char batfile[128];

int swapBytes = 0;

int debuglevel = 0;
int expert_flag = 0;

int fgcolor;
int winsetwidth, winsetheight;
int do_backingstore = 0;

TimeClock timeclock;

int verify_action = 0;		/* request verification of actions if TRUE */
int allow_dc = 1;		/* allow double click ops */

int dobackground = 1;
int first = 1;
int auto_redraw = 1;
int draw_focus_flag;
int focus_policy = CLICK;
int cursource;
int curtype;

int region_pts = 0;
int regiontype = 0;
int regionlinkto = 0;

double regionx[1200];
double regiony[1200];
int nregion;
int region_flag = 0;

int ptofile = 0;		/* flag to indicate destination of hardcopy
				 * output, ptofile = 0 means print to printer
				 * non-zero print to file */
char printstr[1024] = "pout.dat";	/* hardcopy to this file */
char batchprefix[1024] = "xmvis5";

int scrolling_islinked;
double scrollper, shexper;

int running = 0;
int device, tdevice, hdevice;

/* graph definition */
int maxgraph = MAXGRAPH;
graph *g;
int cg = 0;			/* the current graph */
double *ax, *bx, *cx, *dx;

int curaxis = 0;

/* region definition */
region rg[MAXREGION];
int nr = 0;			/* the current region */


int format_types[] = { DECIMAL, EXPONENTIAL, POWER, GENERAL,
    DDMMYY, MMDDYY, MMYY, MMDD,
    MONTHDAY, DAYMONTH, MONTHS, MONTHL, DAYOFWEEKS,
    DAYOFWEEKL, DAYOFYEAR, HMS, MMDDHMS, MMDDYYHMS,
    DEGREESLON, DEGREESMMLON, DEGREESMMSSLON, MMSSLON,
    DEGREESLAT, DEGREESMMLAT, DEGREESMMSSLAT, MMSSLAT, 0
};

plotstr pstr[MAXSTR];		/* strings */
boxtype boxes[MAXBOXES];	/* boxes */
linetype lines[MAXLINES];	/* lines */

/* default string parameters */
plotstr sysstr;
boxtype sysbox;			/* boxes */
linetype sysline;		/* lines */

plotstr timestamp;		/* timestamp */

world_stack ws[MAX_ZOOM_STACK];	/* zoom stack */
int ws_top;			/* stack pointer */
int curw;			/* for cycling through the stack */

/* animation control */
int no_display = 0;
int anim_state = 0;

/*
 * control display items on main panel
 */
int display_flags[30], display_order[30], ndisplay = NDISPLAY;
int display_strings = 1;


FILE *logfp = NULL;
int logfile = 0;

double page_per = 0.1;

int curset;
int curteanl;
int curteanlem;
int curadcirc;
int curadcircem;
int curadc3dvm;
int curadc3dem;
int curadc3d;
int curadc3dcnc;
int curela;
int curelaem;
int curvelhist;
int curdrog;
int curhist;
int curchist;
int curslice;
int curflux;
int curzoom;
int curtrans;

Isolparms curip;
int curplaceitem;

Grid grid[MAXGRIDS];
Gridt gridt[MAXGRIDTS];
Flow_tct flowf[MAXTEANL];
Flow2d flowt[MAXADCIRC];
Flow2d flowh[MAXADCIRC];
ADCIRC3D adc3d[MAXADCIRC3D];	/* ADCIRC 3D */
Scalar2d elaconc[MAXELA];
Particle_field drogues[MAXPATHLINES];
Track track[MAXTRACK];
TimeHistory hist[MAXHISTMARKERS];
TimeHistory timehist[MAXTIMEHIST];
TideStation **tidestat;
int ntidestat;
Station *sta;
int nsta;
Boundaries bounds[MAXBOUNDS];
Transect trans[MAXTRANS];

char buf[1024];			/* a string used here and there */
char *curprint;			/* the default printer */

char statusstr[1024];

#else				/* in other than main.c */

extern int change_gno;
extern int change_type;
extern int hardcopyflag;
extern int cursortype;
extern int cursor_oldx;
extern int cursor_oldy;
extern char batfile[];

extern int swapBytes;

extern int debuglevel;
extern int expert_flag;

extern int fgcolor;
extern int winsetwidth, winsetheight;
extern int do_backingstore;

extern TimeClock timeclock;

extern int verify_action;	/* request verification of actions if TRUE */
extern int allow_dc;		/* allow double click ops */

extern int dobackground;
extern int first;
extern int auto_redraw;
extern int draw_focus_flag;
extern int focus_policy;
extern int cursource;
extern int curtype;

extern int ptofile;
extern char printstr[];
extern char batchprefix[];

extern int scrolling_islinked;
extern double scrollper, shexper;

extern int running;
extern int device, tdevice, hdevice;

/* graph definition */
extern int maxgraph;
extern graph *g;
extern int cg;			/* the current graph */
extern double *ax, *bx, *cx, *dx;

extern int curaxis;

/* region definition */
extern region rg[MAXREGION];
extern int nr;			/* the current region */

extern int format_types[];

extern plotstr pstr[];		/* strings */
extern boxtype boxes[];		/* boxes */
extern linetype lines[];	/* lines */

/* default string parameters */
extern plotstr sysstr;
extern boxtype sysbox;		/* boxes */
extern linetype sysline;	/* lines */

extern plotstr timestamp;	/* timestamp */

extern world_stack ws[];	/* zoom stack */
extern int ws_top;		/* stack pointer */
extern int curw;		/* for cycling through the stack */

/* animation control */
extern int no_display;
extern int anim_state;

/*
 * control display items on main panel
 */
extern int display_flags[], display_order[], ndisplay;
extern int display_strings;

extern FILE *logfp;
extern int logfile;

extern double page_per;

extern int curset, curteanl, curteanlem;
extern int curadcirc, curadcircem, curela, curelaem, curvelhist, curdrog;
extern int curadc3dvm;
extern int curadc3dem;
extern int curadc3d;
extern int curadc3dcnc;
extern int curhist, curchist;
extern int curslice, curflux, curzoom, curtrans;

extern int region_pts;
extern int regiontype;
extern int regionlinkto;

extern Isolparms curip;
extern int curplaceitem;

extern Grid grid[];
extern Gridt gridt[];
extern Boundary boundary[];
extern Boundaries bounds[];
extern Station *sta;
extern int nsta;

extern Flow_tct flowf[];
extern Flow2d flowt[];
extern Flow2d flowh[];
extern ADCIRC3D adc3d[];
extern Transect trans[MAXTRANS];

extern Scalar2d elaconc[];

extern Particle_field drogues[];
extern Track track[];

extern TimeHistory hist[];
extern TimeHistory timehist[];
extern TideStation **tidestat;
extern int ntidestat;

extern char buf[];		/* a string used here and there */
extern char *curprint;		/* the default printer */

extern double regionx[];
extern double regiony[];
extern int nregion;
extern int region_flag;

extern char statusstr[];

#endif				/* ifdef MAIN */
