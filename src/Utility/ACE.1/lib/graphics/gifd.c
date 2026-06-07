/* $Id: gifd.c,v 1.5 2005/10/12 16:49:08 pturner Exp $
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 * driver for GIF format
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "externs.h"

#include "gd.h"
#include "gdfontl.h"
#include "gdfonts.h"

extern double charsize;
extern double devcharsize;
extern int ptofile;
extern char printstr[];
extern char *curprint;		/* curprint = gd_prstr */
int maxcolors = 48;

int gifdtrans = 0, gifdinterlace = 0; /* transparency and interlace */

/*
 * spool using these
 */
#ifndef GD_PRSTR
char gd_prstr[128] = "/bin/mv ";

#else
char gd_prstr[128] = GD_PRSTR1;

#endif

static char tmpbuf[64];

#define GDXMIN 0
#define GDXMAX 1024
#define GDYMIN 0
#define GDYMAX 1024
#define DXGD 1024
#define DYGD 1024

#define MINCOLOR 0
#define MAXCOLOR 16
#define MAXLINEWIDTH 1
#define MAXLINESTYLE 14

typedef struct {
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    double charsize;
    int symsize;
    int xticl;
    int yticl;
    int arrowlength;
} deviceparms;

deviceparms dp;

deviceparms gdparms;

/* user parameters */
static int gdxmin = GDXMIN;
static int gdxmax = GDXMAX;
static int gdymin = GDYMIN;
static int gdymax = GDYMAX;
static double gdcharsize = 0.5;
static int gdsymsize = 5;
static int gdxticl = 10;
static int gdyticl = 10;
static int gdarrowlength = 15;

static int gddx = DXGD;
static int gddy = DYGD;

static int gdcolor;
static int gddmode;
static int gdfont = 0;
static int gdlinestyle;
static int gdfillpat = 0;

static int gdcolors[256];

double xconv(double x), yconv(double y);

static FILE *gdout;

static char *fname;

static int orientflag = 0;

static gdImagePtr im_in, im_out;

/*
 * the following function is for debug purposes, testing for in/out
 */
void cgif(int x, int y)
{
    if (x < 0 || y < 0 || x >= gdxmax || y >= gdxmax) {
    }
}

void gdsetdefaults(void)
{
    gdxmin = GDXMIN;
    gdxmax = GDXMAX;
    gdymin = GDYMIN;
    gdymax = GDYMAX;
    gdcharsize = 0.5;
    gdsymsize = 5;
    gdxticl = 10;
    gdyticl = 10;
    gdarrowlength = 15;
    gddx = gdxmax - gdxmin;
    gddy = gdymax - gdymin;
}

void gdsetdevice(int xmin,
		 int xmax,
		 int ymin,
		 int ymax,
		 double charsize,
		 int symsize,
		 int xticl,
		 int yticl,
		 int arrowlength)
{
    gdxmin = xmin;
    gdxmax = xmax;
    devwidth = xmax - xmin;
    gdymin = ymin;
    gdymax = ymax;
    devheight = ymax - ymin;
    gdcharsize = charsize;
    gdsymsize = symsize;
    gdxticl = xticl;
    gdyticl = yticl;
    gdarrowlength = arrowlength;
    gddx = xmax - xmin;
    gddy = ymax - ymin;
}

int gdsetmode(int mode)
{
    int i;
    static char tbuf[256];
    static int first = 1;
    static int printToStdout = 0;
    static int curfile = 0;
    char sysbuf[256];
    char *s;

    if (mode % 2) {
	if (!ptofile) {
	    sprintf(tbuf, "xmgr%06d.gif", curfile++);
	    fname = tbuf;
	} else {
	    fname = printstr;
	    if (strcmp(fname, "stdout") == 0) {
		printToStdout = 1;
	    }
	}
	im_out = gdImageCreate(gddx, gddy);
    }

    for (i = 0; i < maxcolors; i++) {
	gdcolors[i] = gdImageColorAllocate(im_out, red[i], green[i], blue[i]);
    }
    switch (mode) {
    case 3:			/* GD portrait */
	gddx = gdxmax;
	gddy = gdymax;
	orientflag = 1;
    case 1:			/* GD landscape */
	gddx = gdxmax;
	gddy = gdymax;
	break;
    case 2:			/* */
    case 4:			/* */
	if (gifdinterlace) {
	    gdImageInterlace(im_out, 1);
	}
	if (gifdtrans) {
	    gdImageColorTransparent(im_out, 0);
	}
	if (printToStdout) {
	    gdout = stdout;
	} else {
	    gdout = fopen(fname, "wb");
	}
	gdImageGif(im_out, gdout);
	if (!printToStdout) {
	    fclose(gdout);
	} else {
	    printToStdout = 0;
	}
	gdImageDestroy(im_out);

	if (!ptofile) {
/*
	    sprintf(sysbuf, "/bin/mv %s out.gif", fname);
	    system(sysbuf);
*/
	}
	orientflag = 0;
	break;
    }
    return 1;
}

static int x1 = 99999, y1 = 99999;

void drawgd(int x2, int y2, int mode)
{

    if (mode) {
	cgif(x1, gdymax - y1);
	cgif(x2, gdymax - y2);
	if (gdlinestyle == 1) {
	    gdImageLine(im_out, x1, gdymax - y1, x2, gdymax - y2, gdcolors[gdcolor]);
	} else if (gdlinestyle > 1) {
	    gdImageLine(im_out, x1, gdymax - y1, x2, gdymax - y2, gdStyled);
	}
    } else {
	if (!(x1 == x2 && y1 == y2)) {
	}
    }
    x1 = x2;
    y1 = y2;
}

int xconvgd(double x)
{
    if (orientflag) {
	return ((int) (gdymin + gddy * xconv(x)));
    } else {
	return ((int) (gdxmin + gddx * xconv(x)));
    }
}

int yconvgd(double y)
{
    if (orientflag) {
	return ((int) (gdxmin + gddx * yconv(y)));
    } else {
	return ((int) (gdymin + gddy * yconv(y)));
    }
}

void gdsetfont(int n)
{
    hselectfont(gdfont = n);
}

int gdsetcolor(int c)
{
    if (c) {
	c = (c - 1) % maxcolors + 1;
    }
    gdcolor = c;
    return c;
}

int gdsetlinewidth(int c)
{
    gdImageSetThickness(im_out, c);
    return c;
}

void gddrawtic(int x, int y, int dir, int updown)
{
    switch (dir) {
    case 0:
	switch (updown) {
	case 0:
	    drawgd(x, y, 0);
	    drawgd(x, y + devxticl, 1);
	    break;
	case 1:
	    drawgd(x, y, 0);
	    drawgd(x, y - devxticl, 1);
	    break;
	}
	break;
    case 1:
	switch (updown) {
	case 0:
	    drawgd(x, y, 0);
	    drawgd(x + devyticl, y, 1);
	    break;
	case 1:
	    drawgd(x, y, 0);
	    drawgd(x - devyticl, y, 1);
	    break;
	}
	break;
    }
}

int gdsetlinestyle(int style)
{
    char stmp[20];
    int gdStyleArray[20];

    switch (style) {
    case 1:
	gdStyleArray[0] = gdcolor;
	gdImageSetStyle(im_out, gdStyleArray, 1);
	break;
    case 2:
	gdStyleArray[0] = gdcolor;
	gdStyleArray[1] = gdcolor;
	gdStyleArray[2] = gdTransparent;
	gdStyleArray[3] = gdTransparent;
	gdImageSetStyle(im_out, gdStyleArray, 4);
	break;
    case 3:
	gdStyleArray[0] = gdcolor;
	gdStyleArray[1] = gdcolor;
	gdStyleArray[2] = gdcolor;
	gdStyleArray[3] = gdcolor;
	gdStyleArray[4] = gdTransparent;
	gdStyleArray[5] = gdTransparent;
	gdStyleArray[6] = gdTransparent;
	gdStyleArray[7] = gdTransparent;
	gdImageSetStyle(im_out, gdStyleArray, 8);
	break;
    case 4:
	gdStyleArray[0] = gdcolor;
	gdStyleArray[1] = gdcolor;
	gdStyleArray[2] = gdcolor;
	gdStyleArray[3] = gdTransparent;
	gdStyleArray[4] = gdTransparent;
	gdStyleArray[5] = gdTransparent;
	gdStyleArray[6] = gdTransparent;
	gdStyleArray[7] = gdTransparent;
	gdStyleArray[8] = gdTransparent;
	gdImageSetStyle(im_out, gdStyleArray, 9);
	break;
    case 5:
	gdStyleArray[0] = gdcolor;
	gdStyleArray[1] = gdcolor;
	gdStyleArray[2] = gdTransparent;
	gdStyleArray[3] = gdTransparent;
	gdStyleArray[4] = gdcolor;
	gdStyleArray[5] = gdcolor;
	gdStyleArray[6] = gdcolor;
	gdStyleArray[7] = gdcolor;
	gdStyleArray[8] = gdcolor;
	gdStyleArray[9] = gdcolor;
	gdStyleArray[10] = gdTransparent;
	gdStyleArray[11] = gdTransparent;
	gdImageSetStyle(im_out, gdStyleArray, 12);
	break;
    default:
	break;
    }
    return (gdlinestyle = style);
}

void dispstrgd(int x, int y, int rot, char *s, int just, int fudge)
{
    puthersh(x, y, gdcharsize * charsize, rot, just, gdcolor, drawgd, s);
}

int gdsetpat(int pat)
{
    return (gdfillpat = pat % 7);
}

int setpatgd(int pat)
{
    switch (pat) {
    case 0:
	return (0);
    case 1:
	break;
    case 2:
	break;
    case 3:
	break;
    case 4:
	break;
    case 5:
	break;
    case 6:
	break;
    default:
	return (0);
    }
    return (pat);
}

void gdfill(int n, int *px, int *py)
{
    if (n < 3) {
	return;
    }
}


void gdfillcolor(int n, int *px, int *py)
{
    gdPoint *points;
    int i, x, y, *ytmp, cnt;
    if (n < 3) {
        return;
    }
    points = (gdPoint *) malloc(n * sizeof(gdPoint));
    if (points == NULL) {
        fprintf(stderr, "error in gdfillcolor(): can't malloc ytmp\n");
	return;
    }
    for (i = 0; i < n; i++) {
        points[i].x = px[i];
	points[i].y = gdymax - py[i];
    }
    gdImageFilledPolygon(im_out, points, n, gdcolors[gdcolor]);
    free(points);
}

void gdleavegraphics(void)
{
    gdsetmode(gddmode + 1);
}

void gddrawarc(int x, int y, int r)
{
    gdImageArc(im_out, x, gdymax - y, 2 * r, 2 * r, 0, 360, gdcolor);
}

void gdfillarc(int x, int y, int r)
{
    gdImageFilledArc(im_out, x, gdymax - y, 2 * r, 2 * r, 0, 360, gdcolor, gdArc);
}

void gddrawellipse(int x, int y, int xm, int ym)
{
    gdImageArc(im_out, x, gdymax - y, xm, ym, 0, 360, gdcolor);
}

void gdfillellipse(int x, int y, int xm, int ym)
{
    gdImageFilledArc(im_out, x, gdymax - y, xm, ym, 0, 360, gdcolor, gdArc);
}

int gdinitgraphics(int dmode)
{
    gddmode = dmode;
    if (!gdsetmode(gddmode)) {
	return -1;
    }
/* TODO temporary */
    devwidthmm = (int) (gddx / 72.0 * 25.4);
    devheightmm = (int) (gddy / 72.0 * 25.4);
    devwidth = gdxmax - gdxmin;
    devheight = gdymax - gdymin;

    devconvx = xconvgd;
    devconvy = yconvgd;
    vector = drawgd;
    devwritestr = dispstrgd;
    devsetcolor = gdsetcolor;
    devsetfont = gdsetfont;
    devsetline = gdsetlinestyle;
    devsetlinew = gdsetlinewidth;
    devdrawtic = gddrawtic;
    devsetpat = gdsetpat;
    devfill = gdfill;
    devdrawarc = gddrawarc;
    devfillarc = gdfillarc;
    devfillcolor = gdfillcolor;
    devdrawellipse = gddrawellipse;
    devfillellipse = gdfillellipse;
    devleavegraphics = gdleavegraphics;
    devcharsize = gdcharsize;
    devsymsize = gdsymsize;
    devxticl = gdxticl;
    devyticl = gdyticl;
    devarrowlength = gdarrowlength;
    setfont(2);
    setcolor(1);
    setlinestyle(0);
    return 0;
}
