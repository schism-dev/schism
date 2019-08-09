/*
 * interface to device drivers
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 */
#ifndef lint
static char RCSid[] = "$Id: draw.c,v 1.4 2006/03/13 19:46:54 pturner Exp $";
#endif

/*
   as of 2.05, the default line clipping algorithm is liang-barsky
   This was done as there were rare cases where the S-C algorithm
   would go into an infinite loop. In case there are problems with
   L-B, comment out this define to go back to the S-C clip

   PJT - 12-7-91
 */

#define LIANG_BARSKY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symdefs.h"

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define NORMAL 0
#define LOG 1
#define POLAR 2

/* maximum size of polygon for fills and patterns */
#define MAX_POLY 15000

double charsize = 1.0;
static double xg1, xg2, yg1, yg2;	/* world coordinates */
static double rxg1, rxg2, ryg1, ryg2;	/* real world coordinates (different
					 * if log plots) */
static double dxg, dyg;		/* delta x, y of world co-ordinates */
static double xv1, xv2, yv1, yv2;	/* viewpoint coordinates */
static double dxv, dyv;		/* distance covered by viewport */
static double rx, ry;		/* dxv/dxg & dyv/dyg */
static double crx, cry;		/* constants for world -> normal */
static double curx, cury;	/* current location of pen */
static int color;		/* current color */
static int lines;		/* current linestyle */
static int linew;		/* current line width */
static int ylabpos;		/* info for positioning y-axis label */
static int devtype;		/* current device */
static int clipflag = 1;	/* clip in my_draw2 */
static int scaletypex = NORMAL;	/* NORMAL, LOG, POLAR */
static int scaletypey = NORMAL;	/* NORMAL, LOG, POLAR */
static int swappedx = FALSE;	/* if axis is in reverse */
static int swappedy = FALSE;	/* if axis is in reverse */

static int do_video = 0;	/* do a video mode wherever possible */

/*
 * globals. for device drivers
 */
double devcharsize;		/* device charsize */
int (*devsetcolor) ();		/* device set color */
int (*devsetfont) ();		/* device set font */
int (*devsetline) ();		/* device set line style */
int (*devsetlinew) ();		/* device set line width */
int (*devsetpat) ();		/* device set fill pattern */
void (*devdrawarc) ();		/* device arc routine */
void (*devfillarc) ();		/* device fill arc routine */
void (*devdrawellipse) ();	/* device ellipse routine */
void (*devfillellipse) ();	/* device ellipse arc routine */
void (*devfill) ();		/* device fill polygon with a pattern */
void (*devfillcolor) ();	/* device fill polygon with color */
int (*devconvx) ();		/* viewport to device mapping for x */
int (*devconvy) ();		/* viewport to device mapping for y */
int (*vector) ();		/* device draw line */
int (*devwritestr) ();		/* device write string */
int (*devdrawtic) ();		/* device draw tic */
int (*devleavegraphics) ();	/* device exit */
int devsymsize;			/* device default symbol size */
int devxticl;			/* device tic length for x-axis */
int devyticl;			/* device tic length for y-axis */
int devarrowlength;		/* device arrowlength */
int devfontsrc;			/* device or hershey font selector */
int devwidthmm;			/* device width in mm */
int devheightmm;		/* device height in mm */
int devwidth;			/* device number of points in width */
int devheight;			/* device number of points in height */
int devoffsx;			/* device offset in x */
int devoffsy;			/* device offset in y */
int devorient;			/* = 1 if y starts at top and increases
				 * downward (terminal devices) */

static int savexticl;		/* save device tic length for x-axis */
static int saveyticl;		/* save device tic length for y-axis */

static int curfontd = 0;	/* current font */
static int current_pattern;
static double hdelta;

void my_move2();
void writestr();

extern int stringextentx();	/* in chersh.c */
extern int stringextenty();

double xconv(), yconv();
void device2view(), view2world();

void drawpolysym(double *x, double *y, int len, int sym,
		 int skip, int fill, double size);

/*
 * device2world - given (x,y) in screen coordinates, return the world coordinates
 *             in (wx,wy) used for the display only
 */
void device2world(x, y, wx, wy)
int x, y;
double *wx, *wy;
{
    double vx, vy;

    device2view(x, y, &vx, &vy);
    view2world(vx, vy, wx, wy);
}

/*
 * device2view - given (x,y) in screen coordinates, return the viewport coordinates
 */
void device2view(x, y, vx, vy)
int x, y;
double *vx, *vy;
{
    *vx = (double) x / devwidth;
    if (devorient) {
	*vy = (double) (y - devheight) / (-devheight);
    } else {
	*vy = (double) y / devheight;
    }
}

/*
 * view2world - given (vx,vy) in viewport coordinates, return world coordinates
 *            in (x,y)
 */
void view2world(vx, vy, x, y)
double vx, vy;
double *x, *y;
{
    *x = (vx - crx) / rx;
    if (scaletypex == LOG) {
	*x = pow(10.0, *x);
    }
    *y = (vy - cry) / ry;
    if (scaletypey == LOG) {
	*y = pow(10.0, *y);
    }
}

/*
 * world2deviceabs - given world coordinates (wx,wy) return the device coordinates,
 *              for the display only
 */
int world2deviceabs(wx, wy, x, y)
double wx, wy;
int *x, *y;
{
    if (scaletypex == LOG) {
	if (wx > 0.0) {
	    wx = log10(wx);
	} else {
	    return 0;
	}
    }
    *x = (rx * wx + crx) * devwidth;
    if (scaletypey == LOG) {
	if (wy > 0.0) {
	    wy = log10(wy);
	} else {
	    return 0;
	}
    }
    if (devorient) {
	*y = devheight - (ry * wy + cry) * devheight;
    } else {
	*y = (ry * wy + cry) * devheight;
    }
    return 1;
}

/*
 * world2device - given world coordinates (wx,wy) return the device coordinates,
 *              for the display only
 */
int world2device(wx, wy, x, y)
double wx, wy;
int *x, *y;
{
    if (scaletypex == LOG) {
	if (wx > 0.0) {
	    wx = log10(wx);
	} else {
	    return 0;
	}
    }
    *x = (rx * wx + crx) * (devwidth + devoffsx);
    if (scaletypey == LOG) {
	if (wy > 0.0) {
	    wy = log10(wy);
	} else {
	    return 0;
	}
    }
    *y = (ry * wy + cry) * (devheight + devoffsy);
    return 1;
}

int world2view(x, y, vx, vy)
double x, y, *vx, *vy;
{
    if (scaletypex == LOG) {
	if (x > 0.0) {
	    x = log10(x);
	} else {
	    return 0;
	}
    }
    *vx = rx * x + crx;
    if (scaletypey == LOG) {
	if (y > 0.0) {
	    y = log10(y);
	} else {
	    return 0;
	}
    }
    *vy = ry * y + cry;
    return 1;
}

/*
 * map world co-ordinates to viewport
 */
double xconv(x)
double x;
{
    if (scaletypex == LOG) {
	if (x > 0.0) {
	    x = log10(x);
	} else {
	    return 0;
	}
    }
    return (rx * x + crx);
}

double yconv(y)
double y;
{
    if (scaletypey == LOG) {
	if (y > 0.0) {
	    y = log10(y);
	} else {
	    return 0;
	}
    }
    return (ry * y + cry);
}

int set_coordmap(mapx, mapy)
int mapx, mapy;
{
    scaletypex = mapx;
    scaletypey = mapy;
}

void setfixedscale(xv1, yv1, xv2, yv2, xg1, yg1, xg2, yg2)
double xv1, yv1, xv2, yv2;
double *xg1, *yg1, *xg2, *yg2;
{
    double dx, dy, dxv, dyv, dxi, dyi, scalex, scaley;

    dx = *xg2 - *xg1;
    dy = *yg2 - *yg1;
    dxv = xv2 - xv1;
    dyv = yv2 - yv1;
    dxi = dxv * devwidthmm;
    dyi = dyv * devheightmm;
    scalex = dx / dxi;
    scaley = dy / dyi;
    if (scalex > scaley) {
	*yg2 = *yg1 + scalex * dyi;
    } else {
	*xg2 = *xg1 + scaley * dxi;
    }
}

void setmapscale(mapscale, xv1, yv1, xv2, yv2, xg1, yg1, xg2, yg2)
double mapscale, xv1, yv1, xv2, yv2;
double *xg1, *yg1, *xg2, *yg2;
{
    double dx, dy, dxv, dyv, dxi, dyi, scalex, scaley;

    dx = *xg2 - *xg1;
    dy = *yg2 - *yg1;
    dxv = xv2 - xv1;
    dyv = yv2 - yv1;
    dxi = dxv * devwidthmm * 0.001;	/* convert to meters */
    dyi = dyv * devheightmm * 0.001;	/* convert to meters */
    *xg2 = *xg1 + mapscale * dxi;
    *yg2 = *yg1 + mapscale * dyi;
}

/*
 * defineworld - really should be called definewindow, defines the scaling
 *               of the plotting rectangle to be used for clipping
 */
void defineworld(x1, y1, x2, y2, mapx, mapy)
double x1, y1, x2, y2;
int mapx, mapy;
{
    rxg1 = x1;
    ryg1 = y1;
    rxg2 = x2;
    ryg2 = y2;
    set_coordmap(mapx, mapy);
    switch (scaletypex) {
    case NORMAL:
	xg1 = x1;
	xg2 = x2;
	dxg = xg2 - xg1;
	break;
    case LOG:
	xg1 = log10(x1);
	xg2 = log10(x2);
	dxg = xg2 - xg1;
	break;
    }
    switch (scaletypey) {
    case NORMAL:
	yg1 = y1;
	yg2 = y2;
	dyg = yg2 - yg1;
	break;
    case LOG:
	yg1 = log10(y1);
	yg2 = log10(y2);
	dyg = yg2 - yg1;
	break;
    }
}

/*
 * viewport - define the location of the clipping rectangle defined by
 *            defineworld on the display device
 */
void viewport(x1, y1, x2, y2)
double x1, x2, y1, y2;
{
    xv1 = x1;
    xv2 = x2;
    yv1 = y1;
    yv2 = y2;
    dxv = x2 - x1;
    dyv = y2 - y1;
    rx = dxv / dxg;
    ry = dyv / dyg;
    crx = -rx * xg1 + xv1;
    cry = -ry * yg1 + yv1;
}

/*
 * clip if clipflag = TRUE
 */
void setclipping(fl)
{
    clipflag = fl;
}

/*
 * clip lines within rectangle defined in world coordinates
 */
static int region(xend, yend)
double xend, yend;
{
    int endpoint;

    endpoint = 0;
    if (xend < xg1)
	endpoint = 1;		/* left */
    else {
	if (xend > xg2)
	    endpoint = 2;	/* right */
    }
    if (yend < yg1)
	endpoint |= 4;		/* bottom */
    else {
	if (yend > yg2)
	    endpoint |= 8;	/* top */
    }
    return (endpoint);
}

/*
 * liang-barsky line clipping, the default as of v2.05
 */
int clipt(d, n, te, tl)
double d, n;
double *te, *tl;
{
    double t;
    int ac = 1;

    if (d > 0.0) {
	t = n / d;
	if (t > *tl) {
	    ac = 0;
	} else if (t > *te) {
	    *te = t;
	}
    } else if (d < 0.0) {
	t = n / d;
	if (t < *te) {
	    ac = 0;
	} else if (t < *tl) {
	    *tl = t;
	}
    } else {
	if (n > 0.0) {
	    ac = 0;
	}
    }
    return ac;
}

static int goodpointxy(x, y)
double x, y;
{
    return ((x >= xg1) && (x <= xg2) && (y >= yg1) && (y <= yg2));
}

static int clip2d(x0, y0, x1, y1)
double *x0, *y0, *x1, *y1;
{
    double te = 0.0, tl = 1.0, dx = (*x1 - *x0), dy = (*y1 - *y0);

    if (dx == 0.0 && dy == 0.0 && goodpointxy(*x1, *y1)) {
	return 1;
    }
    if (clipt(dx, xg1 - *x0, &te, &tl)) {
	if (clipt(-dx, *x0 - xg2, &te, &tl)) {
	    if (clipt(dy, yg1 - *y0, &te, &tl)) {
		if (clipt(-dy, *y0 - yg2, &te, &tl)) {
		    if (tl < 1.0) {
			*x1 = *x0 + tl * dx;
			*y1 = *y0 + tl * dy;
		    }
		    if (te > 0.0) {
			*x0 = *x0 + te * dx;
			*y0 = *y0 + te * dy;
		    }
		    return 1;
		}
	    }
	}
    }
    return 0;
}

/*
 * my_draw2 - draw a line from the current point (curx,cury) to (x2,y2) clipping the
 *         line to the viewport providing clipflag is TRUE
 */
void my_draw2(x2, y2)
double x2, y2;
{
    double m, minverse, x, y, x1, y1, xtmp, ytmp;
    int dir, dir1, dir2, region();

    x1 = curx;
    y1 = cury;
    switch (scaletypex) {
    case NORMAL:
	break;
    case LOG:
	x2 = log10(x2);
	break;
    }
    switch (scaletypey) {
    case NORMAL:
	break;
    case LOG:
	y2 = log10(y2);
	break;
    }
    xtmp = x2;
    ytmp = y2;
    if (!clipflag) {
	switch (scaletypex) {
	case NORMAL:
	    break;
	case LOG:
	    x1 = pow(10.0, x1);
	    x2 = pow(10.0, x2);
	    break;
	}
	switch (scaletypey) {
	case NORMAL:
	    break;
	case LOG:
	    y1 = pow(10.0, y1);
	    y2 = pow(10.0, y2);
	    break;
	}
	my_move2(x1, y1);
	(*vector) ((*devconvx) (x2), (*devconvy) (y2), 1);
	curx = xtmp;
	cury = ytmp;
	return;
    }
#ifdef LIANG_BARSKY
    if (!clip2d(&x1, &y1, &x2, &y2)) {
	curx = xtmp;
	cury = ytmp;
	return;
    }
#else
    dir1 = region(x1, y1);
    dir2 = region(x2, y2);
    if (x1 != x2) {
	m = (y2 - y1) / (x2 - x1);
    }
    if (y2 != y1) {
	minverse = (x2 - x1) / (y2 - y1);
    }
    while ((dir1 != 0) || (dir2 != 0)) {
	if (dir1 & dir2) {
	    curx = xtmp;
	    cury = ytmp;
	    return;
	}
	if (dir1 == 0) {
	    dir = dir2;
	    x = x2;
	    y = y2;
	} else {
	    dir = dir1;
	    x = x1;
	    y = y1;
	}
	if (1 & dir) {
	    y = m * (xg1 - x1) + y1;
	    x = xg1;
	} else {
	    if (2 & dir) {
		y = m * (xg2 - x1) + y1;
		x = xg2;
	    } else {
		if (4 & dir) {
		    x = minverse * (yg1 - y1) + x1;
		    y = yg1;
		} else {
		    if (8 & dir) {
			x = minverse * (yg2 - y1) + x1;
			y = yg2;
		    }
		}
	    }
	}
	if (dir == dir1) {
	    x1 = x;
	    y1 = y;
	    dir1 = region(x1, y1);
	} else {
	    x2 = x;
	    y2 = y;
	    dir2 = region(x2, y2);
	}
    }
#endif
    switch (scaletypex) {
    case NORMAL:
	break;
    case LOG:
	x1 = pow(10.0, x1);
	x2 = pow(10.0, x2);
	break;
    }
    switch (scaletypey) {
    case NORMAL:
	break;
    case LOG:
	y1 = pow(10.0, y1);
	y2 = pow(10.0, y2);
	break;
    }
    my_move2(x1, y1);
    (*vector) ((*devconvx) (x2), (*devconvy) (y2), 1);
    curx = xtmp;
    cury = ytmp;
}

/*
 * my_move2 - make (x,y) the current point (curx,cury)
 */
void my_move2(x, y)
double x, y;
{
    (*vector) ((*devconvx) (x), (*devconvy) (y), 0);
    switch (scaletypex) {
    case NORMAL:
	break;
    case LOG:
	x = log10(x);
	break;
    }
    switch (scaletypey) {
    case NORMAL:
	break;
    case LOG:
	y = log10(y);
	break;
    }
    curx = x;
    cury = y;
}

#define EDGE_LEFT       0
#define EDGE_RIGHT      1
#define EDGE_ABOVE      2
#define EDGE_BELOW      3

static int inside(x, y, which)
double x, y;
int which;
{
    switch (which) {
    case EDGE_ABOVE:
	return (y <= ryg2);
    case EDGE_BELOW:
	return (y >= ryg1);
    case EDGE_LEFT:
	return (x >= rxg1);
    case EDGE_RIGHT:
	return (x <= rxg2);
    }
}

static void intersect(px, py, x1, y1, x2, y2, which)
double *px, *py;
double x1, x2, y1, y2;
int which;
{
    double m = (x2 - x1), b = 0.0;

    if (m != 0.0) {
	m = (y2 - y1) / m;
	b = y2 - m * x2;
    }
    switch (which) {
    case EDGE_ABOVE:
	*py = ryg2;
	if (m != 0.0) {
	    *px = (ryg2 - b) / m;
	} else {
	    *px = x1;
	}
	break;
    case EDGE_BELOW:
	*py = ryg1;
	if (m != 0.0) {
	    *px = (ryg1 - b) / m;
	} else {
	    *px = x1;
	}
	break;
    case EDGE_LEFT:
	*px = rxg1;
	*py = m * rxg1 + b;
	break;
    case EDGE_RIGHT:
	*px = rxg2;
	*py = m * rxg2 + b;
	break;
    }
}

static void clip_edge(px1, py1, px2, py2, n1, n2, which)
double px1[], py1[], px2[], py2[];
int n1, *n2;
int which;

{
    int i, k, n = n1;
    double px, py, pxprev, pyprev;

    i = 0;
    *n2 = 0;
    k = 0;
    pxprev = px1[n - 1];
    pyprev = py1[n - 1];
    for (i = 0; i < n; i++) {
	if (inside(px1[i], py1[i], which)) {
	    if (inside(pxprev, pyprev, which)) {
		px2[k] = px1[i];
		py2[k] = py1[i];
		k++;
	    } else {
		intersect(&px, &py, pxprev, pyprev, px1[i], py1[i], which);
		px2[k] = px;
		py2[k] = py;
		k++;
		px2[k] = px1[i];
		py2[k] = py1[i];
		k++;
	    }
	} else if (inside(pxprev, pyprev, which)) {
	    intersect(&px, &py, pxprev, pyprev, px1[i], py1[i], which);
	    px2[k] = px;
	    py2[k] = py;
	    k++;
	} else {
	    /* do nothing */
	}
	pxprev = px1[i];
	pyprev = py1[i];
    }
    *n2 = k;
}

static void polygon_clip(px1, py1, n)
double px1[], py1[];
int *n;

{
    int n1 = *n, n2, i;
    double px2[MAX_POLY], py2[MAX_POLY];

    clip_edge(px1, py1, px2, py2, n1, &n2, EDGE_ABOVE);
    clip_edge(px2, py2, px1, py1, n2, &n1, EDGE_BELOW);
    clip_edge(px1, py1, px2, py2, n1, &n2, EDGE_LEFT);
    clip_edge(px2, py2, px1, py1, n2, &n1, EDGE_RIGHT);
    *n = n1;
}

int setpattern(k)
int k;
{
    int savepat = current_pattern;

    (*devsetpat) (current_pattern = k);
    return (savepat);
}

void fillpattern(n, px, py)
int n;
double px[], py[];

{
    int i;
    int ipx[MAX_POLY], ipy[MAX_POLY];
    double pxtmp[MAX_POLY], pytmp[MAX_POLY];

    if (clipflag) {
	for (i = 0; i < n; i++) {
	    pxtmp[i] = px[i];
	    pytmp[i] = py[i];
	}
	polygon_clip(pxtmp, pytmp, &n);
	if (n == 0) {
	    return;
	}
    }
    for (i = 0; i < n; i++) {
	ipx[i] = (*devconvx) (pxtmp[i]);
	ipy[i] = (*devconvy) (pytmp[i]);
    }
    ipx[n] = ipx[0];
    ipy[n] = ipy[0];
    n++;
    (*devfill) (n, ipx, ipy);
}

void fillcolor(n, px, py)
int n;
double px[], py[];

{
    int i;
    int ipx[MAX_POLY], ipy[MAX_POLY];
    double pxtmp[MAX_POLY], pytmp[MAX_POLY];

    if (devtype == 77) {
	for (i = 0; i < n; i++) {
	    ipx[i] = (*devconvx) (px[i]);
	    ipy[i] = (*devconvy) (py[i]);
	}
	(*devfillcolor) (n, ipx, ipy);
    } else {
	if (clipflag) {
	    for (i = 0; i < n; i++) {
		pxtmp[i] = px[i];
		pytmp[i] = py[i];
	    }
	    polygon_clip(pxtmp, pytmp, &n);
	    if (n == 0) {
		return;
	    }
	}
	for (i = 0; i < n; i++) {
	    ipx[i] = (*devconvx) (pxtmp[i]);
	    ipy[i] = (*devconvy) (pytmp[i]);
	}
	ipx[n] = ipx[0];
	ipy[n] = ipy[0];
	n++;
	(*devfillcolor) (n, ipx, ipy);
    }
}

void fillrectcolor(x1, y1, x2, y2)
double x1, y1, x2, y2;

{
    int i, n = 4;
    int ipx[5], ipy[5];
    double pxtmp[5], pytmp[5];

    pxtmp[0] = x1;
    pytmp[0] = y1;
    pxtmp[1] = x1;
    pytmp[1] = y2;
    pxtmp[2] = x2;
    pytmp[2] = y2;
    pxtmp[3] = x2;
    pytmp[3] = y1;

    if (clipflag) {
	polygon_clip(pxtmp, pytmp, &n);
	if (n == 0) {
	    return;
	}
    }
    for (i = 0; i < n; i++) {
	ipx[i] = (*devconvx) (pxtmp[i]);
	ipy[i] = (*devconvy) (pytmp[i]);
    }
    ipx[n] = ipx[0];
    ipy[n] = ipy[0];
    n++;
    (*devfillcolor) (n, ipx, ipy);
}

void fillrectpat(x1, y1, x2, y2)
double x1, y1, x2, y2;

{
    int i, n = 4;
    int ipx[5], ipy[5];
    double pxtmp[5], pytmp[5];

    pxtmp[0] = x1;
    pytmp[0] = y1;
    pxtmp[1] = x1;
    pytmp[1] = y2;
    pxtmp[2] = x2;
    pytmp[2] = y2;
    pxtmp[3] = x2;
    pytmp[3] = y1;

    if (clipflag) {
	polygon_clip(pxtmp, pytmp, &n);
	if (n == 0) {
	    return;
	}
    }
    for (i = 0; i < n; i++) {
	ipx[i] = (*devconvx) (pxtmp[i]);
	ipy[i] = (*devconvy) (pytmp[i]);
    }
    ipx[n] = ipx[0];
    ipy[n] = ipy[0];
    n++;
    (*devfill) (n, ipx, ipy);
}

/*
 * setfont - make f the current font to use for writing strings
 */
int setfont(f)
int f;
{
    int itmp = curfontd;

    curfontd = f;
    (*devsetfont) (f);
    return (itmp);
}

/*
 * rect - draw a rectangle using the current color and linestyle
 */
void rect(x1, y1, x2, y2)
double x1, x2, y1, y2;
{
    my_move2(x1, y1);
    my_draw2(x1, y2);
    my_draw2(x2, y2);
    my_draw2(x2, y1);
    my_draw2(x1, y1);
}

void rectstr(x1, y1, x2, y2, px, py)
double x1, x2, y1, y2;
int px, py;
{
    char buf[28];
    int ifx = stringextentx(devcharsize, "M");
    int ify = stringextenty(devcharsize, "M");
    if (px >= 0) {
	sprintf(buf, "%.*lf", px, x1);
	(*devwritestr) ((*devconvx) (x1), (*devconvy) (y1) - ify, 0, buf, 2, 1);
	sprintf(buf, "%.*lf", px, x2);
	(*devwritestr) ((*devconvx) (x2), (*devconvy) (y1) - ify, 0, buf, 2, 1);
    }
    if (py >= 0) {
	sprintf(buf, "%.*lf", py, y1);
	(*devwritestr) ((*devconvx) (x1) - ifx, (*devconvy) (y1), 0, buf, 1, 1);
	sprintf(buf, "%.*lf", py, y2);
	(*devwritestr) ((*devconvx) (x1) - ifx, (*devconvy) (y2), 0, buf, 1, 1);
    }
}

/*
 * symok - return TRUE if (x,y) lies within the currently defined window in
 *         world coordinates. used to clip symbols that bypass the clipping
 *         done in my_draw2.
 */
int symok(x, y)
double x, y;
{
    if (clipflag) {
	return ((x >= rxg1) && (x <= rxg2) && (y >= ryg1) && (y <= ryg2));
    } else {
	return 1;
    }
}

/*
 * histbox - draw a box of width 2*hdelta centered at x with height y starting
 *           from y = 0.0
 */
static void xhistbox(x, y)
double x, y;
{
    double tmpx[4];
    double tmpy[4];
    int i;
    int hist_pattern = 5;

    tmpx[0] = x - hdelta;
    tmpy[0] = 0.0;
    tmpx[1] = x - hdelta;
    tmpy[1] = y;
    tmpx[2] = x + hdelta;
    tmpy[2] = y;
    tmpx[3] = x + hdelta;
    tmpy[3] = 0.0;
    for (i = 0; i < 4; i++) {
	if (tmpx[i] < rxg1)
	    tmpx[i] = rxg1;
	else if (tmpx[i] > rxg2)
	    tmpx[i] = rxg2;
	if (tmpy[i] < ryg1)
	    tmpy[i] = ryg1;
	else if (tmpy[i] > ryg2)
	    tmpy[i] = ryg2;
    }
    /*
     * setpattern(hist_pattern); fillpattern(4, tmpx, tmpy); setpattern(0);
     */
    my_move2(tmpx[0], tmpy[0]);
    for (i = 0; i < 4; i++) {
	my_draw2(tmpx[(i + 1) % 4], tmpy[(i + 1) % 4]);
    }
}

/*
 * histbox3 - draw a box of width 2*hdelta centered at y with length x starting
 *           from x = 0.0
 */
static void yhistbox(x, y)
double x, y;
{
    int i;

    double tmpx[4];
    double tmpy[4];

    tmpy[0] = y - hdelta;
    tmpx[0] = 0.0;
    tmpy[1] = y - hdelta;
    tmpx[1] = x;
    tmpy[2] = y + hdelta;
    tmpx[2] = x;
    tmpy[3] = y + hdelta;
    tmpx[3] = 0.0;
    for (i = 0; i < 4; i++) {
	if (tmpx[i] < xg1)
	    tmpx[i] = xg1;
	else if (tmpx[i] > xg2)
	    tmpx[i] = xg2;
	if (tmpy[i] < yg1)
	    tmpy[i] = yg1;
	else if (tmpy[i] > yg2)
	    tmpy[i] = yg2;
    }
    /*
     * setpattern(current_pattern); fillpattern(4, tmpx, tmpy);
     */
    my_move2(tmpx[0], tmpy[0]);
    for (i = 0; i < 4; i++) {
	my_draw2(tmpx[(i + 1) % 4], tmpy[(i + 1) % 4]);
    }
}

/*
 * lengthpoly - return the length of a polyline in device coords
 */
int lengthpoly(x, y, n)
double x[], y[];
int n;

{
    int i;
    double dist = 0.0, xtmp, ytmp;

    for (i = 1; i < n; i++) {
	xtmp = (*devconvx) (x[i]) - (*devconvx) (x[i - 1]);
	ytmp = (*devconvy) (y[i]) - (*devconvy) (y[i - 1]);
	dist = dist + hypot(xtmp, ytmp);
    }
    return ((int) dist);
}

/*
 * drawpoly - draw a connected line in the current color and linestyle
 *            with nodes given by (x[],y[])
 */
void drawpoly(x, y, n)
double x[], y[];
int n;

{
    int i;

    my_move2(x[0], y[0]);
    for (i = 1; i < n; i++) {
	my_draw2(x[i], y[i]);
    }
}

/*
 * drawpolyseg - draw segments, treating each successive pairs of points
 *               as a line segment
 */
void drawpolyseg(x, y, n)
double x[], y[];
int n;

{
    int i, itmp;

    for (i = 0; i < (n / 2); i++) {
	itmp = i * 2;
	my_move2(x[itmp], y[itmp]);
	my_draw2(x[itmp + 1], y[itmp + 1]);
    }
}

void openclose(x, y1, y2, ebarlen, xy)
double x, y1, y2, ebarlen;
int xy;
{
    int ilen;

    if (xy) {
	ilen = stringextentx(ebarlen * devcharsize, "M");
	if (symok(x, y1)) {
	    (*vector) ((*devconvx) (x), (*devconvy) (y1), 0);
	    (*vector) ((*devconvx) (x) - ilen, (*devconvy) (y1), 1);
	}
	if (symok(x, y2)) {
	    (*vector) ((*devconvx) (x), (*devconvy) (y2), 0);
	    (*vector) ((*devconvx) (x) + ilen, (*devconvy) (y2), 1);
	}
    } else {
	ilen = stringextenty(ebarlen * devcharsize, "M");
	if (symok(y1, x)) {
	    (*vector) ((*devconvx) (y1), (*devconvy) (x), 0);
	    (*vector) ((*devconvx) (y1), (*devconvy) (x) - ilen, 1);
	}
	if (symok(y2, x)) {
	    (*vector) ((*devconvx) (y2), (*devconvy) (x), 0);
	    (*vector) ((*devconvx) (y2), (*devconvy) (x) + ilen, 1);
	}
    }
}

void errorbar(x, y, ebarlen, xy)
double x, y, ebarlen;
int xy;
{
    int ilen;

    if (symok(x, y)) {
	if (xy) {
	    ilen = stringextentx(ebarlen * devcharsize, "M");
	    (*vector) ((*devconvx) (x) - ilen, (*devconvy) (y), 0);
	    (*vector) ((*devconvx) (x) + ilen, (*devconvy) (y), 1);
	} else {
	    ilen = stringextenty(ebarlen * devcharsize, "M");
	    (*vector) ((*devconvx) (x), (*devconvy) (y) - ilen, 0);
	    (*vector) ((*devconvx) (x), (*devconvy) (y) + ilen, 1);
	}
    }
}

/*
 * make the current color col
 */
int setcolor(col)
int col;
{
    int scol = color;

    color = (*devsetcolor) (col);
    return scol;
}

/*
 * make the current linestyle style
 */
int setlinestyle(style)
int style;
{
    int slin = lines;

    lines = (*devsetline) (style);
    return slin;
}

/*
 * make the current line width wid
 */
int setlinewidth(wid)
int wid;
{
    int slinw = linew;

    if (devsetlinew != NULL) {
	linew = (*devsetlinew) (wid);
    }
    return slinw;
}


/*
 * drawtic - interface to device driver routines, done this way
 *           as there were problems with low resolution devices
 *           (this is probably not necessary now)
 *      dir = 0 = draw in up (for x-axis) or right (for y-axis)
 *      axis = 0 = draw x-axis tic,   1 = y-axis
 */
void drawtic(x, y, dir, axis)
double x, y;
int dir, axis;
{
    (*devdrawtic) ((*devconvx) (x), (*devconvy) (y), dir, axis);
}

/*
 * set the current character size to size
 */
double setcharsize(size)
double size;
{
    double s = charsize;
    if (devtype == 1 || devtype == 2) {
	pssetfontsize(size);
    }
    charsize = size;
    return s;
}

/*
 * set the current length of the tick marks
 */

void setticksize(sizex, sizey)
double sizex, sizey;
{
    devxticl = (int) savexticl *sizex;
    devyticl = (int) saveyticl *sizey;
}

/*
 * writestr - user interface to the current device text drawing routine
 */
void writestr(x, y, dir, just, s)
double x, y;
int dir, just;
char *s;
{
    if (s == NULL) {
	return;
    }
    if (strlen(s) == 0) {
	return;
    }
    (*devwritestr) ((*devconvx) (x), (*devconvy) (y), dir, s, just, devtype);
}

/*
 * writestrbox - user interface to the current device text drawing routine
 */
void writestrbox(x, y, dir, just, s)
double x, y;
int dir, just;
char *s;
{
    int ilenx;
    int ileny;
    int ix, iy;
    int ipx[5], ipy[5];
    int c;

    if (s == NULL) {
	return;
    }
    if (strlen(s) == 0) {
	return;
    }
    ilenx = stringextentx(charsize * devcharsize, s);
    ileny = stringextenty(charsize * devcharsize, s);
    ix = (*devconvx) (x);
    iy = (*devconvy) (y);
    ipx[0] = ix;
    ipy[0] = iy - ileny;
    ipx[1] = ix + ilenx;
    ipy[1] = iy - ileny;
    ipx[2] = ix + ilenx;
    ipy[2] = iy + ileny;
    ipx[3] = ix;
    ipy[3] = iy + ileny;
    ipx[4] = ix;
    ipy[4] = iy - ileny;
    c = setcolor(5);
    (*devfillcolor) (5, ipx, ipy);
    setcolor(1);
    (*vector) (ipx[0], ipy[0], 0);
    (*vector) (ipx[1], ipy[1], 1);
    (*vector) (ipx[2], ipy[2], 1);
    (*vector) (ipx[3], ipy[3], 1);
    (*vector) (ipx[0], ipy[0], 1);
    setcolor(c);
    (*devwritestr) (ix, iy, dir, s, just, devtype);
}

void writesymstr(double x, double y, int rot, int just, char *s,
		 int sym, int color, int fill, int loc, double symsize)
{
    int c, ifudgex, ifudgey;
    if (s == NULL) {
	return;
    }
    c = setcolor(color);
    drawpolysym(&x, &y, 1, sym, 0, fill, symsize);
    setcolor(c);
    ifudgex = stringextentx(devcharsize, "M");
    ifudgey = stringextenty(devcharsize, "My");
    switch (loc) {
    case 0:			/* right */
	(*devwritestr) ((*devconvx) (x) + ifudgex, (*devconvy) (y), rot, s, 0, devtype);
	break;
    case 1:			/* left */
	(*devwritestr) ((*devconvx) (x) - ifudgex, (*devconvy) (y), rot, s, 1, devtype);
	break;
    case 2:			/* above */
	(*devwritestr) ((*devconvx) (x), (*devconvy) (y) + ifudgey, rot, s, 2, devtype);
	break;
    case 3:			/* below */
	(*devwritestr) ((*devconvx) (x), (*devconvy) (y) - ifudgey, rot, s, 2, devtype);
	break;
    default:
	(*devwritestr) ((*devconvx) (x) + ifudgex, (*devconvy) (y), rot, s, 1, devtype);
	break;
    }
}

/*
 * draw the title, subtitle
 */
void drawtitle(title, which)
char *title;
int which;

{
    int fudge, ix, iy;
    int tmp;

    fudge = 4;
    tmp = ((*devconvx) (rxg2) + (*devconvx) (rxg1)) / 2;	/* center x-axis label
								 * and title */
    if (title[0] && !which) {
	iy = (int) (3.0 * stringextenty(charsize * devcharsize, "X"));
    } else {
	iy = stringextenty(charsize * devcharsize, "Xy");
    }
    ix = stringextentx(charsize * devcharsize, title);
    (*devwritestr) (tmp, (*devconvy) (ryg2) + iy, 0, title, 2, 0);
}

/*
   draw grid lines rather than ticmarks
 */
void drawgridlines(dir, start, end, y1, y2, step, cy, ly, wy)
int dir, ly, cy, wy;
double start, end, y1, y2, step;
{
    int slin, swid, scol;

    slin = lines;
    swid = linew;
    scol = color;
    setcolor(cy);
    setlinestyle(ly);
    setlinewidth(wy);
    while (start <= end) {
	if (dir) {
	    my_move2(y1, start);
	    my_draw2(y2, start);
	} else {
	    my_move2(start, y1);
	    my_draw2(start, y2);
	}
	start += step;
    }
    setcolor(scol);
    setlinestyle(slin);
    setlinewidth(swid);
}


/*
 * draw a circle
 */
void drawcircle(xc, yc, s, f)
double xc, yc, s;
int f;
{
    int x = (*devconvx) (xc), y = (*devconvy) (yc), xm, ym, cs;

    xm = fabs((double) (x - (*devconvx) (xc + s)));
    ym = fabs((double) (y - (*devconvy) (yc + s)));
    switch (f) {
    case 0:
	(*devdrawellipse) (x, y, xm, ym);
	break;
    case 1:
	(*devfillellipse) (x, y, xm, ym);
	if (devtype == 0) {
	    (*devdrawellipse) (x, y, xm, ym);
	}
	break;
    case 2:
	cs = setcolor(0);
	(*devfillellipse) (x, y, xm, ym);
	setcolor(cs);
	(*devdrawellipse) (x, y, xm, ym);
	break;
    }
}

/*
   place symbols at the vertices of a polygon specified by x & y.
 */

double barwid = 0.01;

void symtriangle1();
void symtriangle2();
void symtriangle3();
void symtriangle4();

void symcircle(x, y, s, f)
int x, y, f;
double s;
{
    int cs;
    int side = (int) devsymsize * s;

    switch (f) {
    case 0:
	(*devdrawarc) (x, y, side);
	break;
    case 1:
	(*devfillarc) (x, y, side);
	if (devtype == 0) {
	    (*devdrawarc) (x, y, side);
	}
	break;
    case 2:
	cs = setcolor(0);
	(*devfillarc) (x, y, side);
	setcolor(cs);
	(*devdrawarc) (x, y, side);
	break;
    }
}

void symfish1(x, y, s, f)
int x, y, f;
double s;
{
    int cs;
    int side = (int) devsymsize * s;

    switch (f) {
    case 0:
	(*devdrawarc) (x, y, side);
	break;
    case 1:
	(*devfillarc) (x, y, side);
	if (devtype == 0) {
	    (*devdrawarc) (x, y, side);
	}
	break;
    case 2:
	cs = setcolor(0);
	(*devfillarc) (x, y, side);
	setcolor(cs);
	(*devdrawarc) (x, y, side);
	break;
    }
    symtriangle2(x + 2 * side, y, s, f);
}

void symfish2(x, y, s, f)
int x, y, f;
double s;
{
    int cs;
    int side = (int) devsymsize * s;

    switch (f) {
    case 0:
	(*devdrawarc) (x, y, side);
	break;
    case 1:
	(*devfillarc) (x, y, side);
	if (devtype == 0) {
	    (*devdrawarc) (x, y, side);
	}
	break;
    case 2:
	cs = setcolor(0);
	(*devfillarc) (x, y, side);
	setcolor(cs);
	(*devdrawarc) (x, y, side);
	break;
    }
    symtriangle4(x - 2 * side, y, s, f);
}

void symsquare(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s * 0.85;
    int sx[4], sy[4], sc;

    switch (f) {
    case 0:
	(*vector) (x - side, y - side, 0);
	(*vector) (x - side, y + side, 1);
	(*vector) (x + side, y + side, 1);
	(*vector) (x + side, y - side, 1);
	(*vector) (x - side, y - side, 1);
	break;
    case 1:
	sx[0] = x - side;
	sx[1] = sx[0];
	sx[2] = x + side;
	sx[3] = sx[2];
	sy[0] = y - side;
	sy[1] = y + side;
	sy[2] = sy[1];
	sy[3] = sy[0];
	(*devfillcolor) (4, sx, sy);
	if (devtype == 0) {
	    (*vector) (x - side, y - side, 0);
	    (*vector) (x - side, y + side, 1);
	    (*vector) (x + side, y + side, 1);
	    (*vector) (x + side, y - side, 1);
	    (*vector) (x - side, y - side, 1);
	}
	break;
    case 2:
	sc = setcolor(0);
	sx[0] = x - side;
	sx[1] = sx[0];
	sx[2] = x + side;
	sx[3] = sx[2];
	sy[0] = y - side;
	sy[1] = y + side;
	sy[2] = sy[1];
	sy[3] = sy[0];
	(*devfillcolor) (4, sx, sy);
	setcolor(sc);
	(*vector) (sx[0], sy[0], 0);
	(*vector) (sx[0], sy[1], 1);
	(*vector) (sx[2], sy[1], 1);
	(*vector) (sx[2], sy[0], 1);
	(*vector) (sx[0], sy[0], 1);
    }
}

void symtriangle1(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s;
    int sx[3], sy[3], sc;

    switch (f) {
    case 0:
	(*vector) (x, y + side, 0);
	(*vector) (x - side, y - side, 1);
	(*vector) (x + side, y - side, 1);
	(*vector) (x, y + side, 1);
	break;
    case 1:
	sx[0] = x;
	sx[1] = x - side;
	sx[2] = x + side;
	sy[0] = y + side;
	sy[1] = y - side;
	sy[2] = sy[1];
	(*devfillcolor) (3, sx, sy);
	if (devtype == 0) {
	    (*vector) (x, y + side, 0);
	    (*vector) (x - side, y - side, 1);
	    (*vector) (x + side, y - side, 1);
	    (*vector) (x, y + side, 1);
	}
	break;
    case 2:
	sc = setcolor(0);
	sx[0] = x;
	sx[1] = x - side;
	sx[2] = x + side;
	sy[0] = y + side;
	sy[1] = y - side;
	sy[2] = sy[1];
	(*devfillcolor) (3, sx, sy);
	setcolor(sc);
	(*vector) (sx[0], sy[0], 0);
	(*vector) (sx[1], sy[1], 1);
	(*vector) (sx[2], sy[2], 1);
	(*vector) (sx[0], sy[0], 1);
    }
}

void symtriangle2(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s;
    int sx[3], sy[3], sc;

    switch (f) {
    case 0:
	(*vector) (x - side, y, 0);
	(*vector) (x + side, y - side, 1);
	(*vector) (x + side, y + side, 1);
	(*vector) (x - side, y, 1);
	break;
    case 1:
	sx[0] = x - side;
	sx[1] = x + side;
	sx[2] = x + side;
	sy[0] = y;
	sy[1] = y - side;
	sy[2] = y + side;
	(*devfillcolor) (3, sx, sy);
	if (devtype == 0) {
	    (*vector) (x - side, y, 0);
	    (*vector) (x + side, y - side, 1);
	    (*vector) (x + side, y + side, 1);
	    (*vector) (x - side, y, 1);
	}
	break;
    case 2:
	sc = setcolor(0);
	sx[0] = x - side;
	sx[1] = x + side;
	sx[2] = x + side;
	sy[0] = y;
	sy[1] = y - side;
	sy[2] = y + side;
	(*devfillcolor) (3, sx, sy);
	setcolor(sc);
	(*vector) (sx[0], sy[0], 0);
	(*vector) (sx[1], sy[1], 1);
	(*vector) (sx[2], sy[2], 1);
	(*vector) (sx[0], sy[0], 1);
    }
}

void symtriangle3(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s;
    int sx[3], sy[3], sc;

    switch (f) {
    case 0:
	(*vector) (x - side, y + side, 0);
	(*vector) (x, y - side, 1);
	(*vector) (x + side, y + side, 1);
	(*vector) (x - side, y + side, 1);
	break;
    case 1:
	sx[0] = x - side;
	sx[1] = x;
	sx[2] = x + side;
	sy[0] = y + side;
	sy[1] = y - side;
	sy[2] = y + side;
	(*devfillcolor) (3, sx, sy);
	if (devtype == 0) {
	    (*vector) (x - side, y + side, 0);
	    (*vector) (x, y - side, 1);
	    (*vector) (x + side, y + side, 1);
	    (*vector) (x - side, y + side, 1);
	}
	break;
    case 2:
	sc = setcolor(0);
	sx[0] = x - side;
	sx[1] = x;
	sx[2] = x + side;
	sy[0] = y + side;
	sy[1] = y - side;
	sy[2] = y + side;
	(*devfillcolor) (3, sx, sy);
	setcolor(sc);
	(*vector) (sx[0], sy[0], 0);
	(*vector) (sx[1], sy[1], 1);
	(*vector) (sx[2], sy[2], 1);
	(*vector) (sx[0], sy[0], 1);
    }
}

void symtriangle4(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s;
    int sx[3], sy[3], sc;

    switch (f) {
    case 0:
	(*vector) (x - side, y + side, 0);
	(*vector) (x - side, y - side, 1);
	(*vector) (x + side, y, 1);
	(*vector) (x - side, y + side, 1);
	break;
    case 1:
	sx[0] = x - side;
	sx[1] = x - side;
	sx[2] = x + side;
	sy[0] = y + side;
	sy[1] = y - side;
	sy[2] = y;
	(*devfillcolor) (3, sx, sy);
	if (devtype == 0) {
	    (*vector) (x - side, y + side, 0);
	    (*vector) (x - side, y - side, 1);
	    (*vector) (x + side, y, 1);
	    (*vector) (x - side, y + side, 1);
	}
	break;
    case 2:
	sc = setcolor(0);
	sx[0] = x - side;
	sx[1] = x - side;
	sx[2] = x + side;
	sy[0] = y + side;
	sy[1] = y - side;
	sy[2] = y;
	(*devfillcolor) (3, sx, sy);
	setcolor(sc);
	(*vector) (sx[0], sy[0], 0);
	(*vector) (sx[1], sy[1], 1);
	(*vector) (sx[2], sy[2], 1);
	(*vector) (sx[0], sy[0], 1);
    }
}

void symdiamond(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s;
    int sx[4], sy[4], sc;

    switch (f) {
    case 0:
	(*vector) (x, y + side, 0);
	(*vector) (x - side, y, 1);
	(*vector) (x, y - side, 1);
	(*vector) (x + side, y, 1);
	(*vector) (x, y + side, 1);
	break;
    case 1:
	sx[0] = x;
	sx[1] = x - side;
	sx[2] = sx[0];
	sx[3] = x + side;
	sy[0] = y + side;
	sy[1] = y;
	sy[2] = y - side;
	sy[3] = sy[1];
	(*devfillcolor) (4, sx, sy);
	if (devtype == 0) {
	    (*vector) (x, y + side, 0);
	    (*vector) (x - side, y, 1);
	    (*vector) (x, y - side, 1);
	    (*vector) (x + side, y, 1);
	}
	break;
    case 2:
	sc = setcolor(0);
	sx[0] = x;
	sx[1] = x - side;
	sx[2] = sx[0];
	sx[3] = x + side;
	sy[0] = y + side;
	sy[1] = y;
	sy[2] = y - side;
	sy[3] = sy[1];
	(*devfillcolor) (4, sx, sy);
	setcolor(sc);
	(*vector) (sx[0], sy[0], 0);
	(*vector) (sx[1], sy[1], 1);
	(*vector) (sx[2], sy[2], 1);
	(*vector) (sx[3], sy[3], 1);
	(*vector) (sx[0], sy[0], 1);
    }
}

void symplus(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) devsymsize * s;

    (*vector) (x - side, y, 0);
    (*vector) (x + side, y, 1);
    (*vector) (x, y + side, 0);
    (*vector) (x, y - side, 1);
}

void symx(x, y, s, f)
int x, y, f;
double s;
{
    int side = (int) (devsymsize * s * 0.707);

    (*vector) (x - side, y - side, 0);
    (*vector) (x + side, y + side, 1);
    (*vector) (x - side, y + side, 0);
    (*vector) (x + side, y - side, 1);
}

void symstar(x, y, s, f)
int x, y, f;
double s;
{
}

void symsplat(x, y, s, f)
int x, y, f;
double s;
{
    symplus(x, y, s, f);
    symx(x, y, s, f);
}

/*
   case SYM_HILOX:
   case SYM_HILOY:
   case SYM_OPENCLOSEX:
   case SYM_OPENCLOSEY:
   case SYM_CLOSEX:
   case SYM_CLOSEY:
   case SYM_HILO_OPCLX:
   case SYM_HILO_OPCLX:
   case SYM_ERRORX1:
   case SYM_ERRORY1:
   case SYM_ERRORX2:
   case SYM_ERRORY2:
   case SYM_ERRORXY:
 */
void drawpolysym(double *x, double *y, int len, int sym,
		 int skip, int fill, double size)
{
    int i;
    char s[10];
    double xtmp, ytmp;

    skip += 1;
    switch (sym) {
    case SYM_DOT:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		my_move2(x[i], y[i]);
		my_draw2(x[i], y[i]);
	    }
	}
	break;
    case SYM_CIRCLE:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symcircle((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_FISH1:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symfish1((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_FISH2:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symfish2((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_SQUARE:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symsquare((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_DIAMOND:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symdiamond((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_TRIANG1:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symtriangle1((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_TRIANG2:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symtriangle2((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_TRIANG3:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symtriangle3((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_TRIANG4:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symtriangle4((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_PLUS:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symplus((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_X:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symx((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_SPLAT:
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		symsplat((*devconvx) (x[i]), (*devconvy) (y[i]), size, fill);
	    }
	}
	break;
    case SYM_IMPULSEY:
	for (i = 0; i < len; i += skip) {
	    if (rxg1 < 0.0 && rxg2 > 0.0) {
		xtmp = 0.0;
	    } else {
		xtmp = rxg1;
	    }
	    my_move2(xtmp, y[i]);
	    my_draw2(x[i], y[i]);
	}
	break;
    case SYM_IMPULSEX:
	for (i = 0; i < len; i += skip) {
	    if (ryg1 < 0.0 && ryg2 > 0.0) {
		ytmp = 0.0;
	    } else {
		ytmp = ryg1;
	    }
	    my_move2(x[i], ytmp);
	    my_draw2(x[i], y[i]);
	}
	break;
    case SYM_VERTX:
	for (i = 0; i < len; i += skip) {
	    my_move2(x[i], ryg1);
	    my_draw2(x[i], ryg2);
	}
	break;
    case SYM_VERTY:
	for (i = 0; i < len; i += skip) {
	    my_move2(rxg1, y[i]);
	    my_draw2(rxg2, y[i]);
	}
	break;
    case SYM_HISTOX:		/* histogram x */
	for (i = 0; i < len - 1; i += skip) {
	    rect(x[i], y[i], x[i + 1], 0.0);
	}
	break;
    case SYM_BARX:		/* bar x */
	hdelta = barwid * (xg2 - xg1);
	for (i = 0; i < len; i++) {
	    xhistbox(x[i], y[i]);
	}
	break;
    case SYM_HISTOY:		/* histogram y */
	for (i = 0; i < len - 1; i += skip) {
	    rect(x[i], y[i], 0.0, y[i + 1]);
	}
	break;
    case SYM_BARY:		/* bar y */
	hdelta = barwid * (yg2 - yg1);
	for (i = 0; i < len; i++) {
	    yhistbox(x[i], y[i]);
	}
	break;
    case SYM_STAIRX:		/* stairstep */
	my_move2(x[0], y[0]);
	for (i = 0; i < len - 1; i += skip) {
	    my_draw2(x[i + 1], y[i]);
	    my_draw2(x[i + 1], y[i + 1]);
	}
	break;
    case SYM_STAIRY:
	my_move2(x[0], y[0]);
	for (i = 0; i < len - 1; i += skip) {
	    my_draw2(x[i], y[i + 1]);
	    my_draw2(x[i + 1], y[i + 1]);
	}
	break;
    case SYM_LOC:		/* index of point */
	for (i = 0; i < len; i += skip) {
	    if (symok(x[i], y[i])) {
		sprintf(s, "%d", i + 1);
		writestr(x[i], y[i], 0, 0, s);
	    }
	}
	break;
    }
}

/*
 * draw the head of an arrow
 */
void draw_head(ix1, iy1, ix2, iy2, sa, type)
int ix1, iy1, ix2, iy2, sa, type;
{
    double dist, ang, dx, dy;
    int xr, yr, xl, yl, s1, s2;
    int varrx[4], varry[4];

    dx = ix2 - ix1;
    dy = iy2 - iy1;
    if (dx == 0.0 && dy == 0.0) {
/*
   errwin("Can't draw arrow, dx = dy = 0.0");
 */
	return;
    }
    dist = hypot(dx, dy);
    ang = atan2(dy, dx);
    s1 = (int) ((dist - sa) * cos(ang));
    s2 = (int) ((dist - sa) * sin(ang));
    xr = (int) (sa * cos(ang + M_PI / 2.0) / 2.0);
    yr = (int) (sa * sin(ang + M_PI / 2.0) / 2.0);
    varrx[0] = s1 + ix1 - xr;
    varry[0] = s2 + iy1 - yr;
    varrx[1] = ix2;
    varry[1] = iy2;
    xl = (int) (sa * cos(ang - M_PI / 2.0) / 2.0);
    yl = (int) (sa * sin(ang - M_PI / 2.0) / 2.0);
    varrx[2] = s1 + ix1 - xl;
    varry[2] = s2 + iy1 - yl;
    switch (type) {
    case 0:
	(*vector) (varrx[0], varry[0], 0);
	(*vector) (varrx[1], varry[1], 1);
	(*vector) (varrx[2], varry[2], 1);
	break;
    case 1:
	(*devfillcolor) (3, varrx, varry);
	break;
    case 2:
	(*vector) (varrx[0], varry[0], 0);
	(*vector) (varrx[1], varry[1], 1);
	(*vector) (varrx[2], varry[2], 1);
	(*vector) (varrx[0], varry[0], 1);
	break;
    }
}

/*
 * draw an arrow
 */
void draw_arrow(x1, y1, x2, y2, end, asize, type)
double x1, y1, x2, y2, asize;
int end;			/* 0 = none 1 = arrow at x1, y1  2 = arrow at
				 * x2, y2 3 arrow at both ends */
{
    int ix1, iy1, ix2, iy2;

    int sa = (int) (asize * devarrowlength);

    ix1 = (*devconvx) (x1);
    ix2 = (*devconvx) (x2);
    iy1 = (*devconvy) (y1);
    iy2 = (*devconvy) (y2);
    (*vector) (ix1, iy1, 0);
    (*vector) (ix2, iy2, 1);
    switch (end) {
    case 0:
	break;
    case 1:
	draw_head(ix2, iy2, ix1, iy1, sa, type);
	break;
    case 2:
	draw_head(ix1, iy1, ix2, iy2, sa, type);
	break;
    case 3:
	draw_head(ix2, iy2, ix1, iy1, sa, type);
	draw_head(ix1, iy1, ix2, iy2, sa, type);
	break;
    }
}

/*
 * flow fields are reserved for local use
 */

#define ADJ  (8.0*M_PI/180.0)

void velplt(xx, yy, u, v, vscale)
double u, v, vscale;
double xx, yy;
{
    double x, y, s, xl, yl, xr, yr, adj;
    double theta, i1, i2;
    int ix1, iy1;
    int ix2, iy2;
    double asize = 1.0;

    int sa = (int) (asize * devarrowlength);

    if (!symok(xx, yy)) {
	return;
    }
    if ((fabs(u) > 1e-10) || (fabs(v) > 1e-10)) {
	ix1 = (*devconvx) (xx);
	iy1 = (*devconvy) (yy);
	theta = atan2(v, u);
	x = u * vscale * devwidth / (double) devwidthmm;
	y = v * vscale * devheight / (double) devheightmm;
/*
printf("%lf %lf %lf %lf %lf %d %d %d %d\n", u, v, vscale, x, y, devwidth, devheight, devwidthmm, devheightmm);
*/
	ix2 = ix1 + x;
	iy2 = iy1 + y;
	s = hypot(x, y);
	(*vector) (ix1, iy1, 0);
	(*vector) (ix2, iy2, 1);
	i1 = theta + ADJ;
	i2 = theta - ADJ;
	xl = ix1 + s * 0.6 * cos(i1);
	yl = iy1 + s * 0.6 * sin(i1);
	xr = ix1 + s * 0.6 * cos(i2);
	yr = iy1 + s * 0.6 * sin(i2);
	(*vector) ((int) xl, (int) yl, 0);
	(*vector) ((int) (ix1 + x), (int) (iy1 + y), 1);
	(*vector) ((int) xr, (int) yr, 1);
/*
   draw_head(ix1, iy1, ix2, iy2, sa, 0);
 */
    }
}

/*
 * connect the ends of 2 vectors - for drawing
 * tidal ellipses
 */
void vellineplt(xx, yy, u1, v1, u2, v2, vscale)
double u1, v1, u2, v2, vscale;
double xx, yy;
{
    double x, y;
    int ix, iy;
    int ix1, iy1;
    int ix2, iy2;

    if (!symok(xx, yy)) {
	return;
    }
    if ((fabs(u1) > 1e-10) || (fabs(v1) > 1e-10)) {
	ix = (*devconvx) (xx);
	iy = (*devconvy) (yy);
	x = u1 * vscale * devwidth / (double) devwidthmm;
	y = v1 * vscale * devheight / (double) devheightmm;
	ix1 = ix + x;
	iy1 = iy + y;
	x = u2 * vscale * devwidth / (double) devwidthmm;
	y = v2 * vscale * devheight / (double) devheightmm;
	ix2 = ix + x;
	iy2 = iy + y;
	(*vector) (ix1, iy1, 0);
	(*vector) (ix2, iy2, 1);
    }
}

void fluxplt(xx, yy, u, v, vscale)
double u, v, vscale;
double xx, yy;
{
    double x, y, s, xl, yl, xr, yr, adj;
    double theta, i1, i2;
    int ix1, iy1;
    int ix2, iy2;
    double asize = 1.0;

    int sa = (int) (asize * devarrowlength);

    if (!symok(xx, yy)) {
	return;
    }
    if ((fabs(u) > 1e-10) || (fabs(v) > 1e-10)) {
	ix1 = (*devconvx) (xx);
	iy1 = (*devconvy) (yy);
	theta = atan2(v, u);
	x = u * vscale * devwidth / (double) devwidthmm;
	y = v * vscale * devheight / (double) devheightmm;
	ix2 = ix1 + x;
	iy2 = iy1 + y;
	s = hypot(x, y);
	(*vector) (ix1, iy1, 0);
	(*vector) (ix2, iy2, 1);
    }
}

void drawsym(x, y, sym, size, fill)
int x, y, sym, fill;
double size;
{
    switch (sym) {
    case SYM_DOT:
	(*vector) (x, y, 0);
	(*vector) (x, y, 1);
	break;
    case SYM_CIRCLE:
	symcircle(x, y, size, fill);
	break;
    case SYM_SQUARE:
	symsquare(x, y, size, fill);
	break;
    case SYM_DIAMOND:
	symdiamond(x, y, size, fill);
	break;
    case SYM_TRIANG1:
	symtriangle1(x, y, size, fill);
	break;
    case SYM_TRIANG2:
	symtriangle2(x, y, size, fill);
	break;
    case SYM_TRIANG3:
	symtriangle3(x, y, size, fill);
	break;
    case SYM_TRIANG4:
	symtriangle4(x, y, size, fill);
	break;
    case SYM_PLUS:
	symplus(x, y, size, fill);
	break;
    case SYM_X:
	symx(x, y, size, fill);
	break;
    case SYM_SPLAT:
	symsplat(x, y, size, fill);
	break;
    }
}

static int xm1, xm2, ym1, ym2;

/*
 * draw the legend at world coordinates (x,y)
 */
void putlegend(i, d, xlen, ylen, size, x, y, sy, ly, cy, wy, s, barflag, fill, fu, fc, fp)
int i, d, xlen, ylen, sy, ly, cy, wy, barflag, fill, fu, fc, fp;
double size, x, y;
char *s;
{
    int ipx[4], ipy[4], itmp, scol, slins, slinw, xtmp, ytmp;
    static int maxx = 0, ifudgex, ifudgey;

    xtmp = (*devconvx) (x);
    ytmp = (*devconvy) (y) - i * ylen * ifudgey;
    if (i == 0) {
	ifudgey = stringextenty(charsize * devcharsize, "N");
	ifudgex = stringextentx(charsize * devcharsize, "N");
	xm1 = xtmp;
	xm2 = xtmp + (xlen + 1) * ifudgex + stringextentx(charsize * devcharsize, s);
	ym2 = ytmp;
	ym1 = ytmp;
    } else {
	ym1 = ytmp;
	itmp = xtmp + (xlen + 1) * ifudgex + stringextentx(charsize * devcharsize, s);
	if (xm2 < itmp) {
	    xm2 = itmp;
	}
    }
    if (d) {
	return;
    }
    scol = setcolor(cy);
    slins = setlinestyle(ly);
    slinw = setlinewidth(wy);
    if (barflag) {
	ipx[0] = xtmp;
	ipy[0] = ytmp - ifudgey;
	ipx[1] = xtmp + xlen * ifudgex;
	ipy[1] = ytmp - ifudgey;
	ipx[2] = xtmp + xlen * ifudgex;
	ipy[2] = ytmp + ifudgey;
	ipx[3] = xtmp;
	ipy[3] = ytmp + ifudgey;
	if (fu) {
	    setpattern(fp);
	    (*devfill) (4, ipx, ipy);
	} else {
	    setcolor(fc);
	    (*devfillcolor) (4, ipx, ipy);
	}
	if (ly && wy) {
	    setcolor(cy);
	    setlinestyle(ly);
	    setlinewidth(wy);
	    (*vector) (ipx[0], ipy[0], 0);
	    (*vector) (ipx[1], ipy[1], 1);
	    (*vector) (ipx[2], ipy[2], 1);
	    (*vector) (ipx[3], ipy[3], 1);
	    (*vector) (ipx[0], ipy[0], 1);
	}
    } else {
	(*vector) (xtmp, ytmp, 0);
	if (ly) {
	    (*vector) (xtmp + xlen * ifudgex, ytmp, 1);
	}
	if ((sy > 0) && (sy <= SYM_SPLAT)) {
	    setlinestyle(slins);
	    if (ly) {
		drawsym(xtmp, ytmp, sy, size, fill);
		drawsym(xtmp + xlen * ifudgex, ytmp, sy, size, fill);
	    } else {
		drawsym(xtmp + xlen * ifudgex, ytmp, sy, size, fill);
	    }
	}
    }
    setcolor(scol);
    setlinestyle(slins);
    setlinewidth(slinw);
    (*devwritestr) (xtmp + (xlen + 1) * ifudgex, ytmp, 0, s, 0, 1);
}

putlegendrect(
		 fill,
		 fillusing,
		 fillcolor,
		 fillpat,
		 cy,
		 wy,
		 ly)
int fill, fillusing, fillcolor, fillpat, cy, wy, ly;
{
    int ifudgex, ifudgey;
    int ipx[4], ipy[4];
    int scol, slins, slinw;

    scol = setcolor(cy);
    slins = setlinestyle(ly);
    slinw = setlinewidth(wy);
    ifudgey = stringextenty(charsize * devcharsize, "N");
    ifudgex = stringextentx(charsize * devcharsize, "N");
    if (fill) {
	ipx[0] = xm1 - ifudgex;
	ipy[0] = ym1 - ifudgey;
	ipx[1] = xm1 - ifudgex;
	ipy[1] = ym2 + ifudgey;
	ipx[2] = xm2 + ifudgex;
	ipy[2] = ym2 + ifudgey;
	ipx[3] = xm2 + ifudgex;
	ipy[3] = ym1 - ifudgey;
	if (fillusing) {
	    setcolor(fillcolor);
	    (*devfillcolor) (4, ipx, ipy);
	} else {
	    setpattern(fillpat);
	    (*devfill) (4, ipx, ipy);
	}
    }
    setcolor(cy);
    (*vector) (xm1 - ifudgex, ym1 - ifudgey, 0);
    (*vector) (xm1 - ifudgex, ym2 + ifudgey, 1);
    (*vector) (xm2 + ifudgex, ym2 + ifudgey, 1);
    (*vector) (xm2 + ifudgex, ym1 - ifudgey, 1);
    (*vector) (xm1 - ifudgex, ym1 - ifudgey, 1);
    setcolor(scol);
    setlinestyle(slins);
    setlinewidth(slinw);
}

/*
 * draw the legend at world coordinates (x,y)
 */
void putpatlegend(i, xlen, ylen, xgap, ygap, size, x, y, val, pat, ly, cy, ib, prec, labtype, lab)
int i, xlen, ylen, xgap, ygap, pat, ly, cy, ib, prec, labtype;
double size, x, y, val;
char *lab;
{
    int ifudgex, ifudgey, scol, slin, xtmp, ytmp;
    int ix1, ix2, iy1, iy2, px[5], py[5];
    double wx1, wy1;
    char buf[128];
    static int maxx = 0;
    static double valp;

    size *= devcharsize;
    ifudgey = stringextenty(size, "N");
    ifudgex = stringextentx(size, "N");
    view2world(x, y, &wx1, &wy1);
    ix1 = (*devconvx) (wx1);
    iy1 = (*devconvy) (wy1) - i * ygap * ifudgey;
    ix2 = ix1 - ifudgex * xlen;
    iy2 = iy1 - ifudgey * ylen;
    px[0] = ix1;
    py[0] = iy1;
    px[1] = ix1;
    py[1] = iy2;
    px[2] = ix2;
    py[2] = iy2;
    px[3] = ix2;
    py[3] = iy1;
    px[4] = ix1;
    py[4] = iy1;
    setcolor(cy);
    (*devfillcolor) (5, px, py);
    (*vector) (px[0], py[0], 0);
    (*vector) (px[1], py[1], 1);
    (*vector) (px[2], py[2], 1);
    (*vector) (px[3], py[3], 1);
    (*vector) (px[4], py[4], 1);
    switch (labtype) {
    case 0:
	sprintf(buf, "< %.*lf ", prec, val);
	break;
    case 1:
	if (i == 0) {
	    sprintf(buf, "< %.*lf ", prec, val);
	    valp = val;
	} else {
	    sprintf(buf, "%.*lf - %.*lf ", prec, valp, prec, val);
	    valp = val;
	}
	break;
    case 2:
	if (lab != NULL) {
	    strcpy(buf, lab);
	}
	break;

    }
    (*devwritestr) (ix1 + (xlen + xgap) * ifudgex, iy1 - ifudgey / 2, 0, buf, 0, 1);
}

void my_doublebuffer(mode)
int mode;
{
    switch (devtype) {
    case 0:
	xlibdoublebuffer(mode);
	break;
    }
}

void my_frontbuffer(mode)
int mode;
{
    switch (devtype) {
    case 0:
	xlibfrontbuffer(mode);
	break;
    }
}

void my_backbuffer(mode)
int mode;
{
    switch (devtype) {
    case 0:
	xlibbackbuffer(mode);
	break;
    }
}

void my_swapbuffer()
{
    switch (devtype) {
    case 0:
	xlibswapbuffer();
	break;
    }
}

set_colormapdata(i, rr, gg, bb)
int i, rr, gg, bb;
{
    switch (devtype) {
    case 0:
	xlibsetcmap(i, rr, gg, bb);
	break;
    }
}

set_colormap(cm, i, rr, gg, bb)
int cm, i, rr, gg, bb;
{
    switch (cm) {
    case 0:
	break;
    case 1:
	i = i + 16;
	break;
    case 2:
	i = i + 32;
	break;
    }
    xlibsetcmap(i, rr, gg, bb);
}

/*
 * initialize the graphics device device
 * return -1 if unable to open device
 */
int initgraphics(device)
int device;
{
    int retval;

    switch (device) {
    case 0:
	retval = xlibinitgraphics(1);
	break;
    case 1:
	retval = psinitgraphics(1);
	break;
    case 2:
	retval = psinitgraphics(3);
	break;
    case 7:
	retval = gdinitgraphics(1);
	break;
    case 8:
	retval = gdinitgraphics(3);
	break;
    default:
	retval = -1;
    }
    if (retval != -1) {
	devtype = device;
	savexticl = devxticl;
	saveyticl = devyticl;
    }
    return retval;
}

void leavegraphics()
{
    (*devleavegraphics) ();
}

static int ew = 2, eh = 5;

void putelevmarker(x, y, ex, ey, emin, emax, e)
double x, y, ex, ey, emin, emax, e;
{
    int ilen, ix1, iy1, ix2, iy2, ie, ifudgex, ifudgey;
    char buf[10];
    double tx, ty, txx, tyy, te, de = emax - emin;
    ifudgey = stringextenty(devcharsize, "N");
    ifudgex = stringextentx(devcharsize, "N");
    ilen = ifudgey * eh;
    te = (e - emin) / de;
    iy1 = (int) (te * ilen);
    my_move2(x, y);
    my_draw2(ex, ey);
    symcircle((*devconvx) (x), (*devconvy) (y), 0.8, 1);
    (*vector) ((*devconvx) (ex), (*devconvy) (ey), 0);
    (*vector) ((*devconvx) (ex), (*devconvy) (ey) + ilen, 1);
    (*vector) ((*devconvx) (ex) + ifudgex, (*devconvy) (ey) + ilen, 1);
    (*vector) ((*devconvx) (ex) + ifudgex, (*devconvy) (ey), 1);
    (*vector) ((*devconvx) (ex), (*devconvy) (ey), 1);
    (*vector) ((*devconvx) (ex), (*devconvy) (ey) + iy1, 0);
    (*vector) ((*devconvx) (ex) + ifudgex, (*devconvy) (ey) + iy1, 1);
    sprintf(buf, "%.1lf", emin);
    (*devwritestr) ((*devconvx) (ex) + (int) (1.5 * ifudgex), (*devconvy) (ey), 0, buf, 0, 1);
    sprintf(buf, "%.1lf", emax);
    (*devwritestr) ((*devconvx) (ex) + (int) (1.5 * ifudgex), (*devconvy) (ey) + ilen, 0, buf, 0, 1);
    flush_pending();
}

void set_video(vflag)
int vflag;
{
    do_video = vflag;
}

void drawclock(x, y, tott, t, clock_color, clock_fillcolor)
double x, y;
double tott, t;
int clock_color, clock_fillcolor;
{
    int xf = (*devconvx) (x), yf = (*devconvy) (y), r = 2;
    int it = (int) t % (int) (tott * 3600);
    int xt, yt, ifudgex, ifudgey;
    double tt;
    char buf[20];
    if (do_video) {
	r = 4;
    }
    ifudgey = stringextenty(devcharsize, "N");
    ifudgex = stringextentx(devcharsize, "M");

    r *= ifudgex;

    tt = M_PI / 2.0 + (it / (tott * 3600.0) * 2.0 * M_PI);
    xt = (int) (r * cos(tt));
    yt = (int) (r * sin(tt));
    setcolor(clock_fillcolor);
    (*devfillarc) (xf, yf, r);
    setcolor(clock_color);
    (*devdrawarc) (xf, yf, r);
    (*vector) (xf, yf, 0);
    (*vector) (xf - xt, yf + yt, 1);

    (*devwritestr) (xf, (int) (yf + r + ifudgey), 0, "0.0", 2, 1);
    sprintf(buf, "%.1lf", tott / 2.0);
    (*devwritestr) (xf, (int) (yf - r - ifudgey), 0, buf, 2, 1);
}

void drawmapscale(mapx, mapy, maplen, buf)
double mapx, mapy, maplen;
char *buf;
{
    double wx, wy;
    int x1, y1, x2, y2;
    int ifx = stringextentx(devcharsize, "M");
    int ify = stringextenty(devcharsize, "M");
    wx = maplen + mapx;
    wy = mapy;
    x1 = (*devconvx) (mapx);
    x2 = (*devconvx) (wx);
    y1 = (*devconvy) (mapy);
    y2 = (*devconvy) (wy);
    (*vector) (x1, y1, 0);
    (*vector) ((x2 + x1) / 2, y1, 1);
    (*vector) ((x2 + x1) / 2, y1 - ify / 2, 0);
    (*vector) ((x2 + x1) / 2, y1 + ify / 2, 1);
    (*vector) (x1, y1 - ify / 2, 0);
    (*vector) (x1, y1 + ify / 2, 1);
    (*vector) (x2, y1 + ify / 2, 1);
    (*vector) (x2, y1 - ify / 2, 1);
    (*vector) (x1, y1 - ify / 2, 1);
    (*devwritestr) (x1, y1 - (int) (1.5 * ify), 0, "0", 0, 1);
    (*devwritestr) (x2, y1 - (int) (1.5 * ify), 0, buf, 0, 1);
}

void drawvellegend(velx, vely, vellenu, vellenv, vscale, buf)
double velx, vely, vellenu, vellenv, vscale;
char *buf;
{
    int ify = stringextenty(devcharsize, "M");
    setclipping(0);
    velplt(velx, vely, vellenu, vellenv, vscale);
    (*devwritestr) ((*devconvx) (velx), (*devconvy) (vely) + ify, 0, buf, 0, 1);
    setclipping(1);
}

void drawfluxlegend(velx, vely, vellenu, vellenv, vscale, buf)
double velx, vely, vellenu, vellenv, vscale;
char *buf;
{
    int ifx = stringextentx(devcharsize, "M");
    fluxplt(velx, vely, vellenu, vellenv, vscale);
    (*devwritestr) ((*devconvx) (velx) + ifx, (*devconvy) (vely), 0, buf, 0, 1);
}

void drawisollegend(legx, legy, cis, n, vgap, hgap, form, prec, lcol, cmap, nmap, vid)
double legx, legy, cis[];
int n, vgap, hgap, form, prec, lcol, cmap[], nmap, vid;
{
    char buf[50];
    int i, j;
    int ifx = stringextentx(devcharsize, "M");
    int ify = stringextenty(devcharsize, "M");
    int wx1, wy1, wx2, wy2;
    int px[4], py[4];
    for (i = 0; i < n; i++) {
	wx1 = (*devconvx) (legx);
	wy1 = (*devconvy) (legy) - i * 2 * ify;
	wx2 = (*devconvx) (legx) + 2 * ifx;
	wy2 = (*devconvy) (legy) - (i + 1) * 2 * ify;
	if (i != n - 1) {
	    px[0] = wx1;
	    py[0] = wy1;
	    px[1] = wx1;
	    py[1] = wy2;
	    px[2] = wx2;
	    py[2] = wy2;
	    px[3] = wx2;
	    py[3] = wy1;
	    setcolor(cmap[i]);
	    (*devfillcolor) (4, px, py);
	    setcolor(lcol);
	    (*vector) (px[0], py[0], 0);
	    for (j = 1; j < 5; j++) {
		(*vector) (px[j % 4], py[j % 4], 1);
	    }
	}
	switch (form) {
	case 0:
	    sprintf(buf, " %.*lf ", prec, cis[i]);
	    break;
	case 1:
	    sprintf(buf, " %.*le ", prec, cis[i]);
	    break;
	default:
	    sprintf(buf, " %.*lg ", prec, cis[i]);
	    break;
	}
	(*devwritestr) (wx2, wy1, 0, buf, 0, 1);
    }
}

void drawisollegendh(legx, legy, cis, n, vgap, hgap, form, prec, lcol, cmap, nmap, vid)
double legx, legy, cis[];
int n, vgap, hgap, form, prec, lcol, cmap[], nmap, vid;
{
    char buf[50];
    int i, j;
    int ifx = stringextentx(devcharsize, "M");
    int ify = stringextenty(devcharsize, "M");
    int wx1, wy1, wx2, wy2;
    int px[4], py[4];
    for (i = 0; i < n; i++) {
	wx1 = (*devconvx) (legx) + i * 2 * ifx;
	wy1 = (*devconvy) (legy);
	wx2 = (*devconvx) (legx) + (i + 1) * 2 * ifx;
	wy2 = (*devconvy) (legy) - 2 * ify;
	if (i != n - 1) {
	    px[0] = wx1;
	    py[0] = wy1;
	    px[1] = wx1;
	    py[1] = wy2;
	    px[2] = wx2;
	    py[2] = wy2;
	    px[3] = wx2;
	    py[3] = wy1;
	    setcolor(cmap[i]);
	    (*devfillcolor) (4, px, py);
	    setcolor(lcol);
	    (*vector) (px[0], py[0], 0);
	    for (j = 1; j < 5; j++) {
		(*vector) (px[j % 4], py[j % 4], 1);
	    }
	}
	switch (form) {
	case 0:
	    sprintf(buf, " %.*lf ", prec, cis[i]);
	    break;
	case 1:
	    sprintf(buf, " %.*le ", prec, cis[i]);
	    break;
	default:
	    sprintf(buf, " %.*lg ", prec, cis[i]);
	    break;
	}
	(*devwritestr) (wx1, wy1 + ify, 45, buf, 0, 1);
    }
}


double convt(t, units)
double t;
int units;
{
    double retval;
    switch (units) {
    case 0:
	retval = (t);		/* seconds */
	break;
    case 1:
	retval = (t / 60.0);	/* minutes */
	break;
    case 2:
	retval = (t / 3600.0);	/* hours */
	break;
    case 3:
	retval = (t / (3600.0 * 24.0));		/* days */
	break;
    case 4:
	retval = (t / (3600.0 * 24.0 * 7.0));	/* weeks */
	break;
    }
    return retval;
}

char *format_time(t, units)
double t;
int units;
{
    static char buf[256];
    switch (units) {
    case 0:
	sprintf(buf, "%0.lfs", convt(t, units));
	break;
    case 1:
	sprintf(buf, "%0.lfm", convt(t, units));
	break;
    case 2:
	sprintf(buf, "%0.lfh", convt(t, units));
	break;
    case 3:
	sprintf(buf, "%0.lfd", convt(t, units));
	break;
    case 4:
	sprintf(buf, "%0.lfw", convt(t, units));
	break;
    case 5:
	sprintf(buf, "%0.lfM", convt(t, units));
	break;
    case 6:
	sprintf(buf, "%0.lfy", convt(t, units));
	break;
    }
    return (buf);
}


void timeline(t, start, end, tdiv, units, prec, locx, locy, len, layout, c1, c2, c3)
double t, start, end, tdiv, locx, locy;
int units, len, layout, c1, c2, c3, prec;
{
    char buf[100];
    int i, nticks, tdis, px[4], py[4];
    int ticl, itic;
    double ct;
    double tstart, tstop;
    double dist;

    int ifx = stringextentx(devcharsize, "M");
    int ify = stringextenty(devcharsize, "M");
    ticl = ifx;
    len = len * ify;

    dist = end - start;

    if (tdiv == 0.0) {
	return;
    }
    nticks = (end - start) / tdiv;
    if (nticks == 0) {
	return;
    }
    tdis = len / nticks;

    tstart = start;
    tstop = end;

    setcolor(c1);

    switch (units) {
    case 0:
	sprintf(buf, "%.*lf Seconds", prec, start);
	break;
    case 1:
	sprintf(buf, "%.*lf Minutes", prec, start);
	break;
    case 2:
	sprintf(buf, "%.*lf Hours", prec, start);
	break;
    case 3:
	sprintf(buf, "%.*lf Days", prec, start);
	break;
    case 4:
	sprintf(buf, "%.*lf Weeks", prec, start);
	break;
    }
    (*devwritestr) ((*devconvx) (locx) + 2 * ticl, (*devconvy) (locy), 0, buf, 0, 1);

    setcolor(c2);

    px[0] = (*devconvx) (locx);
    py[0] = (*devconvy) (locy);
    px[1] = (*devconvx) (locx) + ticl;
    py[1] = (*devconvy) (locy);
    px[2] = (*devconvx) (locx) + ticl;
    py[2] = (*devconvy) (locy) - len;
    px[3] = (*devconvx) (locx);
    py[3] = (*devconvy) (locy) - len;

    (*devfillcolor) (4, px, py);

    setcolor(c1);

    (*vector) ((*devconvx) (locx), (*devconvy) (locy), 0);
    (*vector) ((*devconvx) (locx) + ticl, (*devconvy) (locy), 1);
    (*vector) ((*devconvx) (locx) + ticl, (*devconvy) (locy) - len, 1);
    (*vector) ((*devconvx) (locx), (*devconvy) (locy) - len, 1);
    (*vector) ((*devconvx) (locx), (*devconvy) (locy), 1);

    sprintf(buf, "%.1lf", end);
    (*devwritestr) ((*devconvx) (locx) + 2 * ticl, (*devconvy) (locy) - len, 0, buf, 0, 1);
    for (i = 1; i < nticks; i++) {
	(*vector) ((*devconvx) (locx), (*devconvy) (locy) - i * tdis, 0);
	(*vector) ((*devconvx) (locx) + ticl, (*devconvy) (locy) - i * tdis, 1);
	sprintf(buf, "%.*lf", prec, start + tdiv * i);
	(*devwritestr) ((*devconvx) (locx) + 2 * ticl, (*devconvy) (locy) - i * tdis, 0, buf, 0, 1);
    }

    ct = (convt(t, units) - tstart) / (tstop - tstart);
    if (ct < 0.0 || ct > 1.0) {
	return;
    }
    setcolor(c3);
    itic = (int) (len * ct);
    px[1] = (*devconvx) (locx) + ticl;
    py[1] = (*devconvy) (locy);
    px[2] = (*devconvx) (locx) + ticl;
    py[2] = (*devconvy) (locy) - itic;
    px[3] = (*devconvx) (locx);
    py[3] = (*devconvy) (locy) - itic;
    (*devfillcolor) (4, px, py);
}

set_saveimages(s)
int s;
{
    extern int save_images;
    save_images = s;
}

int mapisolbath[16];
int mapisolconc[16];
static int r[256];
static int g[256];
static int b[256];
extern unsigned char red[], green[], blue[];

void initialize_cms_data()
{
    int i, j, k, rr, bb, gg;
/* standard colors */
    r[0] = 255;
    g[0] = 255;
    b[0] = 255;
    r[1] = 0;
    g[1] = 0;
    b[1] = 0;
    r[2] = 255;
    g[2] = 0;
    b[2] = 0;
    r[3] = 0;
    g[3] = 255;
    b[3] = 0;
    r[4] = 0;
    g[4] = 0;
    b[4] = 255;
    r[5] = 255;
    g[5] = 255;
    b[5] = 0;
    r[6] = 188;
    g[6] = 143;
    b[6] = 143;
    r[7] = 220;
    g[7] = 220;
    b[7] = 220;
    r[8] = 148;
    g[8] = 0;
    b[8] = 211;
    r[9] = 0;
    g[9] = 255;
    b[9] = 255;
    r[10] = 255;
    g[10] = 0;
    b[10] = 211;
    r[11] = 255;
    g[11] = 158;
    b[11] = 0;
    r[12] = 114;
    g[12] = 33;
    b[12] = 188;
    r[13] = 103;
    g[13] = 7;
    b[13] = 72;
    r[14] = 72;
    g[14] = 209;
    b[14] = 204;
    r[15] = 85;
    g[15] = 192;
    b[15] = 52;
/* bathymetry */
    r[16] = 0;
    g[16] = 0;
    b[16] = 255;
    r[17] = 15;
    g[17] = 15;
    b[17] = 255;
    r[18] = 31;
    g[18] = 31;
    b[18] = 255;
    r[19] = 47;
    g[19] = 47;
    b[19] = 255;
    r[20] = 63;
    g[20] = 63;
    b[20] = 255;
    r[21] = 79;
    g[21] = 79;
    b[21] = 255;
    r[22] = 95;
    g[22] = 95;
    b[22] = 255;
    r[23] = 111;
    g[23] = 111;
    b[23] = 255;
    r[24] = 127;
    g[24] = 127;
    b[24] = 255;
    r[25] = 143;
    g[25] = 143;
    b[25] = 255;
    r[26] = 159;
    g[26] = 159;
    b[26] = 255;
    r[27] = 175;
    g[27] = 175;
    b[27] = 255;
    r[28] = 191;
    g[28] = 191;
    b[28] = 255;
    r[29] = 207;
    g[29] = 207;
    b[29] = 255;
    r[30] = 223;
    g[30] = 223;
    b[30] = 255;
    r[31] = 239;
    g[31] = 239;
    b[31] = 255;
/* concentrations */
    r[32] = 255;
    g[32] = 57;
    b[32] = 57;
    r[33] = 255;
    g[33] = 104;
    b[33] = 55;
    r[34] = 255;
    g[34] = 153;
    b[34] = 54;
    r[35] = 255;
    g[35] = 202;
    b[35] = 53;
    r[36] = 255;
    g[36] = 251;
    b[36] = 51;
    r[37] = 209;
    g[37] = 255;
    b[37] = 50;
    r[38] = 158;
    g[38] = 255;
    b[38] = 49;
    r[39] = 107;
    g[39] = 255;
    b[39] = 47;
    r[40] = 54;
    g[40] = 255;
    b[40] = 46;
    r[41] = 45;
    g[41] = 255;
    b[41] = 88;
    r[42] = 43;
    g[42] = 255;
    b[42] = 139;
    r[43] = 42;
    g[43] = 255;
    b[43] = 190;
    r[44] = 41;
    g[44] = 255;
    b[44] = 242;
    r[45] = 39;
    g[45] = 216;
    b[45] = 255;
    r[46] = 38;
    g[46] = 162;
    b[46] = 255;
    r[47] = 37;
    g[47] = 108;
    b[47] = 255;

    for (i = 0; i < 16; i++) {
	red[i] = r[i];
	green[i] = g[i];
	blue[i] = b[i];
    }
    j = 16;
    for (i = 31; i > 15; i--) {
	red[i] = r[j];
	green[i] = g[j];
	blue[i] = b[j];
	mapisolbath[i - 16] = i;
	j++;
    }
    j = 32;
    for (i = 47; i > 31; i--) {
	red[i] = r[j];
	green[i] = g[j];
	blue[i] = b[j];
	mapisolconc[i - 32] = i;
	j++;
    }
}
