/* $Id: svg.c,v 1.2 2003/07/24 15:44:04 pturner Exp $
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 * driver for GIF format
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "externs.h"

extern double charsize;
extern double devcharsize;
extern int ptofile;
extern char printstr[];
extern char *curprint;		/* curprint = svg_prstr */
int maxcolors = 48;

int gifdtrans = 0, gifdinterlace = 0; /* transparency and interlace */

/*
 * spool using these
 */
#ifndef SVG_PRSTR
char svg_prstr[128] = "/bin/mv ";

#else
char svg_prstr[128] = SVG_PRSTR1;

#endif

static char tmpbuf[64];

#define SVGXMIN 0
#define SVGXMAX 600
#define SVGYMIN 0
#define SVGYMAX 800
#define DXSVG 600
#define DYSVG 800

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

deviceparms svgparms;

/* user parameters */
static int svgxmin = SVGXMIN;
static int svgxmax = SVGXMAX;
static int svgymin = SVGYMIN;
static int svgymax = SVGYMAX;
static double svgcharsize = 0.5;
static int svgsymsize = 5;
static int svgxticl = 10;
static int svgyticl = 10;
static int svgarrowlength = 15;

static int svgdx = DXSVG;
static int svgdy = DYSVG;

static int svgcolor;
static int svgdmode;
static int svgfont = 0;
static int svglinestyle;
static int svgfillpat = 0;

static int svgcolors[256];

double xconv(double x), yconv(double y);

static FILE *svgout;

static char *fname;

static int orientflag = 0;

static svgImagePtr im_in, im_out;

/*
 * the following function is for debug purposes, testing for in/out
 */
void cgif(int x, int y)
{
    if (x < 0 || y < 0 || x >= svgxmax || y >= svgxmax) {
    }
}

void svgsetdefaults(void)
{
    svgxmin = SVGXMIN;
    svgxmax = SVGXMAX;
    svgymin = SVGYMIN;
    svgymax = SVGYMAX;
    svgcharsize = 0.5;
    svgsymsize = 5;
    svgxticl = 10;
    svgyticl = 10;
    svgarrowlength = 15;
    svgdx = svgxmax - svgxmin;
    svgdy = svgymax - svgymin;
}

void svgsetdevice(int xmin,
		 int xmax,
		 int ymin,
		 int ymax,
		 double charsize,
		 int symsize,
		 int xticl,
		 int yticl,
		 int arrowlength)
{
    svgxmin = xmin;
    svgxmax = xmax;
    devwidth = xmax - xmin;
    svgymin = ymin;
    svgymax = ymax;
    devheight = ymax - ymin;
    svgcharsize = charsize;
    svgsymsize = symsize;
    svgxticl = xticl;
    svgyticl = yticl;
    svgarrowlength = arrowlength;
    svgdx = xmax - xmin;
    svgdy = ymax - ymin;
}

int svgsetmode(int mode)
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
	im_out = svgImageCreate(svgdx, svgdy);
    }

    for (i = 0; i < maxcolors; i++) {
	svgcolors[i] = svgImageColorAllocate(im_out, red[i], green[i], blue[i]);
    }
    switch (mode) {
    case 3:			/* SVG portrait */
	svgdx = svgxmax;
	svgdy = svgymax;
	orientflag = 1;
    case 1:			/* SVG landscape */
	svgdx = svgxmax;
	svgdy = svgymax;
	break;
    case 2:			/* */
    case 4:			/* */
	if (gifdinterlace) {
	    svgImageInterlace(im_out, 1);
	}
	if (gifdtrans) {
	    svgImageColorTransparent(im_out, 0);
	}
	if (printToStdout) {
	    svgout = stdout;
	} else {
	    svgout = fopen(fname, "wb");
	}
	svgImageGif(im_out, svgout);
	if (!printToStdout) {
	    fclose(svgout);
	} else {
	    printToStdout = 0;
	}
	svgImageDestroy(im_out);

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

void drawsvg(int x2, int y2, int mode)
{

    if (mode) {
	cgif(x1, svgymax - y1);
	cgif(x2, svgymax - y2);
	if (svglinestyle == 1) {
	} else if (svglinestyle > 1) {
	}
    } else {
	if (!(x1 == x2 && y1 == y2)) {
	}
    }
    x1 = x2;
    y1 = y2;
}

int xconvsvg(double x)
{
    if (orientflag) {
	return ((int) (svgymin + svgdy * xconv(x)));
    } else {
	return ((int) (svgxmin + svgdx * xconv(x)));
    }
}

int yconvsvg(double y)
{
    if (orientflag) {
	return ((int) (svgxmin + svgdx * yconv(y)));
    } else {
	return ((int) (svgymin + svgdy * yconv(y)));
    }
}

void svgsetfont(int n)
{
    hselectfont(svgfont = n);
}

int svgsetcolor(int c)
{
    if (c) {
	c = (c - 1) % maxcolors + 1;
    }
    svgcolor = c;
    return c;
}

int svgsetlinewidth(int c)
{
    return c;
}

void svgdrawtic(int x, int y, int dir, int updown)
{
    switch (dir) {
    case 0:
	switch (updown) {
	case 0:
	    drawsvg(x, y, 0);
	    drawsvg(x, y + devxticl, 1);
	    break;
	case 1:
	    drawsvg(x, y, 0);
	    drawsvg(x, y - devxticl, 1);
	    break;
	}
	break;
    case 1:
	switch (updown) {
	case 0:
	    drawsvg(x, y, 0);
	    drawsvg(x + devyticl, y, 1);
	    break;
	case 1:
	    drawsvg(x, y, 0);
	    drawsvg(x - devyticl, y, 1);
	    break;
	}
	break;
    }
}

int svgsetlinestyle(int style)
{
    char stmp[20];
    int svgStyleArray[20];

    switch (style) {
    case 1:
	svgStyleArray[0] = svgcolor;
	svgImageSetStyle(im_out, svgStyleArray, 1);
	break;
    case 2:
	svgStyleArray[0] = svgcolor;
	svgStyleArray[1] = svgcolor;
	svgStyleArray[2] = svgTransparent;
	svgStyleArray[3] = svgTransparent;
	svgImageSetStyle(im_out, svgStyleArray, 4);
	break;
    case 3:
	svgStyleArray[0] = svgcolor;
	svgStyleArray[1] = svgcolor;
	svgStyleArray[2] = svgcolor;
	svgStyleArray[3] = svgcolor;
	svgStyleArray[4] = svgTransparent;
	svgStyleArray[5] = svgTransparent;
	svgStyleArray[6] = svgTransparent;
	svgStyleArray[7] = svgTransparent;
	svgImageSetStyle(im_out, svgStyleArray, 8);
	break;
    case 4:
	svgStyleArray[0] = svgcolor;
	svgStyleArray[1] = svgcolor;
	svgStyleArray[2] = svgcolor;
	svgStyleArray[3] = svgTransparent;
	svgStyleArray[4] = svgTransparent;
	svgStyleArray[5] = svgTransparent;
	svgStyleArray[6] = svgTransparent;
	svgStyleArray[7] = svgTransparent;
	svgStyleArray[8] = svgTransparent;
	svgImageSetStyle(im_out, svgStyleArray, 9);
	break;
    case 5:
	svgStyleArray[0] = svgcolor;
	svgStyleArray[1] = svgcolor;
	svgStyleArray[2] = svgTransparent;
	svgStyleArray[3] = svgTransparent;
	svgStyleArray[4] = svgcolor;
	svgStyleArray[5] = svgcolor;
	svgStyleArray[6] = svgcolor;
	svgStyleArray[7] = svgcolor;
	svgStyleArray[8] = svgcolor;
	svgStyleArray[9] = svgcolor;
	svgStyleArray[10] = svgTransparent;
	svgStyleArray[11] = svgTransparent;
	svgImageSetStyle(im_out, svgStyleArray, 12);
	break;
    default:
	break;
    }
    return (svglinestyle = style);
}

void dispstrsvg(int x, int y, int rot, char *s, int just, int fudge)
{
    puthersh(x, y, svgcharsize * charsize, rot, just, svgcolor, drawsvg, s);
}

int svgsetpat(int pat)
{
    return (svgfillpat = pat % 7);
}

int setpatsvg(int pat)
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

void svgfill(int n, int *px, int *py)
{
    if (n < 3) {
	return;
    }
}

void svgfillcolor(int n, int *px, int *py)
{
    int i, x, y, *ytmp, cnt;
    if (n < 3) {
        return;
    }
    ytmp = (int *) malloc(n * sizeof(int));
    if (ytmp == NULL) {
        fprintf(stderr, "error in svgfillcolor(): can't malloc ytmp\n");
	return;
    }
    for (i = 0; i < n; i++) {
        ytmp[i] = svgymax - py[i];
    }
    polyfill(im_out, px, ytmp, n, svgcolors[svgcolor]);
    free(ytmp);
}

void svgleavegraphics(void)
{
    fclose(svgout);
    svgsetmode(svgdmode + 1);
}

void svgdrawarc(int x, int y, int r)
{
}

void svgfillarc(int x, int y, int r)
{
}

void svgdrawellipse(int x, int y, int xm, int ym)
{
}

void svgfillellipse(int x, int y, int xm, int ym)
{
}

int svginitgraphics(int dmode)
{
    svgdmode = dmode;
    if (!svgsetmode(svgdmode)) {
	return -1;
    }
/* TODO temporary */
    devwidthmm = (int) (svgdx / 72.0 * 25.4);
    devheightmm = (int) (svgdy / 72.0 * 25.4);

    devconvx = xconvsvg;
    devconvy = yconvsvg;
    vector = drawsvg;
    devwritestr = dispstrsvg;
    devsetcolor = svgsetcolor;
    devsetfont = svgsetfont;
    devsetline = svgsetlinestyle;
    devsetlinew = svgsetlinewidth;
    devdrawtic = svgdrawtic;
    devsetpat = svgsetpat;
    devfill = svgfill;
    devdrawarc = svgdrawarc;
    devfillarc = svgfillarc;
    devfillcolor = svgfillcolor;
    devdrawellipse = svgdrawellipse;
    devfillellipse = svgfillellipse;
    devleavegraphics = svgleavegraphics;
    devcharsize = svgcharsize;
    devsymsize = svgsymsize;
    devxticl = svgxticl;
    devyticl = svgyticl;
    devarrowlength = svgarrowlength;
    setfont(2);
    setcolor(1);
    setlinestyle(0);
    svgout = fopen("test.svg", "wb");
   fprintf(prstream, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(prstream, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"");
     fprintf(prstream, " \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n");
      fprintf(prstream, "<!-- generated by %s -->\n", "xmgr5");
       /* Let do the SVG-Viewer the conversion of coordinates. */
       fprintf(prstream, "<svg xml:space=\"preserve\" ");
        fprintf(prstream, "width=\"%.4fin\" height=\"%.4fin\" viewBox=\"%.4f %.4f %.4f %.4f\">\n",
		                 page_width_in, page_height_in,
			             0.0, 0.0, page_width_pp, page_height_pp);
    fprintf(prstream, " <g transform=\"translate(0,%.4f) scale(1,-1)\">\n",
	                page_height_pp);

        fprintf(prstream, " <desc>%s</desc>\n", "test SVG driver");

    return 0;
}
