/*
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 * driver for xlib for gr
 *
 */
#ifndef lint
static char RCSid[] = "$Id: xvlib.c,v 1.3 2003/10/17 02:38:47 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>

#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <Xm/Xm.h>

#include "externs.h"
#include "patterns.h"

/* publicly available variables */
Display *disp;
GC gc;
GC gcxor;
GC gcclr;
Window cwin;
XGCValues gc_val;

int pathlength = 0;
int use_colors;
int revflag = 0;
int doclear = 1;
int noerase = 1;
int overlay = 0;
int bgcolor = 0;
int inwin;
int backingstore;
int win_h, win_w;
int save_fli;
int save_images;
int save_images_type;
int save_images_count;
char save_images_fname[256];

extern double devcharsize;
extern double charsize;

XColor cmscolors[60];

#define MINCOLOR 0
#define MAXCOLOR 48
#define MAXLINEW 15

double devtoworldx(), devtoworldy();
void xlibsetcmap();

static Widget cc;		/* current canvas */
static Widget bc[10];		/* for backing store */

void flush_pending();

static int xlibcolor = 1;
static int xliblinewidth = 0;
static int xlibdmode;
static int xlibfont = 0;
static int xliblinestyle = 1;
static int doublebuff = 0;	/* no double buffering by default */
static Pixmap backbuff[10], backpix[10];
Pixmap displaybuff, backb;
static int backpixind = 0;

double xlibcharsize = 0.60;

double xconv(), yconv();

void drawxlib();
void xlibdoublebuffer();
void xlibinit_tiles();

void get_xlib_dims(w, h)
int *w, *h;
{
    Arg args;
    Dimension ww, wh;

    XtSetArg(args, XmNwidth, &ww);
    XtGetValues(cc, &args, 1);
    XtSetArg(args, XmNheight, &wh);
    XtGetValues(cc, &args, 1);
    *w = ww;
    *h = wh;
}

void set_window(cc_set)
Widget cc_set;
{
    double wx1, wx2, wy1, wy2;
    Arg args;
    Dimension ww, wh;
    int i, found;
    disp = XtDisplay(cc_set);
    cwin = XtWindow(cc_set);
    displaybuff = cwin;
    XtSetArg(args, XmNwidth, &ww);
    XtGetValues(cc_set, &args, 1);
    XtSetArg(args, XmNheight, &wh);
    XtGetValues(cc_set, &args, 1);
    win_w = ww;
    win_h = wh;
    devwidth = win_w;
    devheight = win_h;
    wx1 = DisplayWidth(disp, DefaultScreen(disp));
    wx2 = DisplayWidthMM(disp, DefaultScreen(disp));
    wy1 = DisplayHeight(disp, DefaultScreen(disp));
    wy2 = DisplayHeightMM(disp, DefaultScreen(disp));
    devwidthmm = (int) (wx2 / wx1 * win_w);
    devheightmm = (int) (wy2 / wy1 * win_h);
    displaybuff = cwin;
    cc = cc_set;
    if (backingstore) {
	found = 0;
	for (i = 0; i < 10; i++) {
	    if (cc == bc[i]) {
		found = i;
	    }
	}
	if (!found) {
	    for (i = 0; i < 10; i++) {
		if (!bc[i]) {
		    found = i;
		}
	    }
	    backpixind = found;
	    bc[found] = cc;
	} else {
	    backpixind = found;
	}
    }
}

void xlibinitgc(Widget w)
{
    disp = XtDisplay(w);
    cwin = XtWindow(w);
    gc = DefaultGC(disp, DefaultScreen(disp));
    gc_val.foreground = WhitePixel(disp, DefaultScreen(disp));
    gc_val.function = GXxor;
    gcxor = XCreateGC(disp, cwin, GCFunction | GCForeground, &gc_val);
    gc_val.foreground = WhitePixel(disp, DefaultScreen(disp));
    gc_val.function = GXcopy;
    gcclr = XCreateGC(disp, cwin, GCFunction | GCForeground, &gc_val);
}

void xlibsetbackingstore(int w)
{
    backingstore = w;
}

static void xlibinit()
{
    double wx1, wx2, wy1, wy2;
    static int inc = 1;

    Arg args;
    Dimension ww, wh;

    disp = XtDisplay(cc);
    cwin = XtWindow(cc);
    XtSetArg(args, XmNwidth, &ww);
    XtGetValues(cc, &args, 1);
    XtSetArg(args, XmNheight, &wh);
    XtGetValues(cc, &args, 1);
    win_w = ww;
    win_h = wh;

    devwidth = win_w;
    devheight = win_h;
    wx1 = DisplayWidth(disp, DefaultScreen(disp));
    wx2 = DisplayWidthMM(disp, DefaultScreen(disp));
    wy1 = DisplayHeight(disp, DefaultScreen(disp));
    wy2 = DisplayHeightMM(disp, DefaultScreen(disp));
    devwidthmm = (int) (wx2 / wx1 * win_w);
    devheightmm = (int) (wy2 / wy1 * win_h);
    if (inc) {
	xlibinit_tiles();
	inc = 0;
    }
    if (backingstore) {
	if (!backpix[backpixind]) {
	    backpix[backpixind] = XCreatePixmap(disp, DefaultRootWindow(disp),
						win_w, win_h,
			       DisplayPlanes(disp, DefaultScreen(disp)));
	}
    }
    if (doublebuff) {
	xlibdoublebuffer(doublebuff);
	displaybuff = backbuff[0];
    } else {
	displaybuff = cwin;
    }
    if (doclear && !overlay) {
	xlibsetcolor(bgcolor);
	XFillRectangle(disp, displaybuff, gc, 0, 0, win_w, win_h);
	if (backingstore && backpix[backpixind]) {
	    XFillRectangle(disp, backpix[backpixind], gc, 0, 0, win_w, win_h);
	}
    }
}

void refresh_from_backpix()
{
    if (backingstore && backpix[backpixind]) {
	XCopyArea(disp, backpix[backpixind], cwin, gc, 0, 0, win_w, win_h, 0, 0);
    }
}

void resize_backpix()
{
    if (backingstore && backpix[backpixind]) {
	XFreePixmap(disp, backpix[backpixind]);
	backpix[backpixind] = XCreatePixmap(disp, DefaultRootWindow(disp), win_w, win_h, DisplayPlanes(disp, DefaultScreen(disp)));
    }
}

void xlibdoublebuffer(int mode)
{
    extern int inwin;

    doublebuff = mode;
    if (!inwin) {
	return;
    }
    if (mode) {
	if (!backbuff[0]) {
	    backbuff[0] = XCreatePixmap(disp,
					DefaultRootWindow(disp),
					win_w,
					win_h,
			       DisplayPlanes(disp, DefaultScreen(disp)));
	}
	displaybuff = backbuff[0];
    } else {
	if (backbuff[0]) {
	    XFreePixmap(disp, backbuff[0]);
	    backbuff[0] = NULL;
	    displaybuff = cwin;
	}
    }
}

void xlibfrontbuffer(int mode)
{
    extern int inwin;

    if (!inwin) {
	return;
    }
    if (mode) {
	displaybuff = cwin;
    } else {
	if (doublebuff && backbuff[0]) {
	    displaybuff = backbuff[0];
	}
    }
}

void xlibbackbuffer(int mode)
{
    extern int inwin;

    if (!inwin) {
	return;
    }
    if (mode && doublebuff && backbuff[0]) {
	displaybuff = backbuff[0];
    } else {
	displaybuff = cwin;
    }
}

void xlibswapbuffer()
{
    extern int inwin;

    if (!inwin) {
	return;
    }
    if (doublebuff && backbuff[0]) {
	XCopyArea(disp, displaybuff, cwin, gc, 0, 0, win_w, win_h, 0, 0);
	if (save_fli) {
	}
    }
}

static int xlib_write_mode = 1;

void set_write_mode(int m)
{
    flush_pending();
    xlib_write_mode = m;
}

void xlibsetmode(int mode)
{
    char name[128];
    char com[128];

    switch (mode) {
    case 1:
	xlibinit();
	break;
    case 2:
	flush_pending();
	if (doublebuff && backbuff[0]) {
	    xlibswapbuffer();
	}
	if (save_images) {
	    save_image_on_disk(disp, cwin, displaybuff, 0, 0, win_w, win_h, "tmp.xwd", NULL);
	    switch (save_images_type) {
	    case 0:
		sprintf(name, "%s%05d.tga", save_images_fname, save_images_count);
		save_images_count++;
		sprintf(com, "xwdtopnm < %s | ppmtotga > %s", "tmp.xwd", name);
		break;
	    case 1:
		sprintf(name, "%s%05d.xwd", save_images_fname, save_images_count);
		save_images_count++;
		sprintf(com, "/bin/mv %s %s", "tmp.xwd", name);
		break;
	    case 2:
		sprintf(name, "%s%05d.rle", save_images_fname, save_images_count);
		save_images_count++;
		sprintf(com, "xwdtopnm < %s | ppmtorle > %s", "tmp.xwd", name);
		break;
	    case 3:
		save_images_count++;
		break;
	    }
	    system(com);
	    unlink("tmp.xwd");
	}
	XFlush(disp);
	break;
    }
}

/*
 * fix for dotted/dashed linestyles
 */
#define MAXL 1024
static int npending;
XPoint polypoints[MAXL];

static unsigned char solid[1] =
{1};
static unsigned char dotted[2] =
{3, 1};
static unsigned char shortdashed[2] =
{3, 3};
static unsigned char longdashed[2] =
{7, 7};
static unsigned char dot_dashed[4] =
{1, 3, 7, 3};

static unsigned char *dash_list[] =
{
    solid,
    dotted,
    shortdashed,
    longdashed,
    dot_dashed
};

static int dash_list_length[] =
{1, 2, 2, 2, 4};

void flush_pending()
{
    int i;

    if (npending > 1) {
	if (xlib_write_mode) {
	    XDrawLines(disp, displaybuff, gc, polypoints, npending, CoordModeOrigin);
	    if (backingstore && backpix[backpixind]) {
		XDrawLines(disp, backpix[backpixind], gc, polypoints, npending, CoordModeOrigin);
	    }
	} else {
	    XDrawLines(disp, displaybuff, gcclr, polypoints, npending, CoordModeOrigin);
	    if (backingstore && backpix[backpixind]) {
		XDrawLines(disp, backpix[backpixind], gcclr, polypoints, npending, CoordModeOrigin);
	    }
	}
    }
    npending = 0;
}

static int x1, y1;

void drawxlib(int x, int y, int mode)
{
    if (mode) {
	polypoints[npending].x = x;
	polypoints[npending].y = win_h - y;
	npending++;
	if (npending == MAXL) {
	    flush_pending();
	    polypoints[npending].x = x;
	    polypoints[npending].y = win_h - y;
	    npending = 1;
	}
    } else {
	if ((x == x1 && y == y1)) {
	    return;
	} else {
	    flush_pending();
	    polypoints[npending].x = x;
	    polypoints[npending].y = win_h - y;
	    npending = 1;
	}
    }
    x1 = x;
    y1 = y;
}

int xconvxlib(double x)
{
    return ((int) (win_w * xconv(x)));
}

int yconvxlib(double y)
{
    return ((int) (win_h * yconv(y)));
}

void xlibsetfont(int n)
{
    flush_pending();
    x1 = y1 = 99999;
    hselectfont(xlibfont = n);
}

#define NUM_COLORS 256

unsigned char red[NUM_COLORS], green[NUM_COLORS], blue[NUM_COLORS];
unsigned long colors[NUM_COLORS];
XColor xc[NUM_COLORS];
int ncolors;
Colormap cmap, mycmap;

/*
 * initialize_cms_data()
 *    Initialize the colormap segment data and setup the RGB values.
 */
void initialize_cms_data2()
{
    int i, j, pm[1];
    /* white  */
    red[0] = 255;
    green[0] = 255;
    blue[0] = 255;
    /* black    */
    red[1] = 0;
    green[1] = 0;
    blue[1] = 0;
    /* red    */
    red[2] = 255;
    green[2] = 0;
    blue[2] = 0;
    /* green  */
    red[3] = 0;
    green[3] = 255;
    blue[3] = 0;
    /* blue   */
    red[4] = 0;
    green[4] = 0;
    blue[4] = 255;
    /* yellow */
    red[5] = 255;
    green[5] = 255;
    blue[5] = 0;
    /* brown  */
    red[6] = 188;
    green[6] = 143;
    blue[6] = 143;
    /* gray   */
    red[7] = 220;
    green[7] = 220;
    blue[7] = 220;
    /* violet  */
    red[8] = 148;
    green[8] = 0;
    blue[8] = 211;
    /* cyan  */
    red[9] = 0;
    green[9] = 255;
    blue[9] = 255;
    /* magenta  */
    red[10] = 255;
    green[10] = 0;
    blue[10] = 211;
    /* orange  */
    red[11] = 255;
    green[11] = 158;
    blue[11] = 0;
    /* blue violet  */
    red[12] = 114;
    green[12] = 33;
    blue[12] = 188;
    /* maroon  */
    red[13] = 103;
    green[13] = 7;
    blue[13] = 72;
    /* turquoise  */
    red[14] = 72;
    green[14] = 209;
    blue[14] = 204;
    /* forest green  */
    red[15] = 85;
    green[15] = 192;
    blue[15] = 52;
    for (i = 31; i > 15; i--) {
	j = (i - 15);
	red[i] = (j - 32) * 8;
	green[i] = (j - 32) * 8;
	blue[i] = 255;
    }
    for (i = 47; i > 31; i--) {
	j = (i - 31);
	blue[i] = (j - 48) * 8;
	green[i] = (j - 48) * 8;
	red[i] = 255;
    }
}

void xlibinitPseudoColormap()
{
    int i;
    unsigned long plane_masks[1];

    ncolors = DisplayCells(disp, DefaultScreen(disp));
    if (ncolors > 256) {
	ncolors = 256;
    }
    if (ncolors > 16) {
	cmap = DefaultColormap(disp, DefaultScreen(disp));
	if (XAllocColorCells(disp, cmap, False, plane_masks, 0, colors, 48)) {
	    mycmap = cmap;
	}
	for (i = 0; i < MAXCOLOR; i++) {
	    xlibsetcmap(i, red[i], green[i], blue[i]);
	}
	if (revflag) {
	    iswap(&colors[0], &colors[1]);
	    cswap(&red[0], &red[1]);
	    cswap(&green[0], &green[1]);
	    cswap(&blue[0], &blue[1]);
	}
    } else {
	if (revflag) {
	    colors[1] = WhitePixel(disp, DefaultScreen(disp));
	    colors[0] = BlackPixel(disp, DefaultScreen(disp));
	} else {
	    colors[0] = WhitePixel(disp, DefaultScreen(disp));
	    colors[1] = BlackPixel(disp, DefaultScreen(disp));
	}
    }
}

void xlibinitTrueColormap(void)
{
    int i;

    ncolors = DisplayCells(disp, DefaultScreen(disp));
    if (ncolors > 256) {
	ncolors = 256;
    }
    if (ncolors > 16) {
	cmap = DefaultColormap(disp, DefaultScreen(disp));
	for (i = 0; i < ncolors; i++) {
	    xc[i].pixel = i;
	    xc[i].flags = DoRed | DoGreen | DoBlue;
	}
	mycmap = cmap;
	for (i = 2; i < MAXCOLOR; i++) {
	    xc[i].red = red[i] << 8;
	    xc[i].green = green[i] << 8;
	    xc[i].blue = blue[i] << 8;
	    if (!XAllocColor(disp, cmap, &xc[i])) {
		fprintf(stderr, " Can't allocate color\n");
	    }
	    colors[i] = xc[i].pixel;
	}
    }
    if (revflag) {
	colors[1] = WhitePixel(disp, DefaultScreen(disp));
	colors[0] = BlackPixel(disp, DefaultScreen(disp));
    } else {
	colors[0] = WhitePixel(disp, DefaultScreen(disp));
	colors[1] = BlackPixel(disp, DefaultScreen(disp));
    }
}

void xlibinitcmap()
{
    int depth;
    depth = DisplayPlanes(disp, DefaultScreen(disp));
    if (depth == 8) {
	xlibinitPseudoColormap();
    } else {
	xlibinitTrueColormap();
    }
}

void xlibsetcmap(int i, int r, int g, int b)
{
    XColor xct;
    int depth;
    if (inwin) {
        depth = DisplayPlanes(disp, DefaultScreen(disp));
    } else {
	depth = 9;
    }
    if (depth > 8) {
	red[i] = r;
	green[i] = g;
	blue[i] = b;
	cmscolors[i].red = red[i];
	cmscolors[i].green = green[i];
	cmscolors[i].blue = blue[i];
	xct.red = red[i] << 8;
	xct.green = green[i] << 8;
	xct.blue = blue[i] << 8;
	xct.flags = DoRed | DoGreen | DoBlue;
	if (inwin && use_colors > 2 && i >= 0) {
	    if (!XAllocColor(disp, cmap, &xct)) {
		fprintf(stderr, " Can't allocate color\n");
	    }
	    colors[i] = xct.pixel;
	}
	return;
    }
    red[i] = r;
    green[i] = g;
    blue[i] = b;
    cmscolors[i].red = red[i];
    cmscolors[i].green = green[i];
    cmscolors[i].blue = blue[i];
    if (inwin && use_colors > 2 && i >= 0) {
	xct.green = g << 8;
	xct.blue = b << 8;
	xct.red = r << 8;
	xct.flags = DoRed | DoGreen | DoBlue;
	xct.pixel = colors[i];
	xct.pad = 0;

	XStoreColor(disp, cmap, &xct);
    }
}

int xlibsetlinewidth(int c)
{
    flush_pending();
    x1 = y1 = 99999;
    if (c) {
	c = c % MAXLINEW;
	if (c == 0)
	    c = 1;
	if (xliblinestyle == 1) {
	    XSetLineAttributes(disp, gc, c - 1 == 0 ? 0 : c, LineSolid, CapButt, JoinRound);
	} else {
	    XSetLineAttributes(disp, gc, c - 1 == 0 ? 0 : c, LineOnOffDash, CapButt, JoinRound);
	}
    }
    return (xliblinewidth = c);
}

int xlibsetlinestyle(int style)
{
    flush_pending();
    x1 = y1 = 99999;
    if (style > 1 && xliblinewidth) {
	XSetLineAttributes(disp, gc, xliblinewidth - 1 ? 0 : xliblinewidth, LineOnOffDash, CapButt, JoinRound);
	XSetDashes(disp, gc, 0, dash_list[style - 1], dash_list_length[style - 1]);
    } else if (style == 1 && xliblinewidth) {
	XSetLineAttributes(disp, gc, xliblinewidth - 1 ? 0 : xliblinewidth, LineSolid, CapButt, JoinRound);
    }
    return (xliblinestyle = style);
}

int xlibsetcolor(int c)
{
    flush_pending();
    x1 = y1 = 99999;
    c = c % MAXCOLOR;
    if (use_colors > 2) {
	XSetForeground(disp, gc, colors[c]);
    } else {
	XSetForeground(disp, gc, colors[c != 0]);
    }
    xlibcolor = c;
    return c;
}

void xlibdrawtic(int x, int y, int dir, int updown)
{
    switch (dir) {
    case 0:
	switch (updown) {
	case 0:
	    drawxlib(x, y, 0);
	    drawxlib(x, y + devxticl, 1);
	    break;
	case 1:
	    drawxlib(x, y, 0);
	    drawxlib(x, y - devxticl, 1);
	    break;
	}
	break;
    case 1:
	switch (updown) {
	case 0:
	    drawxlib(x, y, 0);
	    drawxlib(x + devyticl, y, 1);
	    break;
	case 1:
	    drawxlib(x, y, 0);
	    drawxlib(x - devyticl, y, 1);
	    break;
	}
	break;
    }
}

void dispstrxlib(int x, int y, int rot, char *s, int just, int fudge)
{
    flush_pending();
    x1 = y1 = 99999;
    puthersh(x, y, xlibcharsize * charsize, rot, just, xlibcolor, drawxlib, s);
    flush_pending();
    x1 = y1 = 99999;
}

#define MAXPATTERNS 16

static int patno = 0;

static Pixmap tiles[30];
static Pixmap curtile;

void xlibinit_tiles()
{
    int i;
    Pixmap ptmp;

    for (i = 0; i < MAXPATTERNS; i++) {
	tiles[i] = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    }
    for (i = 0; i < MAXPATTERNS; i++) {
	if (tiles[i] == NULL) {
	    printf("bad tile %d\n", i);
	} else {
	    XFillRectangle(disp, tiles[i], gcclr, 0, 0, 16, 16);
	}
    }
    ptmp = XCreateBitmapFromData(disp, cwin, pat0_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[0], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat1_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[1], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat2_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[2], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat3_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[3], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat4_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[4], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat5_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[5], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat6_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[6], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat7_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[7], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat8_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[8], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat9_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[9], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat10_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[10], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat11_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[11], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat12_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[12], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat13_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[13], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat14_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[14], gc, 0, 0, 16, 16, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, pat15_bits, 16, 16);
    XCopyPlane(disp, ptmp, tiles[15], gc, 0, 0, 16, 16, 0, 0, 1);
    curtile = tiles[0];
}

int xlibsetpat(int k)
{
    patno = k;
    if (k > MAXPATTERNS) {
	k = 1;
    }
    if (patno != 0) {
	curtile = tiles[k - 1];
    }
}

void xlibfill(int n, int *px, int *py)
{
    int i, x, y;
    XPoint *p;

    p = (XPoint *) calloc(n, sizeof(XPoint));
    if (p == NULL) {
	return;
    }
    if (patno == 0) {
	return;
    }
    XSetFillStyle(disp, gc, FillTiled);
    XSetTile(disp, gc, curtile);
    for (i = 0; i < n; i++) {
	p[i].x = px[i];
	p[i].y = win_h - py[i];
    }
    XFillPolygon(disp, displaybuff, gc, p, n, Nonconvex, CoordModeOrigin);
    if (backingstore && backpix[backpixind]) {
	XFillPolygon(disp, backpix[backpixind], gc, p, n, Nonconvex, CoordModeOrigin);
    }
    XSetFillStyle(disp, gc, FillSolid);
    free(p);
}

void xlibfillcolor(int n, int *px, int *py)
{
    int i, x, y;
    XPoint *p;

    p = (XPoint *) calloc(n, sizeof(XPoint));
    if (p == NULL) {
	return;
    }
    for (i = 0; i < n; i++) {
	p[i].x = px[i];
	p[i].y = win_h - py[i];
    }
    XFillPolygon(disp, displaybuff, gc, p, n, Nonconvex, CoordModeOrigin);
    if (backingstore && backpix[backpixind]) {
	XFillPolygon(disp, backpix[backpixind], gc, p, n, Nonconvex, CoordModeOrigin);
    }
    free(p);
}

void xlibdrawarc(int x, int y, int r)
{
    XDrawArc(disp, displaybuff, gc, x - r, win_h - (y + r), 2 * r, 2 * r, 0, 360 * 64);
    if (backingstore && backpix[backpixind]) {
	XDrawArc(disp, backpix[backpixind], gc, x - r, win_h - (y + r), 2 * r, 2 * r, 0, 360 * 64);
    }
}

void xlibfillarc(int x, int y, int r)
{
    XFillArc(disp, displaybuff, gc, x - r, win_h - (y + r), 2 * r, 2 * r, 0, 360 * 64);
    if (backingstore && backpix[backpixind]) {
	XFillArc(disp, backpix[backpixind], gc, x - r, win_h - (y + r), 2 * r, 2 * r, 0, 360 * 64);
    }
}

void xlibdrawellipse(int x, int y, int xm, int ym)
{
    XDrawArc(disp, displaybuff, gc, x - xm, win_h - (y + ym), 2 * xm, 2 * ym, 0, 360 * 64);
    if (backingstore && backpix[backpixind]) {
	XDrawArc(disp, backpix[backpixind], gc, x - xm, win_h - (y + ym), 2 * xm, 2 * ym, 0, 360 * 64);
    }
}

void xlibfillellipse(int x, int y, int xm, int ym)
{
    XFillArc(disp, displaybuff, gc, x - xm, win_h - (y + ym), 2 * xm, 2 * ym, 0, 360 * 64);
    if (backingstore && backpix[backpixind]) {
	XFillArc(disp, backpix[backpixind], gc, x - xm, win_h - (y + ym), 2 * xm, 2 * ym, 0, 360 * 64);
    }
}

void xlibleavegraphics()
{
    flush_pending();
    x1 = y1 = 99999;
    xlibsetmode(2);
    XFlush(disp);
}

int xlibinitgraphics(int dmode)
{
    npending = 0;
    x1 = 99999;
    y1 = 99999;
    xlibdmode = dmode;
    xlibsetmode(1);
    devorient = 1;
    devoffsx = 0;
    devoffsy = 0;
    devconvx = xconvxlib;
    devconvy = yconvxlib;
    vector = drawxlib;
    devwritestr = dispstrxlib;
    devsetcolor = xlibsetcolor;
    devsetfont = xlibsetfont;
    devsetline = xlibsetlinestyle;
    devsetlinew = xlibsetlinewidth;
    devdrawtic = xlibdrawtic;
    devsetpat = xlibsetpat;
    devdrawarc = xlibdrawarc;
    devfillarc = xlibfillarc;
    devdrawellipse = xlibdrawellipse;
    devfillellipse = xlibfillellipse;
    devfill = xlibfill;
    devfillcolor = xlibfillcolor;
    devleavegraphics = xlibleavegraphics;
    devxticl = 12;
    devyticl = 12;
    devarrowlength = 12;
    devsymsize = 6;
    devcharsize = xlibcharsize;
    (*devsetcolor) (1);
    xlibsetlinestyle(1);
    return 0;
}

#define ddashed_width 16
#define ddashed_height 7
static char ddashed_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe3, 0x33, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define dotted_width 16
#define dotted_height 7
static char dotted_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x33, 0x33, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define ldashed_width 16
#define ldashed_height 7
static char ldashed_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0x7e, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define sdashed_width 16
#define sdashed_height 7
static char sdashed_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc7, 0x71, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define solid_width 16
#define solid_height 7
static char solid_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define none_width 16
#define none_height 7
static char none_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define w0_width 16
#define w0_height 7
static char w0_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define w1_width 16
#define w1_height 7
static char w1_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define w2_width 16
#define w2_height 7
static char w2_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00};

#define w3_width 16
#define w3_height 7
static char w3_bits[] =
{
  0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00,
    0x00, 0x00};

#define w4_width 16
#define w4_height 7
static char w4_bits[] =
{
  0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00,
    0x00, 0x00};

#define w5_width 16
#define w5_height 7
static char w5_bits[] =
{
  0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x00, 0x00};

#define w6_width 16
#define w6_height 7
static char w6_bits[] =
{
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x00, 0x00};

#define w7_width 16
#define w7_height 7
static char w7_bits[] =
{
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff};

Pixmap linespm[6];
Pixmap linewpm[10];

void xlibinit_pm()
{
    int i;
    Pixmap ptmp;
    for (i = 0; i < 6; i++) {
	linespm[i] = XCreatePixmap(disp, cwin, 16, 7, DisplayPlanes(disp, DefaultScreen(disp)));
    }
    for (i = 0; i < 8; i++) {
	linewpm[i] = XCreatePixmap(disp, cwin, 16, 7, DisplayPlanes(disp, DefaultScreen(disp)));
    }
    ptmp = XCreateBitmapFromData(disp, cwin, none_bits, 16, 7);
    XCopyPlane(disp, ptmp, linespm[0], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, solid_bits, 16, 7);
    XCopyPlane(disp, ptmp, linespm[1], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, dotted_bits, 16, 7);
    XCopyPlane(disp, ptmp, linespm[2], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, sdashed_bits, 16, 7);
    XCopyPlane(disp, ptmp, linespm[3], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, ddashed_bits, 16, 7);
    XCopyPlane(disp, ptmp, linespm[4], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, ldashed_bits, 16, 7);
    XCopyPlane(disp, ptmp, linespm[5], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w0_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[0], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w1_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[1], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w2_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[2], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w3_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[3], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w4_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[4], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w5_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[5], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w6_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[6], gc, 0, 0, 16, 7, 0, 0, 1);
    ptmp = XCreateBitmapFromData(disp, cwin, w7_bits, 16, 7);
    XCopyPlane(disp, ptmp, linewpm[7], gc, 0, 0, 16, 7, 0, 0, 1);
}

/*
 * cursors
 */

#define MOTIF

#include <X11/cursorfont.h>

static Cursor wait_cursor;
static Cursor line_cursor;
static Cursor find_cursor;
static Cursor move_cursor;
static Cursor text_cursor;
static Cursor kill_cursor;
static int cur_cursor = -1;
void set_cursor();

void set_wait_cursor(Widget w)
{
    XDefineCursor(disp, cwin, wait_cursor);
    if (w != NULL) {
	XDefineCursor(disp, cwin, wait_cursor);
    }
    XFlush(disp);
}

void unset_wait_cursor(Widget w)
{
    if (cur_cursor == -1) {
	XUndefineCursor(disp, cwin);
	if (w != NULL) {
	    XUndefineCursor(disp, cwin);
	}
    } else {
	set_cursor(cur_cursor);
    }
}

void set_cursor(int c)
{
    XUndefineCursor(disp, cwin);
    cur_cursor = -1;
    switch (c) {
    case 0:
	XDefineCursor(disp, cwin, line_cursor);
	cur_cursor = 0;
	break;
    case 1:
	XDefineCursor(disp, cwin, find_cursor);
	cur_cursor = 1;
	break;
    case 2:
	XDefineCursor(disp, cwin, text_cursor);
	cur_cursor = 2;
	break;
    case 3:
	XDefineCursor(disp, cwin, kill_cursor);
	cur_cursor = 3;
	break;
    case 4:
	XDefineCursor(disp, cwin, move_cursor);
	cur_cursor = 4;
	break;
    }
}

void init_cursors()
{
    wait_cursor = XCreateFontCursor(disp, XC_watch);
    line_cursor = XCreateFontCursor(disp, XC_crosshair);
    find_cursor = XCreateFontCursor(disp, XC_hand2);
    text_cursor = XCreateFontCursor(disp, XC_xterm);
    kill_cursor = XCreateFontCursor(disp, XC_pirate);
    move_cursor = XCreateFontCursor(disp, XC_fleur);
    cur_cursor = -1;
}

void clearscreen(int c)
{
    setcolor(c);
    XFillRectangle(disp, displaybuff, gc, 0, 0, win_w, win_h);
}
