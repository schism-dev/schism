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
 *
 * Overlay one or more ESRI shapefiles (.shp) on top of the grid.
 *
 * This is a self-contained reader for the ESRI shapefile main file
 * (.shp). It understands the 2D geometry of Point, MultiPoint,
 * PolyLine and Polygon records, including their Z and M variants
 * (only the X,Y portion is used). No external library is required.
 *
 * The shapefile is assumed to be in the SAME coordinate system as the
 * grid (e.g. both lon/lat, or both the same projected meters); no
 * re-projection is performed. Geometry is drawn in world coordinates
 * so it pans and zooms together with the mesh.
 *
 * Added 2026.
 *
 */
#ifndef lint
static char RCSid[] = "$Id: shapefile.c,v 1.0 2026/06/05 00:00:00 ace Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "graphics.h"
#include "symdefs.h"

extern Widget app_shell;
extern Widget canvas;
extern XmStringCharSet charset;

/* ESRI shapefile shape type codes */
#define SHP_NULL         0
#define SHP_POINT        1
#define SHP_POLYLINE     3
#define SHP_POLYGON      5
#define SHP_MULTIPOINT   8
#define SHP_POINTZ      11
#define SHP_POLYLINEZ   13
#define SHP_POLYGONZ    15
#define SHP_MULTIPOINTZ 18
#define SHP_POINTM      21
#define SHP_POLYLINEM   23
#define SHP_POLYGONM    25
#define SHP_MULTIPOINTM 28

#define MAX_SHAPEFILES  50

/* one polyline / ring, or a bag of points (for point type shapefiles) */
typedef struct {
    int npts;
    int maxpts;
    double *x;
    double *y;
} ShpPart;

typedef struct {
    char fname[1024];
    int loaded;
    int ispoint;		/* 1 -> draw as point symbols, 0 -> polylines */
    int shapetype;
    int nparts;
    int maxparts;
    ShpPart *parts;
    double xmin, ymin, xmax, ymax;
} Shapefile;

static Shapefile shapes[MAX_SHAPEFILES];
static int num_shapes = 0;

/* display settings */
int drawshape_flag = 1;
static int shape_line_color = 1;	/* black */
static int shape_point_color = 2;	/* red */
static int shape_linewidth = 1;
static double shape_point_size = 0.8;

/* ----------------------------------------------------------------- */
/* low level: portable readers for big/little endian fields          */
/* ----------------------------------------------------------------- */

static int host_is_little(void)
{
    unsigned int x = 1;
    return *(char *) &x == 1;
}

static int read_be_i32(FILE * fp, int *ok)
{
    unsigned char b[4];
    if (fread(b, 1, 4, fp) != 4) {
	*ok = 0;
	return 0;
    }
    return (int) (((unsigned) b[0] << 24) | ((unsigned) b[1] << 16) |
		  ((unsigned) b[2] << 8) | (unsigned) b[3]);
}

static int read_le_i32(FILE * fp, int *ok)
{
    unsigned char b[4];
    if (fread(b, 1, 4, fp) != 4) {
	*ok = 0;
	return 0;
    }
    return (int) ((unsigned) b[0] | ((unsigned) b[1] << 8) |
		  ((unsigned) b[2] << 16) | ((unsigned) b[3] << 24));
}

static double read_le_double(FILE * fp, int *ok)
{
    unsigned char b[8];
    union {
	double d;
	unsigned char c[8];
    } u;
    int i;
    if (fread(b, 1, 8, fp) != 8) {
	*ok = 0;
	return 0.0;
    }
    if (host_is_little()) {
	for (i = 0; i < 8; i++)
	    u.c[i] = b[i];
    } else {
	for (i = 0; i < 8; i++)
	    u.c[i] = b[7 - i];
    }
    return u.d;
}

/* ----------------------------------------------------------------- */
/* in-memory geometry construction                                   */
/* ----------------------------------------------------------------- */

/* append a new (empty) part to a shapefile, return its index or -1 */
static int shp_new_part(Shapefile * s, int wantpts)
{
    ShpPart *p;
    if (s->nparts >= s->maxparts) {
	int n = s->maxparts ? s->maxparts * 2 : 16;
	ShpPart *np = (ShpPart *) realloc(s->parts, n * sizeof(ShpPart));
	if (np == NULL)
	    return -1;
	s->parts = np;
	s->maxparts = n;
    }
    p = &s->parts[s->nparts];
    p->npts = 0;
    p->maxpts = wantpts > 0 ? wantpts : 16;
    p->x = (double *) malloc(p->maxpts * sizeof(double));
    p->y = (double *) malloc(p->maxpts * sizeof(double));
    if (p->x == NULL || p->y == NULL)
	return -1;
    return s->nparts++;
}

static int shp_add_point(ShpPart * p, double x, double y)
{
    if (p->npts >= p->maxpts) {
	int n = p->maxpts ? p->maxpts * 2 : 16;
	double *nx = (double *) realloc(p->x, n * sizeof(double));
	double *ny = (double *) realloc(p->y, n * sizeof(double));
	if (nx == NULL || ny == NULL)
	    return 0;
	p->x = nx;
	p->y = ny;
	p->maxpts = n;
    }
    p->x[p->npts] = x;
    p->y[p->npts] = y;
    p->npts++;
    return 1;
}

static void shp_free(Shapefile * s)
{
    int i;
    for (i = 0; i < s->nparts; i++) {
	if (s->parts[i].x)
	    free(s->parts[i].x);
	if (s->parts[i].y)
	    free(s->parts[i].y);
    }
    if (s->parts)
	free(s->parts);
    s->parts = NULL;
    s->nparts = 0;
    s->maxparts = 0;
    s->loaded = 0;
}

void clear_shapefiles(void)
{
    int i;
    for (i = 0; i < num_shapes; i++)
	shp_free(&shapes[i]);
    num_shapes = 0;
}

/* ----------------------------------------------------------------- */
/* the reader                                                        */
/* ----------------------------------------------------------------- */

/*
 * Returns 1 on success, 0 on failure. On success a new entry is added
 * to shapes[].
 */
int read_shapefile(char *fname)
{
    FILE *fp;
    Shapefile *s;
    int ok = 1;
    int filecode, filelen, version, ftype;
    int i, recnum, contentlen, rectype;
    long contentstart;
    double bx[4];
    int pointtype;

    if (num_shapes >= MAX_SHAPEFILES) {
	errwin("Too many shapefiles loaded; clear some first");
	return 0;
    }
    fp = fopen(fname, "rb");
    if (fp == NULL) {
	errwin("Unable to open shapefile");
	return 0;
    }

    /* ---- 100 byte file header ---- */
    filecode = read_be_i32(fp, &ok);	/* bytes 0-3 */
    for (i = 0; i < 5; i++)
	(void) read_be_i32(fp, &ok);	/* bytes 4-23 unused */
    filelen = read_be_i32(fp, &ok);	/* bytes 24-27, 16-bit words */
    version = read_le_i32(fp, &ok);	/* bytes 28-31 */
    ftype = read_le_i32(fp, &ok);	/* bytes 32-35, overall type */
    for (i = 0; i < 4; i++)
	bx[i] = read_le_double(fp, &ok);	/* bbox xmin,ymin,xmax,ymax */
    for (i = 0; i < 4; i++)
	(void) read_le_double(fp, &ok);	/* zmin,zmax,mmin,mmax */

    if (!ok || filecode != 9994) {
	fclose(fp);
	errwin("Not a valid ESRI shapefile (.shp)");
	return 0;
    }

    s = &shapes[num_shapes];
    memset(s, 0, sizeof(Shapefile));
    strncpy(s->fname, fname, sizeof(s->fname) - 1);
    s->shapetype = ftype;
    s->xmin = bx[0];
    s->ymin = bx[1];
    s->xmax = bx[2];
    s->ymax = bx[3];

    pointtype = (ftype == SHP_POINT || ftype == SHP_POINTZ ||
		 ftype == SHP_POINTM || ftype == SHP_MULTIPOINT ||
		 ftype == SHP_MULTIPOINTZ || ftype == SHP_MULTIPOINTM);
    s->ispoint = pointtype;

    /* for point type shapefiles, collect everything into one part */
    if (pointtype) {
	if (shp_new_part(s, 256) < 0) {
	    shp_free(s);
	    fclose(fp);
	    errwin("Out of memory reading shapefile");
	    return 0;
	}
    }

    /* ---- records ---- */
    while (1) {
	ok = 1;
	recnum = read_be_i32(fp, &ok);
	contentlen = read_be_i32(fp, &ok);
	if (!ok)
	    break;		/* normal end of file */
	contentstart = ftell(fp);
	rectype = read_le_i32(fp, &ok);
	if (!ok)
	    break;

	switch (rectype) {
	case SHP_NULL:
	    break;
	case SHP_POINT:
	case SHP_POINTZ:
	case SHP_POINTM:
	    {
		double x = read_le_double(fp, &ok);
		double y = read_le_double(fp, &ok);
		if (ok)
		    shp_add_point(&s->parts[0], x, y);
	    }
	    break;
	case SHP_MULTIPOINT:
	case SHP_MULTIPOINTZ:
	case SHP_MULTIPOINTM:
	    {
		int np, k;
		for (i = 0; i < 4; i++)
		    (void) read_le_double(fp, &ok);	/* box */
		np = read_le_i32(fp, &ok);
		for (k = 0; k < np && ok; k++) {
		    double x = read_le_double(fp, &ok);
		    double y = read_le_double(fp, &ok);
		    if (ok)
			shp_add_point(&s->parts[0], x, y);
		}
	    }
	    break;
	case SHP_POLYLINE:
	case SHP_POLYGON:
	case SHP_POLYLINEZ:
	case SHP_POLYGONZ:
	case SHP_POLYLINEM:
	case SHP_POLYGONM:
	    {
		int nparts, npoints, k, pp;
		int *partidx;
		for (i = 0; i < 4; i++)
		    (void) read_le_double(fp, &ok);	/* box */
		nparts = read_le_i32(fp, &ok);
		npoints = read_le_i32(fp, &ok);
		if (!ok || nparts <= 0 || npoints <= 0)
		    break;
		partidx = (int *) malloc((nparts + 1) * sizeof(int));
		if (partidx == NULL)
		    break;
		for (k = 0; k < nparts; k++)
		    partidx[k] = read_le_i32(fp, &ok);
		partidx[nparts] = npoints;
		for (pp = 0; pp < nparts && ok; pp++) {
		    int start = partidx[pp];
		    int end = partidx[pp + 1];
		    int idx = shp_new_part(s, end - start);
		    if (idx < 0)
			break;
		    for (k = start; k < end && ok; k++) {
			double x = read_le_double(fp, &ok);
			double y = read_le_double(fp, &ok);
			if (ok)
			    shp_add_point(&s->parts[idx], x, y);
		    }
		}
		free(partidx);
	    }
	    break;
	default:
	    /* unsupported type, skip via content length */
	    break;
	}

	/* jump to the next record regardless of trailing Z/M data */
	if (fseek(fp, contentstart + (long) contentlen * 2, SEEK_SET) != 0)
	    break;
    }

    fclose(fp);

    if (s->nparts == 0 || (pointtype && s->parts[0].npts == 0)) {
	shp_free(s);
	errwin("Shapefile contained no usable geometry");
	return 0;
    }
    s->loaded = 1;
    num_shapes++;
    return 1;
}

/* ----------------------------------------------------------------- */
/* drawing - called from do_drawobjects()                            */
/* ----------------------------------------------------------------- */

void draw_shapefile(void)
{
    int sidx, p, i;
    if (!drawshape_flag)
	return;
    for (sidx = 0; sidx < num_shapes; sidx++) {
	Shapefile *s = &shapes[sidx];
	if (!s->loaded)
	    continue;
	if (s->ispoint) {
	    setcolor(shape_point_color);
	    for (p = 0; p < s->nparts; p++) {
		if (s->parts[p].npts > 0)
		    drawpolysym(s->parts[p].x, s->parts[p].y,
				s->parts[p].npts, SYM_CIRCLE, 0, 1,
				shape_point_size);
	    }
	} else {
	    setcolor(shape_line_color);
	    setlinewidth(shape_linewidth);
	    for (p = 0; p < s->nparts; p++) {
		ShpPart *pt = &s->parts[p];
		if (pt->npts < 1)
		    continue;
		my_move2(pt->x[0], pt->y[0]);
		for (i = 1; i < pt->npts; i++)
		    my_draw2(pt->x[i], pt->y[i]);
	    }
	    setlinewidth(1);
	}
    }
    setcolor(1);
}

/* union extent of all loaded shapefiles; returns 1 if any */
static int shapefile_extent(double *x1, double *y1, double *x2, double *y2)
{
    int i, got = 0;
    for (i = 0; i < num_shapes; i++) {
	if (!shapes[i].loaded)
	    continue;
	if (!got) {
	    *x1 = shapes[i].xmin;
	    *y1 = shapes[i].ymin;
	    *x2 = shapes[i].xmax;
	    *y2 = shapes[i].ymax;
	    got = 1;
	} else {
	    if (shapes[i].xmin < *x1)
		*x1 = shapes[i].xmin;
	    if (shapes[i].ymin < *y1)
		*y1 = shapes[i].ymin;
	    if (shapes[i].xmax > *x2)
		*x2 = shapes[i].xmax;
	    if (shapes[i].ymax > *y2)
		*y2 = shapes[i].ymax;
	}
    }
    return got;
}

/* ----------------------------------------------------------------- */
/* Motif user interface                                              */
/* ----------------------------------------------------------------- */

static Widget shape_frame = NULL;
static Widget shape_display_item;
static Widget shape_linecol_item;
static Widget shape_pointcol_item;
static Widget shape_linewidth_item;
static Widget shape_status_item;
static Widget rshape_dialog = NULL;

void create_rshape_popup(Widget w, XtPointer client_data, XtPointer call_data);

static void update_shape(void)
{
    char buf[256];
    if (shape_frame == NULL)
	return;
    XmToggleButtonSetState(shape_display_item, drawshape_flag, False);
    sprintf(buf, "%d", shape_line_color);
    xv_setstr(shape_linecol_item, buf);
    sprintf(buf, "%d", shape_point_color);
    xv_setstr(shape_pointcol_item, buf);
    sprintf(buf, "%d", shape_linewidth);
    xv_setstr(shape_linewidth_item, buf);
    sprintf(buf, "%d shapefile(s) loaded", num_shapes);
    if (shape_status_item)
	xv_setstr(shape_status_item, buf);
}

void do_accept_shape_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int v;
    drawshape_flag = XmToggleButtonGetState(shape_display_item);
    shape_line_color = atoi((char *) xv_getstr(shape_linecol_item));
    shape_point_color = atoi((char *) xv_getstr(shape_pointcol_item));
    v = atoi((char *) xv_getstr(shape_linewidth_item));
    if (v < 1)
	v = 1;
    shape_linewidth = v;
    do_drawgrid();
}

void do_clear_shape_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    clear_shapefiles();
    update_shape();
    do_drawgrid();
}

void do_fit_shape_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    double x1, y1, x2, y2, dx, dy;
    if (!shapefile_extent(&x1, &y1, &x2, &y2)) {
	errwin("No shapefiles loaded");
	return;
    }
    dx = (x2 - x1) * 0.05;
    dy = (y2 - y1) * 0.05;
    if (dx == 0.0)
	dx = 1.0;
    if (dy == 0.0)
	dy = 1.0;
    xg1 = x1 - dx;
    xg2 = x2 + dx;
    yg1 = y1 - dy;
    yg2 = y2 + dy;
    do_drawgrid();
}

void do_rshape_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char *s;
    char fname[1024];
    XmFileSelectionBoxCallbackStruct *cbs =
    (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in do_rshape_proc()");
	return;
    }
    strcpy(fname, s);
    XtFree(s);
    set_wait_cursor();
    read_shapefile(fname);
    unset_wait_cursor();
    XtUnmanageChild(rshape_dialog);
    update_shape();
    do_drawgrid();
}

void close_rshape_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    XtUnmanageChild(rshape_dialog);
}

void create_rshape_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    if (rshape_dialog == NULL) {
	rshape_dialog =
	    XmCreateFileSelectionDialog(app_shell, "rshape_dialog", NULL, 0);
	XtVaSetValues(XtParent(rshape_dialog),
		      XmNtitle, "Read shapefile (.shp)", NULL);
	XtVaSetValues(rshape_dialog,
		      XmNdirMask, XmStringCreate("*.shp", charset), NULL);
	XtAddCallback(rshape_dialog, XmNcancelCallback,
		      (XtCallbackProc) close_rshape_popup, (XtPointer) NULL);
	XtAddCallback(rshape_dialog, XmNokCallback,
		      (XtCallbackProc) do_rshape_proc, (XtPointer) NULL);
    }
    XtRaise(rshape_dialog);
}

void create_shape_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget dialog;
    Widget but4[4];

    set_wait_cursor();
    if (shape_frame == NULL) {
	char *label4[4];
	label4[0] = "Accept";
	label4[1] = "Read shapefile...";
	label4[2] = "Clear all";
	label4[3] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	shape_frame = XmCreateDialogShell(app_shell, "Shapefiles", NULL, 0);
	handle_close(shape_frame);
	XtVaSetValues(shape_frame, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(shape_frame, "dialog_rc", NULL, 0);

	shape_status_item = CreateTextItem2(dialog, 30, "Status: ");
	shape_display_item =
	    XtVaCreateManagedWidget("Display shapefiles",
				    xmToggleButtonWidgetClass, dialog, NULL);
	shape_linecol_item =
	    CreateTextItem2(dialog, 6, "Line color (0-15): ");
	shape_pointcol_item =
	    CreateTextItem2(dialog, 6, "Point color (0-15): ");
	shape_linewidth_item =
	    CreateTextItem2(dialog, 6, "Line width: ");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	{
	    Widget fitbut =
		XtVaCreateManagedWidget("Fit view to shapefiles",
					xmPushButtonWidgetClass, dialog, NULL);
	    XtAddCallback(fitbut, XmNactivateCallback,
			  (XtCallbackProc) do_fit_shape_proc, (XtPointer) NULL);
	}

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 4, but4, label4);
	XtAddCallback(but4[0], XmNactivateCallback,
		      (XtCallbackProc) do_accept_shape_proc, (XtPointer) NULL);
	XtAddCallback(but4[1], XmNactivateCallback,
		      (XtCallbackProc) create_rshape_popup, (XtPointer) NULL);
	XtAddCallback(but4[2], XmNactivateCallback,
		      (XtCallbackProc) do_clear_shape_proc, (XtPointer) NULL);
	XtAddCallback(but4[3], XmNactivateCallback,
		      (XtCallbackProc) destroy_dialog, (XtPointer) shape_frame);

	XtManageChild(dialog);
    }
    XtRaise(shape_frame);
    update_shape();
    unset_wait_cursor();
}
