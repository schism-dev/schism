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

 * canvas event proc and set_action()
 *
 */

#ifndef lint
static char RCSid[] = "$Id: events.c,v 1.6 2004/08/05 21:32:05 pturner Exp $";

#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

void my_proc(Widget w, XtPointer data, XEvent * event);
void draw_focus(int gno);
int stringextentx(double size, char *s);

Widget mode_item;
Widget locate_grid_item;
Widget calc_item;
Widget locate_item;		/* locator on main_panel */
Widget comparea_item;
Widget compperi_item;

int go_locateflag = TRUE;

FILE *maxelevfp = NULL;

extern Display *disp;
extern Window cwin;
extern GC gc;
extern GC gcxor;
extern GC gcclr;

extern XmStringCharSet charset;

extern int curplaceitem;

extern int spreadflag;
extern int cercnum;

double area(int gridno, int elem), area_from_nodes(int gridno, int n1, int n2, int n3), get_depth_element(int gno, int ind, double x, double y);

double xconv(), yconv();
double setcharsize();
void motion(XMotionEvent * e);

extern double charsize, xlibcharsize;	/* declared in draw.c and xlib.c
					 * resp. */

static int action_flag = 0;
static int xs, ys;
static int mbox_flag = 0;	/* moving box attached to cursor */
static int mline_flag = 0;	/* moving line attached to cursor */

static int rectflag = 0;
static int rubber_flag = 0;

static double wx1, wx2, wy1, wy2, wx3, wy3, wx, wy;
static int sx, sy, old_x, old_y;
static int tmpind;
static int bnd1, bnd2, ind1, ind2, ind3;	/* for adding elements */
static double sm, sb;		/* for snap lines */
static int bnof, indf;
static int n0, n1, n2, n3;

static int ipolyx[MAXSCRATCH];
static int ipolyy[MAXSCRATCH];
static int nipoly;
static double polyx[MAXSCRATCH];
static double polyy[MAXSCRATCH];
static int npoly;

/* for 1-d slices */
double slicex[1000], slicey[1000];
int slice_npts;

extern int win_h, win_w;

extern int curhist;		/* current history marker */

/*
 * variables for the text handling routine
 */
static int strx = 0, stry = 0;
static int drawx = 0, drawy = 0;
static char tmpstr[256];
static int justflag = 0;
static double si = 0.0;
static double co = 1.0;

#define MAX_AREA_POLY 200
int narea_pts = 0;

double area_polyx[MAX_AREA_POLY];
double area_polyy[MAX_AREA_POLY];
int iarea_polyx[MAX_AREA_POLY];
int iarea_polyy[MAX_AREA_POLY];
int iax[MAX_AREA_POLY];
int iay[MAX_AREA_POLY];

void refresh(Widget w, XtPointer clid, XtPointer calld)
{
    extern int batchflag;
    extern char batchfile[];

    if (debuglevel == 1) {
	printf("called refresh %d\n", inwin);
    }
    set_window(w);
    if (!inwin) {
	inwin = 1;
	my_doublebuffer(0);
	initgraphics(0);
	set_up_world(cg);
	display_image(0);
	my_doublebuffer(1);
	if (batchflag) {
	    runbatch(batchfile);
/*
   exit(1);
 */
	}
    } else {
	refresh_from_backpix();
    }
}

void resize(Widget w, XtPointer clid, XtPointer calld)
{
    Arg args;
    double wx1, wx2, wy1, wy2;

    if (debuglevel == 1) {
	printf("called resize %d\n", inwin);
    }
    if (inwin) {
	set_window(w);
	resize_backpix();
	my_doublebuffer(0);
	initgraphics(0);
	set_up_world(cg);
	display_image(0);
	my_doublebuffer(1);
    }
}

void setpointer(int x, int y)
{
    XWarpPointer(disp, None, cwin, 0, None, (unsigned int) win_w, (unsigned int) win_h, x, y);
}

/*
 * draw the graph focus indicators
 */
void draw_focus(int gno)
{
    int ix1, iy1, ix2, iy2;
    extern int draw_focus_flag;
    void set_stack_message(void);

    set_stack_message();
    if (draw_focus_flag == ON) {
	world2deviceabs(g[gno].w.xg1, g[gno].w.yg1, &ix1, &iy1);
	world2deviceabs(g[gno].w.xg2, g[gno].w.yg2, &ix2, &iy2);
	XFillRectangle(disp, cwin, gcxor, ix1 - 5, iy1 - 5, 10, 10);
	XFillRectangle(disp, cwin, gcxor, ix1 - 5, iy2 - 5, 10, 10);
	XFillRectangle(disp, cwin, gcxor, ix2 - 5, iy2 - 5, 10, 10);
	XFillRectangle(disp, cwin, gcxor, ix2 - 5, iy1 - 5, 10, 10);
/* TODO
   XFillRectangle(disp, cwin, gcxor, (ix1 + ix2) / 2 - 5, iy1 - 5, 10, 10);
   XFillRectangle(disp, cwin, gcxor, (ix1 + ix2) / 2 - 5, iy2 - 5, 10, 10);
 */
    }
}

/*
 * rubber band line
 */
void select_line(int x1, int y1, int x2, int y2)
{
    XDrawLine(disp, cwin, gcxor, x1, y1, x2, y2);
}

/*
 * draw a box on the display
 */
void draw_rectangle(int x1, int y1, int x2, int y2)
{
    XDrawRectangle(disp, cwin, gc, x1, y1, x2, y2);
}

static unsigned char solidline[1] = { 1 };
static unsigned char longdashed[2] = { 7, 7 };

/*
 * draw an xor'ed box
 */
void select_region(int x1, int y1, int x2, int y2)
{
    int dx, dy;

    dx = x2 - x1;
    dy = y2 - y1;
/*
   XSetLineAttributes(disp, gcxor, 1, LineOnOffDash, CapButt, JoinRound);
   XSetDashes(disp, gcxor, 0, longdashed, 2);
 */

    if (dx < 0) {
	iswap(&x1, &x2);
	dx = -dx;
    }
    if (dy < 0) {
	iswap(&y1, &y2);
	dy = -dy;
    }
    XDrawRectangle(disp, cwin, gcxor, x1, y1, dx, dy);
/*
   XSetLineAttributes(disp, gcxor, 0, LineOnOffDash, CapButt, JoinRound);
   XSetDashes(disp, gcxor, 0, solidline, 1);
 */
}

/*
 * switch on the area calculator
 */
void do_select_area(void)
{
    set_action(0);
    narea_pts = 0;
    set_action(COMP_AREA);
}

/*
 * switch on the perimeter calculator
 */
void do_select_peri(void)
{
    set_action(0);
    narea_pts = 0;
    set_action(COMP_PERIMETER);
}

/*
 * locator on main_panel
 */
void getpoints(int x, int y)
{
    double wx, wy, xtmp, ytmp;
    double dsx = 0.0, dsy = 0.0;
    int newg;
    char buf[256];
    extern char locator_format[];
    extern Widget loclab;
    extern XmString string;

    device2world(x, y, &wx, &wy);
    if (g[cg].pointset) {
	dsx = g[cg].dsx;
	dsy = g[cg].dsy;
    }
    if (focus_policy == FOLLOWS) {
	if ((newg = iscontained(cg, wx, wy)) != cg) {
	    draw_focus(cg);
	    cg = newg;
	    defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, islogx(cg), islogy(cg));
	    viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
	    draw_focus(cg);
	    make_format(cg);
	    device2world(x, y, &wx, &wy);
	    update_all(cg);
	}
    }
    if (!go_locateflag) {
	return;
    }
    switch (g[cg].pt_type) {
    case 0:
	xtmp = wx;
	ytmp = wy;
	{
	    char s1[30], s2[30];
	    int form = g[cg].fx;

	    create_ticklabel(form, g[cg].px, wx, s1);
	    form = g[cg].fy;
	    create_ticklabel(form, g[cg].py, wy, s2);
	    sprintf(buf, "G%1d: X, Y = [%s, %s]", cg, s1, s2);
	}
	break;
    case 1:
	xtmp = wx - dsx;
	ytmp = wy - dsy;
	sprintf(buf, locator_format, cg, xtmp, ytmp);
	break;
    case 2:
	xtmp = hypot(dsx - wx, dsy - wy);
	ytmp = 0.0;
	sprintf(buf, locator_format, cg, xtmp, ytmp);
	break;
    case 3:
	if (dsx - wx != 0.0 || dsy - wy != 0.0) {
	    xtmp = hypot(dsx - wx, dsy - wy);
	    ytmp = 180.0 + 180.0 / M_PI * atan2(dsy - wy, dsx - wx);
	    sprintf(buf, locator_format, cg, xtmp, ytmp);
	} else {
	    sprintf(buf, "ERROR: dx = dy = 0.0");
	}
	break;
    case 4:
	xtmp = xconv(wx);
	ytmp = yconv(wy);
	sprintf(buf, locator_format, cg, xtmp, ytmp);
	break;
    case 5:
	sprintf(buf, locator_format, cg, x, y);
	break;
    }
    XmStringFree(string);
    string = XmStringCreateLtoR(buf, charset);
    XtVaSetValues(loclab, XmNlabelString, string, NULL);
}


/*
 * informative message at the lower left of the canvas
 */
void write_mode_str(char *buf)
{
    Arg al;
    char tmpb[256];
    extern XmString mstring;
    extern Widget write_mode_item;

    if (buf == NULL) {
	strcpy(tmpb, "Idle...");
    } else {
	strcpy(tmpb, buf);
    }
    XmStringFree(mstring);
    mstring = XmStringCreateLtoR(tmpb, charset);
    XtSetArg(al, XmNlabelString, mstring);
    XtSetValues(write_mode_item, &al, 1);
    XmUpdateDisplay(write_mode_item);
}

/*
 * draw a cursor for text writing
 * TODO: fix the rotation problems (cursor doesn't track)
 */
void update_text_cursor(char *s, int x, int y)
{
    int hgt, tx, xtx, ytx, xhgt, yhgt;

    hgt = stringextenty(charsize * xlibcharsize, "N") / 2;
    tx = stringextentx(charsize * xlibcharsize, s);
    xtx = (int) tx *co;
    ytx = (int) tx *si;

    xhgt = (int) -hgt * si;
    yhgt = (int) hgt *co;

/*    select_line(x + tx, win_h - y + hgt, x + tx, win_h - y - hgt); */
    select_line(x + xtx + xhgt, win_h - (y + ytx + yhgt), x + xtx - xhgt, win_h - (y + ytx - yhgt));
}

/*
 * set the action_flag to the desired action (actions are
 * defined in defines.h), if 0 then cleanup the results
 * from previous actions.
 */
set_action(int act)
{
    int i, oldval;
    extern Widget canvas;

    oldval = action_flag;
    set_window(canvas);
    if (action_flag == STR_LOC) {
	double wx, wy;

	update_text_cursor(tmpstr, strx, stry);
	if (tmpstr[0]) {
	    device2world(strx, win_h - stry, &wx, &wy);
	    define_string(tmpstr, wx, wy);
	    tmpstr[0] = 0;
	}
    }
    if ((action_flag = act) == 0) {	/* clean up */
	set_cursor(-1);
	npoly = 0;
	nipoly = 0;
	write_mode_str("Idle ...");
	if (rectflag) {
	    select_region(sx, sy, old_x, old_y);
	    rectflag = 0;
	}
	if (rubber_flag) {
	    select_line(sx, sy, old_x, old_y);
	    rubber_flag = 0;
	}
    } else {
	switch (act) {
	case DEFINE_REGION:
	    write_mode_str("Pick points to define a polygonal region");
	    set_cursor(0);
	    break;
	case QUERY_TRANSECT:
	    set_cursor(1);
	    write_mode_str("Click near the node of interest");
	    break;
	case QUERY_ELA:
	case QUERY_ADCIRC:
	case QUERY_TEANL:
	    set_cursor(1);
	    write_mode_str("Click near the node of interest");
	    break;
	case QUERY_ADCIRC_MAXELEV:
	    if (maxelevfp == NULL) {
		maxelevfp = fopen("maxelev.dat", "a");
	    }
	    set_cursor(1);
	    write_mode_str("Click near the node of interest");
	    break;
	case DEL_SAMPLE_ADCIRC:
	    set_cursor(1);
	    write_mode_str("Click near the node to remove sample");
	    break;
	case DEL_SAMPLE_TEANL:
	    set_cursor(1);
	    write_mode_str("Click near the node to remove sample");
	    break;
	case SAMPLE_ADCIRC:
	    set_cursor(1);
	    write_mode_str("Click near a node to use as a sample");
	    break;
	case SAMPLE_TEANL:
	    set_cursor(1);
	    write_mode_str("Click near a node to use as a sample");
	    break;
	case GET_NEAREST_NODE:
	    set_cursor(1);
	    write_mode_str("Click near a node");
	    break;
	case GET_NEAREST_ELEMENT:
	    write_mode_str("Click near an element");
	    set_cursor(1);
	    break;
	case GET_ELEMENT:
	    write_mode_str("Click inside an element");
	    set_cursor(1);
	    break;
	case GET_GRID_DEPTH:
	    write_mode_str("Click at a point inside the domain");
	    set_cursor(1);
	    break;
	case GET_BACK_DEPTH:
	    set_cursor(1);
	    break;
	case GET_ALL:
	    set_cursor(1);
	    break;
	case ZOOM_1ST:
	    set_cursor(0);
	    write_mode_str("Click on one corner of a rectangle to enlarge");
	    break;
	case ZOOM_2ND:
	    write_mode_str("Click on the opposite corner of rectangle to enlarge");
	    break;
	case SLICE_BOX1ST:
	    set_cursor(0);
	    write_mode_str("Click at the location of the slice marker");
	    break;
	case ZOOM_BOX1ST:
	    set_cursor(0);
	    write_mode_str("Click on one corner the rectangle to enlarge");
	    break;
	case ZOOM_BOX2ND:
	    set_cursor(0);
	    write_mode_str("Click on the opposite corner of rectangle to enlarge");
	    break;
	case ZOOM_BOX3RD:
	    set_cursor(0);
	    write_mode_str("Click at the location of the zoom box");
	    break;
	case SLICE_LINE1ST:
	    set_cursor(0);
	    write_mode_str("Left click at the beginning of the slice line");
	    break;
	case SLICE_LINE2ND:
	    set_cursor(0);
	    write_mode_str("Left click at the end of the slice line");
	    break;
	case SLICE_POLY:
	    set_cursor(0);
	    write_mode_str("Click at a point, middle button to accept, right button to cancel");
	    break;
	case COMP_AREA:
	    write_mode_str("Click on points of a polygon, right button to stop");
	    set_cursor(0);
	    break;
	case COMP_PERIMETER:
	    write_mode_str("Click on points of a polygon, right button to stop");
	    set_cursor(0);
	    break;
	case SEL_POINT:
	    set_cursor(0);
	    break;
	case QUERY_STATION:
	    set_cursor(0);
	    write_mode_str("Click near the station of interest");
	    break;
	case QUERY_TIDESTATION:
	    set_cursor(0);
	    write_mode_str("Click near the tidal station of interest");
	    break;
	case PLACE_TIDESTATION:
	    set_cursor(0);
	    write_mode_str("Click near the tidal station of interest");
	    break;
	case PLACE_VSCALE:
	    set_cursor(0);
	    write_mode_str("Click on the location of the velocity scale legend");
	    break;
	case PLACE_WSCALE:
	    set_cursor(0);
	    write_mode_str("Click on the location of the wind scale legend");
	    break;
	case PLACE_ISOLINES_LEGEND:
	    set_cursor(0);
	    write_mode_str("Click on the location of the top of the legend");
	    break;
	case PLACE_DROGUE:
	    set_cursor(0);
	    write_mode_str("Click in the domain to place a drogue");
	    break;
	case DELETE_DROGUE:
	    set_cursor(0);
	    write_mode_str("Click near the drogue to delete");
	    break;
	case PICK_DROGUE:
	    set_cursor(0);
	    write_mode_str("Click near the drogue to pick");
	    break;
	case PICK_DROGUE_COLOR:
	    set_cursor(0);
	    write_mode_str("Click near the drogue to set color");
	    break;
	case PICK_DROGUE_REGION:
	    set_cursor(0);
	    write_mode_str("Click at a point of the region");
	    break;
	case PLACE_TIMELINE:
	    set_cursor(0);
	    write_mode_str("Click at the location of the Time line");
	    break;
	case TIMEINFO:
	    set_cursor(0);
	    write_mode_str("Click on the location of the time info string");
	    break;
	case PLACE_CLOCK:
	    set_cursor(0);
	    write_mode_str("Click at the location of the Tidal clock");
	    break;
	case PLACE_TEANL_ELEV1ST:
	    set_cursor(0);
	    write_mode_str("Click on the node for the elevation marker");
	    break;
	case PLACE_TEANL_ELEV2ND:
	    write_mode_str("Click on the location of the elevation marker");
	    break;
	case PLACE_ADCIRC_ELEV1ST:
	    set_cursor(0);
	    write_mode_str("Click on the node for the elevation marker");
	    break;
	case PLACE_ADCIRC_ELEV2ND:
	    write_mode_str("Click on the location of the elevation marker");
	    break;
	case EDIT_ADCIRC_ELEV:
	case EDIT_TEANL_ELEV:
	    set_cursor(0);
	    write_mode_str("Click near the location of the elevation marker to edit");
	    break;
	case PICK_ADCIRC3D_TRANSECT:
	    set_cursor(0);
	    write_mode_str("Click at the location to use for the sample");
	    break;
	case PLACE_ADCIRC3D:
	    set_cursor(0);
	    write_mode_str("Click at the location to place the marker");
	    break;
	case PICK_ADCIRC3D_XY:
	    set_cursor(0);
	    write_mode_str("Click at the location to use for the sample");
	    break;
	case PICK_ADCIRC3D_NODE:
	    set_cursor(0);
	    write_mode_str("Click at the location to use for the sample");
	    break;
	case QUERY_ADCIRC3D_XY:
	    set_cursor(0);
	    write_mode_str("Click at the location to query the data set.");
	    break;
	case QUERY_ADCIRC3D_NODE:
	    set_cursor(0);
	    write_mode_str("Click at the location to query the data set.");
	    break;
	case PLACE_HIST1ST:
	    set_cursor(0);
	    write_mode_str("Click on the node for the time history marker");
	    break;
	case PLACE_HIST2ND:
	    write_mode_str("Click on the location of the time history marker");
	    break;
	case VIEW_1ST:
	    set_cursor(0);
	    write_mode_str("Pick first corner of viewport");
	    break;
	case VIEW_2ND:
	    write_mode_str("Pick second corner of viewport");
	    break;
	case MAKE_BOX_1ST:
	    set_cursor(0);
	    write_mode_str("First corner of box");
	    break;
	case MAKE_BOX_2ND:
	    write_mode_str("Second corner of box");
	    break;
	case STR_LOC1ST:
	    set_cursor(0);
	    write_mode_str("Pick start of text line");
	    break;
	case STR_LOC2ND:
	    write_mode_str("Pick end of text line");
	    break;
	case MAKE_LINE_1ST:
	    set_cursor(0);
	    write_mode_str("Pick beginning of line");
	    break;
	case MAKE_LINE_2ND:
	    write_mode_str("Pick end of line");
	    break;
	case STR_EDIT:
	    set_cursor(2);
	    write_mode_str("Edit string");
	    break;
	case STR_LOC:
	    set_cursor(2);
	    write_mode_str("Pick beginning of text");
	    break;
	case EDIT_OBJECT:
	    set_cursor(2);
	    write_mode_str("Click near a text line, box, or line to edit");
	    break;
	case DEL_OBJECT:
	    set_cursor(3);
	    write_mode_str("Click near a text line, box, or line to delete");
	    break;
	default:
	    break;
	}
    }
}

/*
 * update string drawn on the canvas
 */
void do_text_string(int op, int c)
{
    char stmp[2];

    drawx = strx;
    drawy = stry;

    update_text_cursor(tmpstr, drawx, drawy);
    set_write_mode(0);
    dispstrxlib(drawx, drawy, sysstr.rot, tmpstr, justflag, 0);
    switch (op) {
    case 0:
	if (strlen(tmpstr) > 0) {
	    tmpstr[strlen(tmpstr) - 1] = 0;
	}
	break;
    case 1:
	sprintf(stmp, "%c", c);
	strcat(tmpstr, stmp);
	break;
    case 2:
	break;
    }
    set_write_mode(1);
    dispstrxlib(drawx, drawy, sysstr.rot, tmpstr, justflag, 0);
    update_text_cursor(tmpstr, drawx, drawy);
}

/*
 * canvas event proc
 */
void my_proc(Widget w, XtPointer data, XEvent * event)
{
    static int x, y, boxno, lineno, stringno;
    static double wx1, wx2, wy1, wy2;
    static double wx, wy, dx, dy;
    static int ty, no, c;
    double xconv(), yconv();
    static KeySym keys;
    static XComposeStatus compose;

    x = event->xmotion.x;
    y = event->xmotion.y;
/*
 * hot keys
 */
    set_window(w);
    switch (event->type) {
    case KeyPress:
	buf[0] = 0;
	XLookupString((XKeyEvent *) event, buf, 1, &keys, &compose);
	switch (c = buf[0]) {
	case 1:		/* ^A */
	    break;
	case 2:		/* ^B */
	    break;
	case 3:		/* ^C */
	    break;
	case 4:		/* ^D */
	    break;
	case 5:		/* ^E */
	    break;
	case 6:		/* ^F */
	    break;
	case 7:		/* ^G */
	    break;
	    /* stay off 8 (^H) - needed by text routines */
	case 12:		/* ^L */
	    break;
	case 14:		/* ^N */
	    break;
	case 16:		/* ^P */
	    break;
	case 18:		/* ^R */
	    break;
	case 19:		/* ^S */
	    break;
	case 20:		/* ^T */
	    break;
	case 22:		/* ^V */
	    break;
	case 23:		/* ^W */
	    break;
	case 24:		/* ^X */
	    break;
	case 26:		/* ^Z */
	    set_action(0);
	    set_action(ZOOM_1ST);
	    break;
	case 8:
	case 127:
	    if (action_flag == STR_LOC) {
		do_text_string(0, 0);
	    }
	    break;
	case '\r':
	case '\n':
	    if (action_flag == STR_LOC) {
		int itmp;

		update_text_cursor(tmpstr, drawx, drawy);
		if (tmpstr[0]) {
		    device2world(drawx, win_h - drawy, &wx, &wy);
		    define_string(tmpstr, wx, wy);
		}
		itmp = (int) (1.25 * stringextenty(charsize * xlibcharsize, "Ny"));
		strx = strx + si * itmp;
		stry = stry - co * itmp;
		tmpstr[0] = 0;
		update_text_cursor(tmpstr, strx, stry);
	    }
	    break;
	default:
	    if (action_flag == STR_LOC) {
		if (c >= 32 && c < 128) {
		    do_text_string(1, c);
		}
	    }
	    break;
	}
	return;
	break;
    case EnterNotify:
	set_window(w);
	defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, islogx(cg), islogy(cg));
	viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
	break;
    case LeaveNotify:
	break;
    case ButtonPress:
	switch (event->xbutton.button) {

	case Button3:
	    switch (action_flag) {
	    case DEFINE_REGION:
		nregion = 0;
		region_flag = 0;
		break;
	    case QUERY_ADCIRC_MAXELEV:
		if (maxelevfp != NULL) {
		    fclose(maxelevfp);
		    maxelevfp = NULL;
		}
		break;
	    case COMP_AREA:
	    case COMP_PERIMETER:
		if (nipoly >= 3) {
		    int i;

		    for (i = 0; i < nipoly; i++) {
		    }
		}
		break;
	    default:
		break;
	    }
	    npoly = 0;
	    nipoly = 0;
	    slice_npts = 0;
	    set_action(0);
	    break;
	case Button1:
	    if (action_flag == 0) {
		if (focus_policy == CLICK) {
		    int newg;

		    device2world(x, y, &wx, &wy);
		    if ((newg = iscontained(cg, wx, wy)) != cg) {
			draw_focus(cg);
			cg = newg;
			defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, islogx(cg), islogy(cg));
			viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
			draw_focus(cg);
			make_format(cg);
			device2world(x, y, &wx, &wy);
			update_all(cg);
		    }
		}
	    }
	    c = go_locateflag;
	    go_locateflag = TRUE;
	    getpoints(x, y);
	    go_locateflag = c;
	    switch (action_flag) {
	    case DEL_OBJECT:	/* delete a box or a line */
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if (ty >= 0) {
		    switch (ty) {
		    case BOX:
			set_write_mode(0);
			draw_box(-2, no);
			set_write_mode(1);
			kill_box(no);
			break;
		    case LINE:
			set_write_mode(0);
			draw_line(-2, no);
			set_write_mode(1);
			kill_line(no);
			break;
		    case STRING:
			set_write_mode(0);
			draw_string(-2, no);
			set_write_mode(1);
			kill_string(no);
			break;
		    }
		    set_action(DEL_OBJECT);
		} else {
		    set_action(0);
		}
		break;
/*
 * select a box or a line to move
 */
	    case MOVE_OBJECT_1ST:
		set_action(MOVE_OBJECT_2ND);
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if (ty < 0) {
		    set_action(0);
		} else {
		    switch (ty) {
		    case BOX:
			if (boxes[no].loctype == VIEW) {
			    sx = (int) (win_w * boxes[no].x1);
			    sy = (int) (win_h - win_h * boxes[no].y1);
			    xs = (int) (win_w * boxes[no].x2);
			    ys = (int) (win_h - win_h * boxes[no].y2);
			} else {
			    world2deviceabs(boxes[no].x1, boxes[no].y1, &sx, &sy);
			    world2deviceabs(boxes[no].x2, boxes[no].y2, &xs, &ys);
			}
			select_region(sx, sy, xs, ys);
			mbox_flag = 1;
			break;
		    case LINE:
			if (lines[no].loctype == VIEW) {
			    sx = (int) (win_w * lines[no].x1);
			    sy = (int) (win_h - win_h * lines[no].y1);
			    xs = (int) (win_w * lines[no].x2);
			    ys = (int) (win_h - win_h * lines[no].y2);
			} else {
			    world2deviceabs(lines[no].x1, lines[no].y1, &sx, &sy);
			    world2deviceabs(lines[no].x2, lines[no].y2, &xs, &ys);
			}
			select_line(sx, sy, xs, ys);
			mline_flag = 1;
			break;
		    case STRING:
			xs = stringextentx(charsize * xlibcharsize, pstr[no].s);
			ys = stringextenty(charsize * xlibcharsize, pstr[no].s);
			if (pstr[no].loctype == VIEW) {
			    sx = (int) (win_w * pstr[no].x);
			    sy = (int) (win_h - win_h * pstr[no].y);
			} else {
			    world2device(pstr[no].x, pstr[no].y, &sx, &sy);
			}
			xs = sx + xs;
			ys = sy + ys;
			mbox_flag = 1;
			select_region(sx, sy, xs, ys);
			break;
		    }
		}

		break;
/*
 * box has been selected and new position found
 */
	    case MOVE_OBJECT_2ND:
		dx = sx - x;
		dy = sy - y;

		set_action(0);
		sx = x;
		sy = y;
		xs = xs - dx;
		ys = ys - dy;
		device2world(sx, sy, &wx1, &wy1);
		device2world(xs, ys, &wx2, &wy2);
		switch (ty) {
		case BOX:
		    set_write_mode(0);
		    draw_box(-2, no);
		    if (boxes[no].loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			boxes[no].gno = cg;
		    }
		    boxes[no].x1 = wx1;
		    boxes[no].x2 = wx2;
		    boxes[no].y1 = wy1;
		    boxes[no].y2 = wy2;
		    set_write_mode(1);
		    draw_box(-2, no);
		    break;
		case LINE:
		    set_write_mode(0);
		    draw_line(-2, no);
		    if (lines[no].loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			lines[no].gno = cg;
		    }
		    lines[no].x1 = wx1;
		    lines[no].x2 = wx2;
		    lines[no].y1 = wy1;
		    lines[no].y2 = wy2;
		    set_write_mode(1);
		    draw_line(-2, no);
		    break;
		case STRING:
		    set_write_mode(0);
		    draw_string(-2, no);
		    if (pstr[no].loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
		    } else {
			pstr[no].gno = cg;
		    }
		    pstr[no].x = wx1;
		    pstr[no].y = wy1;
		    set_write_mode(1);
		    draw_string(-2, no);
		    break;
		}
		set_action(MOVE_OBJECT_1ST);
		break;
/*
 * select a box or a line to copy
 */
	    case COPY_OBJECT1ST:
		set_action(COPY_OBJECT2ND);
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if (ty < 0) {
		    set_action(0);
		} else {
		    switch (ty) {
		    case BOX:
			if (boxes[no].loctype == VIEW) {
			    sx = (int) (win_w * boxes[no].x1);
			    sy = (int) (win_h - win_h * boxes[no].y1);
			    xs = (int) (win_w * boxes[no].x2);
			    ys = (int) (win_h - win_h * boxes[no].y2);
			} else {
			    world2deviceabs(boxes[no].x1, boxes[no].y1, &sx, &sy);
			    world2deviceabs(boxes[no].x2, boxes[no].y2, &xs, &ys);
			}
			select_region(sx, sy, xs, ys);
			mbox_flag = 1;
			break;
		    case LINE:
			if (lines[no].loctype == VIEW) {
			    sx = (int) (win_w * lines[no].x1);
			    sy = (int) (win_h - win_h * lines[no].y1);
			    xs = (int) (win_w * lines[no].x2);
			    ys = (int) (win_h - win_h * lines[no].y2);
			} else {
			    world2deviceabs(lines[no].x1, lines[no].y1, &sx, &sy);
			    world2deviceabs(lines[no].x2, lines[no].y2, &xs, &ys);
			}
			select_line(sx, sy, xs, ys);
			mline_flag = 1;
			break;
		    case STRING:
			xs = stringextentx(charsize * xlibcharsize, pstr[no].s);
			ys = stringextenty(charsize * xlibcharsize, pstr[no].s);
			if (pstr[no].loctype == VIEW) {
			    sx = (int) (win_w * pstr[no].x);
			    sy = (int) (win_h - win_h * pstr[no].y);
			} else {
			    world2device(pstr[no].x, pstr[no].y, &sx, &sy);
			}
			xs = sx + xs;
			ys = sy + ys;
			mbox_flag = 1;
			select_region(sx, sy, xs, ys);
			break;
		    }
		}

		break;
/*
 * box has been selected and new position found
 */
	    case COPY_OBJECT2ND:
		dx = sx - x;
		dy = sy - y;

		set_action(0);
		sx = x;
		sy = y;
		xs = xs - dx;
		ys = ys - dy;
		device2world(sx, sy, &wx1, &wy1);
		device2world(xs, ys, &wx2, &wy2);
		switch (ty) {
		case BOX:
		    if ((boxno = next_box()) >= 0) {
			copy_object(ty, no, boxno);
			if (boxes[no].loctype == VIEW) {
			    wx1 = xconv(wx1);
			    wy1 = yconv(wy1);
			    wx2 = xconv(wx2);
			    wy2 = yconv(wy2);
			} else {
			    boxes[boxno].gno = cg;
			}
			boxes[boxno].x1 = wx1;
			boxes[boxno].x2 = wx2;
			boxes[boxno].y1 = wy1;
			boxes[boxno].y2 = wy2;
			draw_box(-2, boxno);
		    }
		    break;
		case LINE:
		    if ((lineno = next_line()) >= 0) {
			copy_object(ty, no, lineno);
			if (lines[no].loctype == VIEW) {
			    wx1 = xconv(wx1);
			    wy1 = yconv(wy1);
			    wx2 = xconv(wx2);
			    wy2 = yconv(wy2);
			} else {
			    lines[lineno].gno = cg;
			}
			lines[lineno].x1 = wx1;
			lines[lineno].x2 = wx2;
			lines[lineno].y1 = wy1;
			lines[lineno].y2 = wy2;
			draw_line(-2, lineno);
		    }
		    break;
		case STRING:
		    if ((stringno = next_string()) >= 0) {
			copy_object(ty, no, stringno);
			if (pstr[no].loctype == VIEW) {
			    wx1 = xconv(wx1);
			    wy1 = yconv(wy1);
			} else {
			    pstr[stringno].gno = cg;
			}
			pstr[stringno].x = wx1;
			pstr[stringno].y = wy1;
			draw_string(-2, stringno);
		    }
		    break;
		}
		set_action(COPY_OBJECT1ST);
		break;
/*
 * select a box or a line to move
 */
	    case EDIT_OBJECT:
		set_action(EDIT_OBJECT);
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if (ty < 0) {
		    set_action(0);
		} else {
		    switch (ty) {
		    case BOX:
			break;
		    case LINE:
			break;
		    case STRING:
			break;
		    }
		}

		break;
/*
 * make a new box, select first corner
 */
	    case MAKE_BOX_1ST:
		set_action(MAKE_BOX_2ND);
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(sx, sy, x, y);
		break;
/*
 * make a new box, select opposite corner
 */
	    case MAKE_BOX_2ND:
		set_action(0);
		if ((boxno = next_box()) >= 0) {
		    device2world(sx, sy, &wx1, &wy1);
		    device2world(x, y, &wx2, &wy2);
		    if (sysbox.loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			boxes[boxno].gno = cg;
		    }
		    boxes[boxno].loctype = sysbox.loctype;
		    boxes[boxno].x1 = wx1;
		    boxes[boxno].x2 = wx2;
		    boxes[boxno].y1 = wy1;
		    boxes[boxno].y2 = wy2;
		    boxes[boxno].color = sysbox.color;
		    boxes[boxno].linew = sysbox.linew;
		    boxes[boxno].lines = sysbox.lines;
		    boxes[boxno].fill = sysbox.fill;
		    boxes[boxno].fillcolor = sysbox.fillcolor;
		    boxes[boxno].fillpattern = sysbox.fillpattern;
		    draw_box(-2, boxno);
		}
		break;
/*
 * locate angled string
 */
	    case STR_LOC1ST:
		set_action(STR_LOC2ND);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		break;
	    case STR_LOC2ND:
		device2world(sx, sy, &wx1, &wy1);
		device2world(x, y, &wx2, &wy2);
		wx1 = xconv(wx1);
		wy1 = yconv(wy1);
		wx2 = xconv(wx2);
		wy2 = yconv(wy2);
		sysstr.rot = (int) ((atan2((wy2 - wy1) * win_h, (wx2 - wx1) * win_w) * 180.0 / M_PI) + 360.0) % 360;
		updatestrings();
		strings_loc_proc();
		break;
/*
 * make a new line, select start point
 */
	    case MAKE_LINE_1ST:
		sx = x;
		sy = y;
		set_action(MAKE_LINE_2ND);
		rubber_flag = 1;
		select_line(sx, sy, x, y);
		break;
/*
 * make a new line, select end point
 */
	    case MAKE_LINE_2ND:
		set_action(0);
		if ((lineno = next_line()) >= 0) {
		    device2world(sx, sy, &wx1, &wy1);
		    device2world(x, y, &wx2, &wy2);
		    if (sysline.loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			lines[lineno].gno = cg;
		    }
		    lines[lineno].loctype = sysline.loctype;
		    lines[lineno].x1 = wx1;
		    lines[lineno].x2 = wx2;
		    lines[lineno].y1 = wy1;
		    lines[lineno].y2 = wy2;
		    lines[lineno].color = sysline.color;
		    lines[lineno].lines = sysline.lines;
		    lines[lineno].linew = sysline.linew;
		    lines[lineno].arrow = sysline.arrow;
		    lines[lineno].asize = sysline.asize;
		    lines[lineno].atype = sysline.atype;
		    draw_line(-2, lineno);
		}
		break;
/*
 * Edit an existing string
 */
	    case STR_EDIT:
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if ((ty >= 0) && (ty == STRING)) {
		    int ilenx, ileny;

		    wx1 = pstr[no].x;
		    wy1 = pstr[no].y;
		    if (pstr[no].loctype == VIEW) {	/* in viewport coords */
			view2world(wx1, wy1, &wx2, &wy2);
			wx1 = wx2;
			wy1 = wy2;
		    }
		    world2device(wx1, wy1, &strx, &stry);
		    drawx = strx;
		    drawy = stry;
		    strcpy(tmpstr, pstr[no].s);
		    setcharsize(pstr[no].charsize);
		    setfont(pstr[no].font);
		    setcolor(pstr[no].color);
		    sysstr.just = pstr[no].just;
		    justflag = sysstr.just;
		    sysstr.charsize = pstr[no].charsize;
		    sysstr.font = pstr[no].font;
		    sysstr.color = pstr[no].color;
		    sysstr.linew = pstr[no].linew;
		    sysstr.rot = pstr[no].rot;
		    sysstr.loctype = pstr[no].loctype;
		    updatestrings();
		    kill_string(no);
		    si = sin(M_PI / 180.0 * sysstr.rot) * ((double) win_w) / ((double) win_h);
		    co = cos(M_PI / 180.0 * sysstr.rot);

		    ilenx = stringextentx(charsize * xlibcharsize, tmpstr);
		    ileny = stringextenty(charsize * xlibcharsize, tmpstr);

		    switch (justflag) {
		    case 1:
			strx = drawx + co * ilenx - si * ileny;
			stry = drawy + si * ilenx + co * ileny;
			break;
		    case 2:
			strx = drawx + (co * ilenx - si * ileny) / 2;
			stry = drawy + (si * ilenx + co * ileny) / 2;
			break;
		    }
		    update_text_cursor(tmpstr, drawx, drawy);
		    do_text_string(2, 0);
		    action_flag = STR_LOC;
		} else {
		    set_action(0);
		}
		break;
/*
 * locate a string on the canvas
 */
	    case STR_LOC:
		if (tmpstr[0]) {
		    device2world(strx, win_h - stry, &wx, &wy);
		    define_string(tmpstr, wx, wy);
		}
		strx = x;
		stry = win_h - y;
		drawx = strx;
		drawy = stry;
		tmpstr[0] = 0;
		define_string_defaults();
		justflag = sysstr.just;
		setcharsize(sysstr.charsize);
		xlibsetfont(sysstr.font);
		xlibsetcolor(sysstr.color);
		xlibsetlinewidth(sysstr.linew);
		si = sin(M_PI / 180.0 * sysstr.rot) * ((double) win_w) / ((double) win_h);
		co = cos(M_PI / 180.0 * sysstr.rot);
		update_text_cursor(tmpstr, strx, stry);
		break;
/*
 * select a reference point for the locator in main_panel
 */
	    case SEL_POINT:
		device2world(x, y, &wx, &wy);
		g[cg].pointset = TRUE;
		g[cg].dsx = wx;
		g[cg].dsy = wy;
		my_frontbuffer(1);
		draw_ref_point(cg);
		my_frontbuffer(0);
		update_locator_items(cg);
		set_action(0);
		break;
/*
 * set one corner of zoom
 */
	    case ZOOM_1ST:
		set_action(ZOOM_2ND);
		write_mode_str("Click on the opposite corner of rectangle to enlarge");
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(x, y, x, y);
		break;
/*
 * set opposing corner of zoom
 */
	    case ZOOM_2ND:
		set_action(0);
		select_region(sx, sy, old_x, old_y);
		go_blowup(sx, sy, old_x, old_y);
		break;
/*
 * set one corner of viewport
 */
	    case VIEW_1ST:
		set_action(VIEW_2ND);
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(x, y, x, y);
		break;
/*
 * set opposing corner of viewport
 */
	    case VIEW_2ND:
		{
		    double vx1, vx2, vy1, vy2;

		    set_action(0);
		    select_region(sx, sy, old_x, old_y);
		    if (sx == old_x || sy == old_y) {
			errwin("Viewport size incorrect, not changed");
		    } else {
			device2world(sx, sy, &wx1, &wy1);
			device2world(old_x, old_y, &wx2, &wy2);
			world2view(wx1, wy1, &vx1, &vy1);
			world2view(wx2, wy2, &vx2, &vy2);
			if (vx1 > vx2) {
			    fswap(&vx1, &vx2);
			}
			if (vy1 > vy2) {
			    fswap(&vy1, &vy2);
			}
			g[cg].v.xv1 = vx1;
			g[cg].v.xv2 = vx2;
			g[cg].v.yv1 = vy1;
			g[cg].v.yv2 = vy2;
			if (g[cg].type == XYFIXED) {
			    set_up_world(cg);
			}
			update_view(cg);
			drawgraph();
		    }
		}
		break;
/*
 * locate various things
 */
	    case QUERY_ELA:
		device2world(x, y, &wx, &wy);
		find_nearest_node(elaconc[curela].grid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    int c = get_current_step();
		    wx = grid[elaconc[curela].grid].xord[ind1];
		    wy = grid[elaconc[curela].grid].yord[ind1];
		    sprintf(buf, "Node %d: [%.9lg]", ind1 + 1, elaconc[curela].data[c].c[ind1]);
		} else {
		    strcpy(buf, "No node found!");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(QUERY_ELA);
		break;
	    case QUERY_TRANSECT:
		device2world(x, y, &wx, &wy);
		if (get_current_step() < trans[curtrans].nsteps) {
		    int c = get_current_step();
		FindNearestNode(&trans[curtrans].g[c], wx, wy, &ind1);
		if (ind1 >= 0) {
		    wx = trans[curtrans].g[c].xord[ind1];
		    wy = trans[curtrans].g[c].yord[ind1];
		    if (trans[curtrans].g[c].depth != NULL) {
			sprintf(buf, "Node %d: s = [%.9lg]", 
				ind1 + 1, trans[curtrans].g[c].depth[ind1]);
		    }
		} else {
		    strcpy(buf, "No node found!");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(QUERY_TRANSECT);
		}
		break;
	    case QUERY_ADCIRC:
		device2world(x, y, &wx, &wy);
		FindNearestNode(&flowt[curadcirc].g, wx, wy, &ind1);
		if (ind1 >= 0) {
		    int c = get_current_step();
		    wx = flowt[curadcirc].g.xord[ind1];
		    wy = flowt[curadcirc].g.yord[ind1];
		    if (flowt[curadcirc].f[c].e != NULL && flowt[curadcirc].f[c].u != NULL) {
			sprintf(buf, "Node %d: e, u, v = [%.9lg, %.9lg, %.9lg]", ind1 + 1, flowt[curadcirc].f[c].e[ind1], flowt[curadcirc].f[c].u[ind1], flowt[curadcirc].f[c].v[ind1]);
		    } else if (flowt[curadcirc].f[c].u != NULL) {
			sprintf(buf, "Node %d: u, v = [%.9lg, %.9lg]", ind1 + 1, flowt[curadcirc].f[c].u[ind1], flowt[curadcirc].f[c].v[ind1]);
		    } else if (flowt[curadcirc].f[c].e != NULL) {
			sprintf(buf, "Node %d: e = [%.9lg]", ind1 + 1, flowt[curadcirc].f[c].e[ind1]);
		    }
		} else {
		    strcpy(buf, "No node found!");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(QUERY_ADCIRC);
		break;
	    case QUERY_ADCIRC_MAXELEV:
		device2world(x, y, &wx, &wy);
		FindNearestNode(&flowt[curadcirc].g, wx, wy, &ind1);
		if (ind1 >= 0) {
		    wx = flowt[curadcirc].g.xord[ind1];
		    wy = flowt[curadcirc].g.yord[ind1];
		    if (compute_adcirc_maxelev(curadcirc, 0)) {
			sprintf(buf, "%d:%4.2lf\n", ind1 + 1, flowt[curadcirc].global_emax[ind1]);
			writestr(wx, wy, 0, 0, buf);
			fprintf(maxelevfp, "s %lf %lf %.2lf\n", wx, wy, flowt[curadcirc].global_emax[ind1]);
			region_flag = 1;
			regionx[nregion] = wx;
			regiony[nregion] = wy;
			nregion++;
		    }
		} else {
		    strcpy(buf, "No node found!");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(QUERY_ADCIRC_MAXELEV);
		break;
	    case SAMPLE_ADCIRC:
		device2world(x, y, &wx, &wy);
		find_nearest_node(g[cg].curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    wx = grid[g[cg].curgrid].xord[ind1];
		    wy = grid[g[cg].curgrid].yord[ind1];
		    my_frontbuffer(1);
		    my_circle(wx, wy);
		    my_frontbuffer(0);
		    add_sample_node(ind1);
		} else {
		    sprintf(buf, "Error locating node");
		}
		set_action(SAMPLE_ADCIRC);
		break;
	    case DEL_SAMPLE_ADCIRC:
		device2world(x, y, &wx, &wy);
		find_nearest_node(g[cg].curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    del_sample_node(ind1);
		} else {
		    sprintf(buf, "Error locating node");
		}
		set_action(DEL_SAMPLE_ADCIRC);
		break;
	    case GET_NEAREST_NODE:
		device2world(x, y, &wx, &wy);
		find_nearest_node(g[cg].curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    wx = grid[g[cg].curgrid].xord[ind1];
		    wy = grid[g[cg].curgrid].yord[ind1];
		    sprintf(buf, "Node %d: [%.2lf,%.2lf], depth = %.3lf", ind1 + 1, wx, wy, grid[g[cg].curgrid].depth[ind1]);
		} else {
		    sprintf(buf, "Error locating node");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(GET_NEAREST_NODE);
		break;
/*
 * Report on the nearest element to the mouse location
 */
	    case GET_NEAREST_ELEMENT:
		device2world(x, y, &wx, &wy);
		find_nearest_element(g[cg].curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    n0 = grid[g[cg].curgrid].icon[ind1].nl[0] + 1;
		    n1 = grid[g[cg].curgrid].icon[ind1].nl[1] + 1;
		    n2 = grid[g[cg].curgrid].icon[ind1].nl[2] + 1;
		    n3 = grid[g[cg].curgrid].icon[ind1].nl[3] + 1;
		    wx = area(g[cg].curgrid, ind1);
		    if (wx < 0.0) {
			wy = sqrt(-wx / M_PI);
		    } else {
			wy = sqrt(wx / M_PI);
		    }
		    if (grid[g[cg].curgrid].icon[ind1].ngeom == 4) {
			sprintf(buf, "Element %d: [%d %d %d %d], A = %.2lf, ER = %.2lf", ind1 + 1, n0, n1, n2, n3, wx, wy);
		    } else {
			sprintf(buf, "Element %d: [%d %d %d], A = %.2lf, ER = %.2lf", ind1 + 1, n0, n1, n2, wx, wy);
		    }
		} else {
		    sprintf(buf, "Error locating element");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(GET_NEAREST_ELEMENT);
		break;
/*
 * Report on the element containing the current mouse location
 */
	    case GET_ELEMENT:
		device2world(x, y, &wx, &wy);
		find_element(g[cg].curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    n0 = grid[g[cg].curgrid].icon[ind1].nl[0] + 1;
		    n1 = grid[g[cg].curgrid].icon[ind1].nl[1] + 1;
		    n2 = grid[g[cg].curgrid].icon[ind1].nl[2] + 1;
		    n3 = grid[g[cg].curgrid].icon[ind1].nl[3] + 1;
		    wx = area(g[cg].curgrid, ind1);
		    if (wx < 0.0) {
			wy = sqrt(-wx / M_PI);
		    } else {
			wy = sqrt(wx / M_PI);
		    }
		    if (grid[g[cg].curgrid].icon[ind1].ngeom == 4) {
			sprintf(buf, "Element %d: [%d %d %d %d], A = %.2lf, ER = %.2lf", ind1 + 1, n0, n1, n2, n3, wx, wy);
		    } else {
			sprintf(buf, "Element %d: [%d %d %d], A = %.2lf, ER = %.2lf", ind1 + 1, n0, n1, n2, wx, wy);
		    }
		} else {
		    sprintf(buf, "No element found");
		}
		xv_setstr(locate_grid_item, buf);
		set_action(GET_ELEMENT);
		break;
/*
 * Verbose report on the node and element
 */
	    case GET_ALL:
		device2world(x, y, &wx, &wy);
		break;
/*
 * Get the depth at the mouse location
 */
	    case GET_GRID_DEPTH:
		device2world(x, y, &wx, &wy);
		find_element(g[cg].curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {

		    sprintf(buf, "Element %d: at (%lf, %lf) Depth = %lf", ind1 + 1, wx, wy, get_depth_element(g[cg].curgrid, ind1, wx, wy));
		    xv_setstr(locate_grid_item, buf);
		    set_action(GET_GRID_DEPTH);
		}
		break;
/*
 * Define a region for operations
 */
	    case DEFINE_REGION:
		{
		    device2world(x, y, &regionx[nregion], &regiony[nregion]);
		    ipolyx[npoly] = x;
		    ipolyy[npoly] = y;
		    nipoly++;
		    nregion++;
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(DEFINE_REGION);
		}
		break;
/*
 * define a polygonal region
 */
	    case DEF_REGION:
		device2world(x, y, &area_polyx[region_pts], &area_polyy[region_pts]);
		iax[region_pts] = x;
		iay[region_pts] = y;
		region_pts++;
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(DEF_REGION);
		break;
/*
 * Compute the area of a polygon
 */
	    case COMP_AREA:
		{
		    double area, comp_area(int n, double *x, double *y);
		    extern Widget arealab;
		    extern XmString astring;

		    device2world(x, y, &area_polyx[narea_pts], &area_polyy[narea_pts]);
		    iax[narea_pts] = x;
		    iay[narea_pts] = y;
		    narea_pts++;
		    if (narea_pts <= 2) {
			area = 0.0;
		    } else {
			area = comp_area(narea_pts, area_polyx, area_polyy);
		    }
		    sprintf(buf, "[%lf]", fabs(area));
		    XmStringFree(astring);
		    astring = XmStringCreateLtoR(buf, charset);
		    XtVaSetValues(arealab, XmNlabelString, astring, NULL);
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(COMP_AREA);
		}
		break;
	    case COMP_PERIMETER:
		{
		    double area, comp_perimeter(int n, double *x, double *y);
		    extern Widget perimlab;
		    extern XmString pstring;

		    device2world(x, y, &area_polyx[narea_pts], &area_polyy[narea_pts]);
		    iax[narea_pts] = x;
		    iay[narea_pts] = y;
		    narea_pts++;
		    if (narea_pts <= 1) {
			area = 0.0;
		    } else {
			area = comp_perimeter(narea_pts, area_polyx, area_polyy);
		    }
		    sprintf(buf, "[%lf]", fabs(area));

		    XmStringFree(pstring);
		    pstring = XmStringCreateLtoR(buf, charset);
		    XtVaSetValues(perimlab, XmNlabelString, pstring, NULL);
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(COMP_PERIMETER);
		}
		break;
	    case PICK_BATH:
		rubber_flag = 1;
		sx = x;
		sy = y;
		device2world(x, y, &wx, &wy);
		slicex[slice_npts] = wx;
		slicey[slice_npts] = wy;
		slice_npts++;
		set_action(PICK_BATH);
		break;
	    case SLICE_POLY:
		device2world(x, y, &wx, &wy);
		slicex[slice_npts] = wx;
		slicey[slice_npts] = wy;
		slice_npts++;
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(SLICE_POLY);
		break;
	    case SLICE_LINE1ST:
		device2world(x, y, &wx1, &wy1);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(SLICE_LINE2ND);
		break;
	    case SLICE_LINE2ND:
		device2world(x, y, &wx2, &wy2);
		define_sliceline(cg, curslice, wx1, wy1, wx2, wy2);
		set_action(0);
		break;
	    case SLICE_BOX1ST:
		device2world(x, y, &wx, &wy);
		define_slicebox(cg, curslice, wx, wy);
		set_action(0);
		break;
	    case SLICE_BOX2ND:
		break;
	    case SLICE_BOX3RD:
		break;
/*
 * set one corner of zoom box
 */
	    case ZOOM_BOX1ST:
		set_action(ZOOM_BOX2ND);
		device2world(x, y, &wx1, &wy1);
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(x, y, x, y);
		break;
/*
 * set opposing corner of zoom box
 */
	    case ZOOM_BOX2ND:
		set_action(0);
		set_action(ZOOM_BOX3RD);
		device2world(x, y, &wx2, &wy2);
		break;
/*
 * set location of zoom box
 */
	    case ZOOM_BOX3RD:
		set_action(0);
		device2world(x, y, &wx, &wy);
		define_zoom(cg, curzoom, wx1, wy1, wx2, wy2, wx, wy);
		draw_zoom(cg, curzoom);
		break;
	    case QUERY_STATION:
		set_action(0);
		device2world(x, y, &wx, &wy);
		FindNearestStation(nsta, sta, wx, wy, &ind1);
		if (ind1 >= 0) {
/*
   create_graph2d_popup(ind1);
 */
		}
		break;
	    case QUERY_TIDESTATION:
		device2world(x, y, &wx, &wy);
		find_nearest_tidestation(wx, wy, &ind1);
		if (ind1 >= 0) {
		    set_current_tidestation(ind1);
		    set_action(QUERY_TIDESTATION);
		} else {
		    set_action(0);
		}
		break;
	    case PLACE_TIDESTATION:
		device2world(x, y, &wx, &wy);
		set_current_tidestationloc(wx, wy);
		set_action(PLACE_TIDESTATION);
		break;
	    case PLACE_DROGUE:
		{
		    extern double drogx[], drogy[];
		    extern int ndrogs;

		    set_action(PLACE_DROGUE);
		    device2world(x, y, &wx, &wy);
		    my_frontbuffer(1);
		    my_circle(wx, wy);
		    my_frontbuffer(0);
		    drogx[ndrogs] = wx;
		    drogy[ndrogs++] = wy;
		}
		break;
	    case DELETE_DROGUE:
		{
		    extern double drogx[], drogy[];
		    extern int ndrogs;

		    if (ndrogs > 0) {
			tmpind = find_nearest_drog(curdrog, wx, wy);
			device2world(x, y, &wx, &wy);
			set_action(DELETE_DROGUE);
		    } else {
			set_action(0);
		    }
		}
		break;
	    case PICK_DROGUE:
		break;
	    case PICK_DROGUE_COLOR:
		{
		    extern double drogx[], drogy[];
		    extern int ndrogs, drog_curcolor;
		    int cstep = get_current_step();

		    device2world(x, y, &wx, &wy);
		    tmpind = find_nearest_drogue(curdrog, wx, wy);
		    if (tmpind >= 0) {
			my_frontbuffer(1);
			setcolor(drog_curcolor);
			drawpolysym(&drogues[curdrog].p[cstep].x[tmpind], &drogues[curdrog].p[cstep].y[tmpind], 1, 2, 0, 1, 1.0);
			set_drogue_color(curdrog, tmpind, drog_curcolor);
			setcolor(1);
			my_frontbuffer(0);
			set_action(PICK_DROGUE_COLOR);
		    } else {
			set_action(0);
		    }
		}
		break;
	    case PICK_DROGUE_REGION:
		{
		    device2world(x, y, &area_polyx[region_pts], &area_polyy[region_pts]);
		    iax[region_pts] = x;
		    iay[region_pts] = y;
		    region_pts++;
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(PICK_DROGUE_REGION);
		}
		break;
	    case PLACE_HIST1ST:
		set_action(PLACE_HIST2ND);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		device2world(sx, sy, &wx1, &wy1);

		break;
	    case PLACE_HIST2ND:
		set_action(0);
		select_line(sx, sy, x, y);
		device2world(x, y, &wx2, &wy2);
		g[cg].hbox[curhist].x = wx1;
		g[cg].hbox[curhist].y = wy1;
		g[cg].hbox[curhist].elem = -1; /* to prime find_element */
		g[cg].hbox[curhist].locx = wx2;
		g[cg].hbox[curhist].locy = wy2;
		my_frontbuffer(1);
		drawhist(cg, curhist, timeclock.curstep, timeclock.curtime);
		my_frontbuffer(0);
		break;
	    case PLACE_TEANL_ELEV1ST:
		set_action(PLACE_TEANL_ELEV2ND);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		device2world(sx, sy, &wx1, &wy1);
		break;
	    case PLACE_TEANL_ELEV2ND:
		set_action(0);
		select_line(sx, sy, x, y);
		device2world(x, y, &wx2, &wy2);
		register_elevmarker(cg, curteanl, curteanlem, TEANL, wx1, wy1, wx2, wy2);
		set_action(PLACE_TEANL_ELEV1ST);
		break;
	    case EDIT_TEANL_ELEV:
		set_action(0);
		device2world(x, y, &wx1, &wy1);
		tmpind = getnearest_elevmarker(cg, curteanl, TEANL, wx1, wy1);
		if (tmpind >= 0) {
		    curteanlem = tmpind;
		    update_teanl_flow();
		} else {
		    errwin("No elevation marker found\n");
		}
		break;
	    case PLACE_ADCIRC_ELEV1ST:
		set_action(PLACE_ADCIRC_ELEV2ND);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		device2world(sx, sy, &wx1, &wy1);
		break;
	    case PLACE_ADCIRC_ELEV2ND:
		set_action(0);
		select_line(sx, sy, x, y);
		device2world(x, y, &wx2, &wy2);
		register_elevmarker(cg, curadcirc, curadcircem, ADCIRC, wx1, wy1, wx2, wy2);
		break;
	    case EDIT_ADCIRC_ELEV:
		set_action(0);
		device2world(x, y, &wx1, &wy1);
		tmpind = getnearest_elevmarker(cg, curadcirc, ADCIRC, wx1, wy1);
		if (tmpind >= 0) {
		    curadcircem = tmpind;
		    update_adcirc_flow();
		} else {
		    errwin("No elevation marker found\n");
		}
		break;
	    case PLACE_ADCIRC3D:
		{
		    extern int curadc3d;
		    device2world(x, y, &wx1, &wy1);
		    register_adc3dvelmarker(cg, curadc3d, wx1, wy1);
		    set_action(0);
		}
		break;
	    case PICK_ADCIRC3D_XY:
		device2world(x, y, &wx1, &wy1);
		set_adc3d_xy(wx1, wy1);
		set_action(0);
		break;
	    case PICK_ADCIRC3D_NODE:
		device2world(x, y, &wx1, &wy1);
		set_adc3d_node(wx1, wy1);
		set_action(0);
		break;
	    case QUERY_ADCIRC3D_XY:
		device2world(x, y, &wx1, &wy1);
		query_adc3d_xy(wx1, wy1);
		set_action(0);
		break;
	    case QUERY_ADCIRC3D_NODE:
		device2world(x, y, &wx1, &wy1);
		query_adc3d_node(wx1, wy1);
		set_action(0);
		break;
	    case PICK_ADCIRC3D_TRANSECT:
		device2world(x, y, &wx1, &wy1);
		add_trans_xy(wx1, wy1);
		set_action(PICK_ADCIRC3D_TRANSECT);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		device2world(sx, sy, &wx1, &wy1);
		break;
	    case PLACE_TITLE:
		set_action(0);
		break;
	    case PLACE_CLOCK:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, TIDALCLOCK, curplaceitem, wx, wy);
		break;
	    case TIMEINFO:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, TIMEINFO, curplaceitem, wx, wy);
		break;
	    case PLACE_TIMELINE:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, TIMELINE, curplaceitem, wx, wy);
		break;
	    case PLACE_MAPSCALE:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, MAPSCALE, curplaceitem, wx, wy);
		break;
	    case PLACE_VSCALE:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, VSCALE, curplaceitem, wx, wy);
		break;
	    case PLACE_WSCALE:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, WSCALE, curplaceitem, wx, wy);
		break;
	    case PLACE_FSCALE:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_loc(cg, FLUX, curplaceitem, wx, wy);
		break;
	    case PLACE_ISOLINES_LEGEND:
		device2world(x, y, &wx, &wy);
		set_loc(cg, ISOLINES, curplaceitem, wx, wy);
		set_action(0);
		break;
	    }
	    break;
	case Button2:
	    switch (action_flag) {
	    case DEFINE_REGION:
		if (nregion >= 3) {
		    region_flag = 1;
		    draw_single_region(regionx, regiony, nregion);
		} else {
		    nregion = 0;
		    region_flag = 0;
		}
		break;
	    case DEF_REGION:
		break;
	    case SLICE_POLY:
		if (slice_npts > 0) {
		    define_slicepoly(cg, curslice, slice_npts, slicex, slicey);
		}
		break;
	    case PICK_BATH:
/*
   write_slice(g[cg].curgrid, 1);
 */
		break;
	    case QUERY_ADCIRC_MAXELEV:
		if (maxelevfp != NULL) {
		    fclose(maxelevfp);
		    maxelevfp = NULL;
		}
		break;
	    default:
		break;
	    }
	    npoly = 0;
	    nipoly = 0;
	    slice_npts = 0;
	    set_action(0);
	    break;
	}
	break;
    case MotionNotify:
	x = event->xmotion.x;
	y = event->xmotion.y;
	if (event->xmotion.state & Button1MotionMask) {
	    switch (action_flag) {
	    }
	}
	if (go_locateflag) {
	    getpoints(x, y);
	}
	break;
    default:
	break;
    }
/*
 * some mouse tracking stuff
 */
    switch (action_flag) {
    case MOVE_OBJECT_2ND:
    case COPY_OBJECT2ND:
	dx = sx - x;
	dy = sy - y;

	switch (ty) {
	case BOX:
	    select_region(sx, sy, xs, ys);
	    sx = x;
	    sy = y;
	    xs = xs - dx;
	    ys = ys - dy;
	    select_region(sx, sy, xs, ys);
	    break;
	case LINE:
	    select_line(sx, sy, xs, ys);
	    sx = x;
	    sy = y;
	    xs = xs - dx;
	    ys = ys - dy;
	    select_line(sx, sy, xs, ys);
	    break;
	case STRING:
	    select_region(sx, sy, xs, ys);
	    sx = x;
	    sy = y;
	    xs = xs - dx;
	    ys = ys - dy;
	    select_region(sx, sy, xs, ys);
	    break;
	}
	break;
    case STR_LOC:
	break;
/*
   case LEG_LOC:
   break;
 */
    }
    if (rectflag) {
	select_region(sx, sy, old_x, old_y);
	select_region(sx, sy, x, y);
    }
    if (rubber_flag) {
	select_line(sx, sy, old_x, old_y);
	select_line(sx, sy, x, y);
    }
    old_x = x;
    old_y = y;
}

/*
 * draw a crosshair cursor
 */
void motion(XMotionEvent * e)
{
    if (e->type != MotionNotify)
	return;
    /* Erase the previous crosshair */
    XDrawLine(disp, cwin, gcxor, 0, cursor_oldy, win_w, cursor_oldy);
    XDrawLine(disp, cwin, gcxor, cursor_oldx, 0, cursor_oldx, win_h);

    /* Draw the new crosshair */
    cursor_oldx = e->x;
    cursor_oldy = e->y;
    XDrawLine(disp, cwin, gcxor, 0, cursor_oldy, win_w, cursor_oldy);
    XDrawLine(disp, cwin, gcxor, cursor_oldx, 0, cursor_oldx, win_h);
}

/*
 * double click detection
 */
#define CLICKINT 400

int double_click(XButtonEvent * e)
{
    static Time lastc = 0;

    if (e->time - lastc < CLICKINT) {
	return 1;
    }
    lastc = e->time;
    return 0;
}

/*
 * for world stack
 */
void set_stack_message(void)
{
    extern XmString sdstring, cystring;
    extern Widget stack_depth_item;
    extern Widget curw_item;
    Arg al;

    if (stack_depth_item) {
	sprintf(buf, " SD:%1d ", g[cg].ws_top);
	XmStringFree(sdstring);
	sdstring = XmStringCreateLtoR(buf, charset);
	XtSetArg(al, XmNlabelString, sdstring);
	XtSetValues(stack_depth_item, &al, 1);
	sprintf(buf, " CW:%1d ", g[cg].curw);
	XmStringFree(cystring);
	cystring = XmStringCreateLtoR(buf, charset);
	XtSetArg(al, XmNlabelString, cystring);
	XtSetValues(curw_item, &al, 1);
    }
}

void switch_current_graph(int gfrom, int gto)
{
    draw_focus(gfrom);
    cg = gto;
    defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, islogx(cg), islogy(cg));
    viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
    draw_focus(cg);
    make_format(cg);
    update_all(cg);
}
