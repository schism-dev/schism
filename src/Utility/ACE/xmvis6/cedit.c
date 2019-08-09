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
 * Color map editor
 *
*/

#ifndef lint
static char RCSid[] = "$Id: cedit.c,v 1.3 2006/12/27 01:02:59 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"

extern unsigned char red[], green[], blue[];

#define MAXCOLORS   48
#define MAXCOLORVAL 65535

static Widget cedit_frame, cedit_panel;

extern Display *disp;
#ifdef USE_GL
Colormap cmap;
#else
extern Colormap cmap;
#endif
static int r[256], g[256], b[256];
extern long colors[];

static WidgetList wl;		/* for colorbar labels */

static int or[256], og[256], ob[256];
static int sr[256], sg[256], sb[256];

static int or1[16], og1[16], ob1[16];
static int or2[16], og2[16], ob2[16];
static int or3[16], og3[16], ob3[16];

static int sr1[16], sg1[16], sb1[16];
static int sr2[16], sg2[16], sb2[16];
static int sr3[16], sg3[16], sb3[16];

static XColor current_color;
static int ncolors = 0;

static int ncur = 0;

static double huev, satv, valv;

static Widget red_slider, blue_slider, green_slider, color_number;
static Widget h_slider, s_slider, v_slider;
static Widget canvas1, canvas2;

static void slider_selected();
static void red_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data);
static void green_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data);
static void blue_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data);
static void h_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data);
static void s_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data);
static void v_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data);
static void set_current_color(Widget w, int number, XEvent * event);
static void updatesliders(void);
static Widget make_slider();
static void interp_rgb(void);
static void interp_hsv(void);
static void create_color_bar(Widget parent);
static void update_color(void);

/* used to set up XmStrings */
extern XmStringCharSet charset;
extern Widget app_shell;

static XmString redstr, bluestr, greenstr, hstr, sstr, vstr;

static Widget bb, sliders, rc, bt;
static Widget interpt;
static Widget interpf;
static Widget interptl;
static Widget interpfl;
static Widget rgbhls_fr1;
static Widget rgbhls_rc1;
static Widget rgbhls_tog1[2];

static Widget rcmapdata_dialog;	/* read data popup */
static Widget wcmapdata_dialog;	/* write data popup */

static int maxcolors = MAXCOLORS;

static int rgbflag = 1;

static double max3(double x, double y, double z)
{
    if (y > x)
	x = y;
    if (z > x)
	x = z;
    return x;
}

static double min3(double x, double y, double z)
{
    if (y < x)
	x = y;
    if (z < x)
	x = z;
    return x;
}

#define UNDEFINED -1.0

void hsv_to_rgb(double h, double s, double v, double *r, double *g, double *b)
{
    double i, f;
    double p, q, t;

    if (v < 0.001)
	v = 0.01;
    if (s < 0.001)
	s = 0.01;
    if (s == 0.0) {
	if (h == UNDEFINED) {
	    h = 0;
	    *r = *g = *b = v;
	}
    } else {
	if (h >= 360.0) {
	    h = 359.9;
	}
	h = h / 60;
	i = floor(h);
	f = h - i;
	p = v * (1 - s);
	q = v * (1 - s * f);
	t = v * (1 - s * (1 - f));

	switch ((int) i) {
	case 0:
	    *r = v, *g = t, *b = p;
	    break;
	case 1:
	    *r = q, *g = v, *b = p;
	    break;
	case 2:
	    *r = p, *g = v, *b = t;
	    break;
	case 3:
	    *r = p, *g = q, *b = v;
	    break;
	case 4:
	    *r = t, *g = p, *b = v;
	    break;
	case 5:
	    *r = v, *g = p, *b = q;
	    break;
	}
    }
}

void rgb_to_hsv(double r, double g, double b, double *h, double *s, double *v)
{
    double maxv;
    double minv;
    r /= MAXCOLORVAL;
    g /= MAXCOLORVAL;
    b /= MAXCOLORVAL;
    maxv = max3(r, g, b);
    minv = min3(r, g, b);
    *v = maxv;
    if (maxv != 0.0) {
	*s = (maxv - minv) / maxv;
    } else {
	*s = 0;
    }

    if (*s == 0.0) {
	*h = 0;
    } else {
	double delta = (maxv - minv);
	if (delta == 0.0) {
	    *h = 0;
	} else if (r == maxv) {
	    *h = (g - b) / delta;
	} else if (g == maxv) {
	    *h = 2 + (b - r) / delta;
	} else if (b == maxv) {
	    *h = 4 + (r - g) / delta;
	}
	*h = *h * 60;
	if (*h < 0) {
	    *h += 360;
	}
    }
}

static void savecolors(void)
{
    int i;
    for (i = 0; i < ncolors; i++) {
	sr[i] = r[i];
	sg[i] = g[i];
	sb[i] = b[i];
    }
    update_color();
}

static void restorecolors(void)
{
    int i;
    for (i = 0; i < ncolors; i++) {
	r[i] = sr[i];
	g[i] = sg[i];
	b[i] = sb[i];
	set_colormapdata(i, r[i], g[i], b[i]);
    }
    current_color.red = r[ncur] << 8;
    current_color.green = g[ncur] << 8;
    current_color.blue = b[ncur] << 8;
    updatesliders();
    update_color();
}

static void defaultcolors(void)
{
    int i;
    for (i = 0; i < ncolors; i++) {
	r[i] = or[i];
	g[i] = og[i];
	b[i] = ob[i];
	set_colormapdata(i, r[i], g[i], b[i]);
    }
    current_color.red = r[ncur] << 8;
    current_color.green = g[ncur] << 8;
    current_color.blue = b[ncur] << 8;
    updatesliders();
    update_color();
}

static void interp_proc(void)
{
    if (rgbflag) {
	interp_rgb();

    } else {
	interp_hsv();
    }
    if (DisplayPlanes(disp, DefaultScreen(disp)) > 8) {
	int i;
	for (i = 0; i < ncolors; i++) {
	    XtVaSetValues(wl[i], XtNbackground, colors[i], NULL);
	}
    }
}

static void interp_rgb(void)
{
    int i, cmin, cmax, dn, intf, intt;
    double rf, gf, bf, rt, gt, bt;
    double dr, dg, db;
    XColor cf, ct, cc;
    cmin = intf = atoi(XmTextGetString(interpf));
    cmax = intt = atoi(XmTextGetString(interpt));
    dn = cmax - cmin;
    cf.flags = DoRed | DoGreen | DoBlue;
    cf.pixel = colors[cmin];
    XQueryColor(disp, cmap, &cf);
    rf = cf.red;
    gf = cf.green;
    bf = cf.blue;

    ct.flags = DoRed | DoGreen | DoBlue;
    ct.pixel = colors[cmax];
    XQueryColor(disp, cmap, &ct);
    rt = ct.red;
    gt = ct.green;
    bt = ct.blue;

    dr = rt - rf;
    dg = gt - gf;
    db = bt - bf;
    savecolors();
    for (i = cmin + 1; i < cmax; i++) {
	cc.flags = DoRed | DoGreen | DoBlue;
	cc.pixel = colors[i];
	cc.red = (int) (rf + dr / dn * (i - cmin));
	cc.green = (int) (gf + dg / dn * (i - cmin));
	cc.blue = (int) (bf + db / dn * (i - cmin));
	set_colormapdata(i, cc.red >> 8, cc.green >> 8, cc.blue >> 8);
    }
}

static void interp_hsv(void)
{
    int i, cmin, cmax, dn, intf, intt;
    double r, g, b;
    double rf, gf, bf;
    double rt, gt, bt;
    double dh, ds, dv;
    double hf, sf, vf;
    double hr, sr, vr;
    double ht, st, vt;
    XColor cf, ct, cc;
    cmin = intf = atoi(XmTextGetString(interpf));
    cmax = intt = atoi(XmTextGetString(interpt));
    dn = cmax - cmin;

    cf.flags = DoRed | DoGreen | DoBlue;
    cf.pixel = colors[cmin];
    XQueryColor(disp, cmap, &cf);
    rf = cf.red;
    gf = cf.green;
    bf = cf.blue;
    rgb_to_hsv(rf, gf, bf, &hf, &sf, &vf);

    ct.flags = DoRed | DoGreen | DoBlue;
    ct.pixel = colors[cmax];
    XQueryColor(disp, cmap, &ct);
    rt = ct.red;
    gt = ct.green;
    bt = ct.blue;
    rgb_to_hsv(rt, gt, bt, &ht, &st, &vt);
    dh = ht - hf;
    ds = st - sf;
    dv = vt - vf;

    savecolors();

    for (i = cmin + 1; i < cmax; i++) {
	ht = (hf + dh / dn * (i - cmin));
	st = (sf + ds / dn * (i - cmin));
	vt = (vf + dv / dn * (i - cmin));
	cc.flags = DoRed | DoGreen | DoBlue;
	cc.pixel = colors[i];
	hsv_to_rgb(ht, st, vt, &r, &g, &b);
	cc.red = (int) (r * MAXCOLORVAL);
	cc.green = (int) (g * MAXCOLORVAL);
	cc.blue = (int) (b * MAXCOLORVAL);
	set_colormapdata(i, cc.red >> 8, cc.green >> 8, cc.blue >> 8);
    }
}

static void manage_rcmapdata_dialog(void)
{
    XtManageChild(rcmapdata_dialog);
}

static void close_rfile_popup(void)
{
    XtUnmanageChild(rcmapdata_dialog);
}

static void rcmapdata_proc(void)
{
    int i, rr, gg, bb;
    FILE *fp;
    char *s, buf[256];
    Arg args;
    XmString list_item;
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(rcmapdata_dialog, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    fp = fopen(s, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return;
    }
    while (fgets(buf, 255, fp) != NULL) {
	sscanf(buf, "%d %d %d %d", &i, &rr, &gg, &bb);
	set_colormapdata(i, rr, gg, bb);
	r[i] = rr;
	g[i] = gg;
	b[i] = bb;
    }
    fclose(fp);
    current_color.red = r[ncur] << 8;
    current_color.green = g[ncur] << 8;
    current_color.blue = b[ncur] << 8;
    updatesliders();
    update_color();
}

static void create_rfile_popup(void)
{
    Arg al[10];
    int ac;

    if (rcmapdata_dialog) {
	manage_rcmapdata_dialog();
	return;
    }
    rcmapdata_dialog = XmCreateFileSelectionDialog(app_shell, "Read colormap", NULL, 0);
    XtAddCallback(rcmapdata_dialog, XmNcancelCallback, (XtCallbackProc) close_rfile_popup, 0);
    XtAddCallback(rcmapdata_dialog, XmNokCallback, (XtCallbackProc) rcmapdata_proc, 0);
    manage_rcmapdata_dialog();
}

static void manage_wcmapdata_dialog(void)
{
    XtManageChild(wcmapdata_dialog);
}

static void close_wfile_popup(void)
{
    XtUnmanageChild(wcmapdata_dialog);
}

static void wcmapdata_proc(void)
{
    int i;
    FILE *fp;
    char *s, buf[256];
    Arg args;
    XmString list_item;
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(wcmapdata_dialog, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    fp = fopen(s, "w");
    XtFree(s);
    if (fp == NULL) {
	errwin("Unable to open file");
	return;
    }
    for (i = 0; i < maxcolors; i++) {
	fprintf(fp, "%d %d %d %d\n", i, red[i], green[i], blue[i]);
    }
    fclose(fp);
}

int write_int(int *d, int n, FILE * fout);
int read_int(int *d, int n, FILE * fout);

void ReadCMapBinary(FILE * fp)
{
    int i, cnt, ind, rr, gg, bb;
    char *s, buf[256];
    read_int(&maxcolors, 1, fp);
    read_int(r, maxcolors, fp);
    read_int(g, maxcolors, fp);
    read_int(b, maxcolors, fp);
    for (i = 0; i < maxcolors; i++) {
	set_colormapdata(i, r[i], g[i], b[i]);
	red[i] = r[i];
	green[i] = g[i];
	blue[i] = b[i];
    }
    current_color.red = r[ncur] << 8;
    current_color.green = g[ncur] << 8;
    current_color.blue = b[ncur] << 8;
}

void WriteCMapBinary(FILE * fp)
{
    int i = maxcolors;
    char *s, buf[256];
    write_int(&i, 1, fp);
    for (i = 0; i < maxcolors; i++) {
	r[i] = red[i];
	g[i] = green[i];
	b[i] = blue[i];
    }
    write_int(r, maxcolors, fp);
    write_int(g, maxcolors, fp);
    write_int(b, maxcolors, fp);
}

void WriteCMap(FILE * fp)
{
    int i = maxcolors;
    char *s, buf[256];
    for (i = 0; i < maxcolors; i++) {
	r[i] = red[i];
	g[i] = green[i];
	b[i] = blue[i];
	fprintf(fp, "colormap %d, %d, %d, %d\n", i, r[i], g[i], b[i]);
    }
}

static void create_wfile_popup(void)
{
    Arg al[10];
    int ac;

    if (wcmapdata_dialog) {
	manage_wcmapdata_dialog();
	return;
    }
    wcmapdata_dialog = XmCreateFileSelectionDialog(app_shell, "Write colormap", NULL, 0);
    XtAddCallback(wcmapdata_dialog, XmNcancelCallback, (XtCallbackProc) close_wfile_popup, 0);
    XtAddCallback(wcmapdata_dialog, XmNokCallback, (XtCallbackProc) wcmapdata_proc, 0);
    manage_wcmapdata_dialog();
}

static void readcmap_proc(void)
{
    create_rfile_popup();
}

static void writecmap_proc(void)
{
    create_wfile_popup();
}

static void setspace(Widget w, int cd)
{
    Arg warg;
    rgbflag = cd;
    if (rgbflag) {
    } else {
    }
}

static void set_cmap_proc(void)
{
}

static void quit_proc(void)
{
    XtUnmanageChild(cedit_frame);
}

void create_cedit_popup(void)
{
    Colormap def_colormap;
    int i;
    Widget fr, rb, w[4], s1, s2;
    ncolors = 48;
    if (cedit_frame) {
	XtManageChild(cedit_frame);
	return;
    }
    for (i = 0; i < ncolors; i++) {
	r[i] = or[i] = red[i];
	g[i] = og[i] = green[i];
	b[i] = ob[i] = blue[i];
    }

    cedit_frame = XmCreateDialogShell(app_shell, "Color edit", NULL, 0);
    handle_close(cedit_frame);
    cedit_panel = XmCreateRowColumn(cedit_frame, "panel", NULL, 0);

    create_color_bar(cedit_panel);

    XtCreateManagedWidget("sep", xmSeparatorWidgetClass, cedit_panel, NULL, 0);

    color_number = XtCreateManagedWidget("Current color: 0", xmLabelWidgetClass, cedit_panel, NULL, 0);

    sliders = XtVaCreateManagedWidget("sliderpanel", xmRowColumnWidgetClass, cedit_panel, XmNorientation, XmHORIZONTAL, NULL);

    s1 = XtCreateManagedWidget("s1", xmRowColumnWidgetClass, sliders, NULL, 0);
    XtCreateManagedWidget("RGB", xmLabelWidgetClass, s1, NULL, 0);

    s2 = XtCreateManagedWidget("s2", xmRowColumnWidgetClass, sliders, NULL, 0);
    XtCreateManagedWidget("HSV", xmLabelWidgetClass, s2, NULL, 0);

    red_slider = XtVaCreateManagedWidget("Red         ", xmScaleWidgetClass, s1, XmNminimum, 0, XmNmaximum, MAXCOLORVAL, XmNorientation, XmHORIZONTAL, XmNscaleWidth, 255, NULL);
    XtAddCallback(red_slider, XmNvalueChangedCallback, (XtCallbackProc) red_slider_moved, NULL);
    XtAddCallback(red_slider, XmNdragCallback, (XtCallbackProc) red_slider_moved, NULL);

    green_slider = XtVaCreateManagedWidget("Green       ", xmScaleWidgetClass, s1, XmNminimum, 0, XmNmaximum, MAXCOLORVAL, XmNorientation, XmHORIZONTAL, XmNscaleWidth, 255, NULL);
    XtAddCallback(green_slider, XmNvalueChangedCallback, (XtCallbackProc) green_slider_moved, NULL);
    XtAddCallback(green_slider, XmNdragCallback, (XtCallbackProc) green_slider_moved, NULL);

    blue_slider = XtVaCreateManagedWidget("Blue       ", xmScaleWidgetClass, s1, XmNminimum, 0, XmNmaximum, MAXCOLORVAL, XmNorientation, XmHORIZONTAL, XmNscaleWidth, 255, NULL);
    XtAddCallback(blue_slider, XmNvalueChangedCallback, (XtCallbackProc) blue_slider_moved, NULL);
    XtAddCallback(blue_slider, XmNdragCallback, (XtCallbackProc) blue_slider_moved, NULL);

    h_slider = XtVaCreateManagedWidget("Hue        ", xmScaleWidgetClass, s2, XmNminimum, 0, XmNmaximum, 360, XmNorientation, XmHORIZONTAL, XmNscaleWidth, 360, NULL);
    XtAddCallback(h_slider, XmNvalueChangedCallback, (XtCallbackProc) h_slider_moved, NULL);
    XtAddCallback(h_slider, XmNdragCallback, (XtCallbackProc) h_slider_moved, NULL);

    s_slider = XtVaCreateManagedWidget("Saturation    ", xmScaleWidgetClass, s2, XmNminimum, 0, XmNmaximum, MAXCOLORVAL, XmNorientation, XmHORIZONTAL, XmNscaleWidth, 250, NULL);
    XtAddCallback(s_slider, XmNvalueChangedCallback, (XtCallbackProc) s_slider_moved, NULL);
    XtAddCallback(s_slider, XmNdragCallback, (XtCallbackProc) s_slider_moved, NULL);


    v_slider = XtVaCreateManagedWidget("Value      ", xmScaleWidgetClass, s2, XmNminimum, 0, XmNmaximum, MAXCOLORVAL, XmNorientation, XmHORIZONTAL, XmNscaleWidth, 250, NULL);
    XtAddCallback(v_slider, XmNvalueChangedCallback, (XtCallbackProc) v_slider_moved, NULL);
    XtAddCallback(v_slider, XmNdragCallback, (XtCallbackProc) v_slider_moved, NULL);

    fr = XmCreateFrame(cedit_panel, "frame", NULL, 0);
    s1 = XmCreateRowColumn(fr, "s1", NULL, 0);
    XtVaSetValues(s1, XmNorientation, XmVERTICAL, NULL);

    s2 = XmCreateRowColumn(s1, "s2", NULL, 0);
    XtVaSetValues(s2, XmNorientation, XmHORIZONTAL, NULL);

    rgbhls_fr1 = XmCreateFrame(s2, "RGBvHLS", NULL, 0);
    rgbhls_rc1 = XmCreateRadioBox(rgbhls_fr1, "RGBvHLS", NULL, 0);

    rgbhls_tog1[0] = XmCreateToggleButton(rgbhls_rc1, "Interpolate in RGB", NULL, 0);
    XtVaSetValues(rgbhls_tog1[0], XmNset, True, NULL);
    XtAddCallback(rgbhls_tog1[0], XmNvalueChangedCallback, (XtCallbackProc) setspace, (XtPointer) 1);
    rgbhls_tog1[1] = XmCreateToggleButton(rgbhls_rc1, "Interpolate in HSV", NULL, 0);
    XtAddCallback(rgbhls_tog1[1], XmNvalueChangedCallback, (XtCallbackProc) setspace, 0);
    XtManageChildren(rgbhls_tog1, 2);
    XtManageChild(rgbhls_rc1);
    XtManageChild(rgbhls_fr1);

    rc = XmCreateRowColumn(s2, "rc", NULL, 0);
    interpf = CreateTextItem2(rc, 5, "From color #");
    interpt = CreateTextItem2(rc, 5, "To color #");
    XtManageChild(rc);

    XtManageChild(s2);

    s2 = XmCreateRowColumn(s1, "s2", NULL, 0);
    XtVaSetValues(s2, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtCreateManagedWidget("Interpolate", xmPushButtonWidgetClass, s2, NULL, 0);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) interp_proc, 0);
    bt = XtCreateManagedWidget("Undo", xmPushButtonWidgetClass, s2, NULL, 0);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) restorecolors, 0);

    XtManageChild(s2);
    XtManageChild(s1);
    XtManageChild(fr);


    XtCreateManagedWidget("sep", xmSeparatorWidgetClass, cedit_panel, NULL, 0);

    s1 = XmCreateRowColumn(cedit_panel, "s1", NULL, 0);
    XtVaSetValues(s1, XmNorientation, XmHORIZONTAL, NULL);

    bt = XtCreateManagedWidget("Done", xmPushButtonWidgetClass, s1, NULL, 0);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) quit_proc, 0);
    bt = XtCreateManagedWidget("Read...", xmPushButtonWidgetClass, s1, NULL, 0);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) readcmap_proc, 0);
    bt = XtCreateManagedWidget("Write...", xmPushButtonWidgetClass, s1, NULL, 0);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) writecmap_proc, 0);
    bt = XtCreateManagedWidget("Reset", xmPushButtonWidgetClass, s1, NULL, 0);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) defaultcolors, 0);
    XtManageChild(s1);

    current_color.pixel = colors[0];
    current_color.red = r[0] << 8;
    current_color.green = g[0] << 8;
    current_color.blue = b[0] << 8;
    current_color.flags = DoRed | DoGreen | DoBlue;

    updatesliders();
    update_color();

    XtManageChild(cedit_panel);
    XtManageChild(cedit_frame);
}

static void create_color_bar(Widget parent)
{
    Widget panel1, panel2, panel3, pan;
    int i, n;
    char name[10];
    Arg wargs[10];
    wl = (WidgetList) XtMalloc(ncolors * sizeof(Widget));

    /*
     * Create the row column manager to hold all color buttons.
     */
    n = 0;
    XtSetArg(wargs[n], XmNpacking, XmPACK_COLUMN);
    n++;
    XtSetArg(wargs[n], XmNnumColumns, 1);
    n++;
    pan = XtCreateManagedWidget("colorpanel", xmRowColumnWidgetClass, parent, wargs, n);
    for (i = 0; i < ncolors; i++) {
	if (i == 0) {
	    n = 0;
	    XtSetArg(wargs[n], XmNpacking, XmPACK_TIGHT);
	    n++;
	    XtSetArg(wargs[n], XmNorientation, XmHORIZONTAL);
	    n++;
	    XtSetArg(wargs[n], XmNadjustLast, False);
	    n++;
	    panel1 = XtCreateWidget("colorpanel", xmRowColumnWidgetClass, pan, wargs, n);
	    XtManageChild(panel1);
	}
	if (i == 16) {
	    XtManageChildren(wl, 16);
	    n = 0;
	    XtSetArg(wargs[n], XmNpacking, XmPACK_TIGHT);
	    n++;
	    XtSetArg(wargs[n], XmNorientation, XmHORIZONTAL);
	    n++;
	    XtSetArg(wargs[n], XmNadjustLast, False);
	    n++;
	    panel2 = XtCreateWidget("colorpanel", xmRowColumnWidgetClass, pan, wargs, n);
	    XtManageChild(panel2);
	}
	if (i == 32) {
	    XtManageChildren(wl + 16, 16);
	    n = 0;
	    XtSetArg(wargs[n], XmNpacking, XmPACK_TIGHT);
	    n++;
	    XtSetArg(wargs[n], XmNorientation, XmHORIZONTAL);
	    n++;
	    XtSetArg(wargs[n], XmNadjustLast, False);
	    n++;
	    panel3 = XtCreateWidget("colorpanel", xmRowColumnWidgetClass, pan, wargs, n);
	    XtManageChild(panel3);
	}
	n = 0;
	XtSetArg(wargs[n], XtNbackground, colors[i]);
	n++;
	XtSetArg(wargs[n], XmNwidth, 26);
	n++;
	sprintf(name, "%d", i);
	if (i < 16) {
	    wl[i] = XtCreateWidget(name, xmLabelWidgetClass, panel1, wargs, n);
	} else if (i < 32) {
	    wl[i] = XtCreateWidget(name, xmLabelWidgetClass, panel2, wargs, n);
	} else if (i < 48) {
	    wl[i] = XtCreateWidget(name, xmLabelWidgetClass, panel3, wargs, n);
	}
	XtAddEventHandler(wl[i], ButtonPressMask, False, (XtEventHandler) set_current_color, (XtPointer) i);
    }
    XtCreateManagedWidget("Annotation", xmLabelWidgetClass, panel1, NULL, 0);
    XtCreateManagedWidget("Bathymetry", xmLabelWidgetClass, panel2, NULL, 0);
    XtCreateManagedWidget("Elev/conc", xmLabelWidgetClass, panel3, NULL, 0);
    XtManageChildren(wl + 32, 16);
    XtManageChild(pan);
}

static void h_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data)
{
    double r, g, b;
    huev = call_data->value;
    hsv_to_rgb(huev, satv, valv, &r, &g, &b);
    current_color.red = (int) (r * MAXCOLORVAL);
    current_color.green = (int) (g * MAXCOLORVAL);
    current_color.blue = (int) (b * MAXCOLORVAL);
    updatesliders();
    update_color();
}

static void s_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data)
{
    double r, g, b;
    satv = call_data->value / (double) (MAXCOLORVAL);
    hsv_to_rgb(huev, satv, valv, &r, &g, &b);
    current_color.red = (int) (r * MAXCOLORVAL);
    current_color.green = (int) (g * MAXCOLORVAL);
    current_color.blue = (int) (b * MAXCOLORVAL);
    updatesliders();
    update_color();
}

static void v_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data)
{
    double r, g, b;
    valv = call_data->value / (double) (MAXCOLORVAL);
    hsv_to_rgb(huev, satv, valv, &r, &g, &b);
    current_color.red = (int) (r * MAXCOLORVAL);
    current_color.green = (int) (g * MAXCOLORVAL);
    current_color.blue = (int) (b * MAXCOLORVAL);
    updatesliders();
    update_color();
}

static void red_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data)
{
    double r, g, b;
    current_color.red = call_data->value;
    r = current_color.red;
    g = current_color.green;
    b = current_color.blue;
    rgb_to_hsv(r, g, b, &huev, &satv, &valv);
    updatesliders();
    update_color();
}

static void blue_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data)
{
    double r, g, b;
    current_color.blue = call_data->value;
    r = current_color.red;
    g = current_color.green;
    b = current_color.blue;
    rgb_to_hsv(r, g, b, &huev, &satv, &valv);
    updatesliders();
    update_color();
}

static void green_slider_moved(Widget w, XtPointer client_data, XmScaleCallbackStruct * call_data)
{
    double r, g, b;
    current_color.green = call_data->value;
    r = current_color.red;
    g = current_color.green;
    b = current_color.blue;
    rgb_to_hsv(r, g, b, &huev, &satv, &valv);
    updatesliders();
    update_color();
}

static void update_color(void)
{
    Arg wargs[1];
    char str[256];
    XmString xmstr;
    extern int inwin;
    if (!inwin)
	return;

    sprintf(str, "Current color: %2d", ncur);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNlabelString, xmstr);
    XtSetValues(color_number, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Red: %3d", current_color.red >> 8);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(red_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Green: %3d", current_color.green >> 8);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(green_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Blue: %3d", current_color.blue >> 8);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(blue_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Hue: %3.0lf", huev);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(h_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Saturation: %.2lf", satv);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(s_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Value: %.2lf", valv);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(v_slider, wargs, 1);
    XmStringFree(xmstr);

    set_colormapdata(ncur, current_color.red >> 8, current_color.green >> 8, current_color.blue >> 8);
    if (DisplayPlanes(disp, DefaultScreen(disp)) > 8) {
	XtVaSetValues(wl[ncur], XtNbackground, colors[ncur], NULL);
    }
}

static void set_current_color(Widget w, int number, XEvent * event)
{
    Arg wargs[10];
    double r, g, b;
    int tmp;
    char str[256];
    XmString xmstr;

    current_color.flags = DoRed | DoGreen | DoBlue;
    current_color.pixel = colors[number];
    ncur = number;
    XQueryColor(disp, cmap, &current_color);
    r = current_color.red;
    g = current_color.green;
    b = current_color.blue;
    rgb_to_hsv(r, g, b, &huev, &satv, &valv);

    sprintf(str, "Current color: %2d", ncur);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNlabelString, xmstr);
    XtSetValues(color_number, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Red: %3d", current_color.red >> 8);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(red_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Green: %3d", current_color.green >> 8);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(green_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Blue: %3d", current_color.blue >> 8);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(blue_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Hue: %3.0lf", huev);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(h_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Saturation: %.2lf", satv);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(s_slider, wargs, 1);
    XmStringFree(xmstr);

    sprintf(str, "Value: %.2lf", valv);
    xmstr = XmStringLtoRCreate(str, XmSTRING_DEFAULT_CHARSET);
    XtSetArg(wargs[0], XmNtitleString, xmstr);
    XtSetValues(v_slider, wargs, 1);
    XmStringFree(xmstr);

    updatesliders();
}

static void updatesliders(void)
{
    Arg wargs[10];
    double h, s, v;
    double r, g, b;
    int tmp;

    r = current_color.red;
    g = current_color.green;
    b = current_color.blue;
    XtSetArg(wargs[0], XmNvalue, current_color.red);
    XtSetValues(red_slider, wargs, 1);

    XtSetArg(wargs[0], XmNvalue, current_color.green);
    XtSetValues(green_slider, wargs, 1);

    XtSetArg(wargs[0], XmNvalue, current_color.blue);
    XtSetValues(blue_slider, wargs, 1);

    rgb_to_hsv(r, g, b, &h, &s, &v);
    huev = h;
    satv = s;
    valv = v;
    tmp = (int) h;
    if (tmp < 0) {
	tmp = 0;
    } else if (tmp > 360) {
	tmp = 360;
    }
    if (h >= 360.0) {
	tmp = 360;
    }
    XtSetArg(wargs[0], XmNvalue, tmp);
    XtSetValues(h_slider, wargs, 1);

    tmp = (int) (s * MAXCOLORVAL);
    if (tmp < 0) {
	tmp = 0;
    } else if (tmp > MAXCOLORVAL) {
	tmp = MAXCOLORVAL;
    }
    XtSetArg(wargs[0], XmNvalue, tmp);
    XtSetValues(s_slider, wargs, 1);

    tmp = (int) (v * MAXCOLORVAL);
    if (tmp < 0) {
	tmp = 0;
    } else if (tmp > MAXCOLORVAL) {
	tmp = MAXCOLORVAL;
    }
    XtSetArg(wargs[0], XmNvalue, tmp);
    XtSetValues(v_slider, wargs, 1);
}
