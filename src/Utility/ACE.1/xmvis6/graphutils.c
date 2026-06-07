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
 * graph utilities
 *
 */

#ifndef lint
static char RCSid[] = "$Id: graphutils.c,v 1.4 2004/08/05 21:32:17 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

void push_world(void);
void update_all(int gno);

int graph_active(int k)
{
    return g[k].active == ON;
}

void get_graph_box(int i, boxtype * b)
{
    memcpy(b, &boxes[i], sizeof(boxtype));
}

void get_graph_line(int i, linetype * l)
{
    memcpy(l, &lines[i], sizeof(linetype));
}

void get_graph_string(int i, plotstr * s)
{
    memcpy(s, &pstr[i], sizeof(plotstr));
}

void get_graph_framep(int gno, framep * f)
{
    memcpy(f, &g[gno].f, sizeof(framep));
}

void get_graph_world(int gno, world * w)
{
    memcpy(w, &g[gno].w, sizeof(world));
}

void get_graph_view(int gno, view * v)
{
    memcpy(v, &g[gno].v, sizeof(view));
}

void get_graph_labels(int gno, labels * labs)
{
    memcpy(labs, &g[gno].labs, sizeof(labels));
}

void get_graph_tickmarks(int gno, tickmarks * t, int a)
{
    memcpy(t, &g[gno].t[a], sizeof(tickmarks));
}

/*
 *
 * set graph props
 *
 */

void set_graph_box(int i, boxtype * b)
{
    memcpy(&boxes[i], b, sizeof(boxtype));
}

void set_graph_line(int i, linetype * l)
{
    memcpy(&lines[i], l, sizeof(linetype));
}

void set_graph_string(int i, plotstr * s)
{
    memcpy(&pstr[i], s, sizeof(plotstr));
}

void set_graph_active(int gno)
{
    g[gno].active = ON;
}

void set_graph_framep(int gno, framep * f)
{
    memcpy(&g[gno].f, f, sizeof(framep));
}

void set_graph_world(int gno, world * w)
{
    memcpy(&g[gno].w, w, sizeof(world));
}

void set_graph_view(int gno, view * v)
{
    memcpy(&g[gno].v, v, sizeof(view));
}

void set_graph_labels(int gno, labels * labs)
{
    memcpy(&g[gno].labs, labs, sizeof(labels));
}

void set_graph_tickmarks(int gno, tickmarks * t, int a)
{
    memcpy(&g[gno].t[a], t, sizeof(tickmarks));
}

void set_axis_prop(int whichgraph, int naxis, int prop, double val)
{
    int i, j, startg, stopg;

    if (whichgraph == -1) {
	startg = 0;
	stopg = maxgraph - 1;
    } else {
	startg = whichgraph;
	stopg = whichgraph;
    }
    for (j = startg; j <= stopg; j++) {
	switch (prop) {
	case ON:
	case OFF:
	    switch (naxis) {
	    case 6:
		for (i = 0; i < 6; i++) {
		    g[j].t[i].active = (int) val;
		}
		break;
	    case 7:
		for (i = 0; i < 6; i += 2) {
		    g[j].t[i].active = (int) val;
		}
		break;
	    case 8:
		for (i = 1; i < 6; i += 2) {
		    g[j].t[i].active = (int) val;
		}
		break;
	    }
	    break;
	case COLOR:
	    switch (naxis) {
	    case 6:
		for (i = 0; i < 6; i++) {
		    g[j].t[i].tl_color = (int) val;
		    g[j].t[i].t_drawbarcolor = (int) val;
		    g[j].t[i].t_color = (int) val;
		    g[j].t[i].t_mcolor = (int) val;
		    g[j].t[i].label.color = (int) val;
		}
		break;
	    case 7:
		for (i = 0; i < 6; i += 2) {
		    g[j].t[i].tl_color = (int) val;
		    g[j].t[i].t_drawbarcolor = (int) val;
		    g[j].t[i].t_color = (int) val;
		    g[j].t[i].t_mcolor = (int) val;
		    g[j].t[i].label.color = (int) val;
		}
		break;
	    case 8:
		for (i = 1; i < 6; i += 2) {
		    g[j].t[i].tl_color = (int) val;
		    g[j].t[i].t_drawbarcolor = (int) val;
		    g[j].t[i].t_color = (int) val;
		    g[j].t[i].t_mcolor = (int) val;
		    g[j].t[i].label.color = (int) val;
		}
		break;
	    }
	    break;
	case LINEWIDTH:
	    switch (naxis) {
	    case 6:
		for (i = 0; i < 6; i++) {
		    g[j].t[i].tl_linew = (int) val;
		    g[j].t[i].t_linew = (int) val;
		    g[j].t[i].t_mlinew = (int) val;
		    g[j].t[i].t_drawbarlinew = (int) val;
		}
		break;
	    case 7:
		for (i = 0; i < 6; i += 2) {
		    g[j].t[i].tl_linew = (int) val;
		    g[j].t[i].t_linew = (int) val;
		    g[j].t[i].t_mlinew = (int) val;
		    g[j].t[i].t_drawbarlinew = (int) val;
		}
		break;
	    case 8:
		for (i = 1; i < 6; i += 2) {
		    g[j].t[i].tl_linew = (int) val;
		    g[j].t[i].t_linew = (int) val;
		    g[j].t[i].t_mlinew = (int) val;
		    g[j].t[i].t_drawbarlinew = (int) val;
		}
		break;
	    }
	    break;
	case FONTP:
	    switch (naxis) {
	    case 6:
		for (i = 0; i < 6; i++) {
		    g[j].t[i].tl_font = (int) val;
		    g[j].t[i].label.font = (int) val;
		}
		break;
	    case 7:
		for (i = 0; i < 6; i += 2) {
		    g[j].t[i].tl_font = (int) val;
		    g[j].t[i].label.font = (int) val;
		}
		break;
	    case 8:
		for (i = 1; i < 6; i += 2) {
		    g[j].t[i].tl_font = (int) val;
		    g[j].t[i].label.font = (int) val;
		}
		break;
	    }
	    break;
	case CHAR:
	    switch (naxis) {
	    case 6:
		for (i = 0; i < 6; i++) {
		    g[j].t[i].tl_charsize = val;
		    g[j].t[i].label.charsize = val;
		}
		break;
	    case 7:
		for (i = 0; i < 6; i += 2) {
		    g[j].t[i].tl_charsize = val;
		    g[j].t[i].label.charsize = val;
		}
		break;
	    case 8:
		for (i = 1; i < 6; i += 2) {
		    g[j].t[i].tl_charsize = val;
		    g[j].t[i].label.charsize = val;
		}
		break;
	    }
	    break;
	}
    }
}

int iscontained(int gno, double wx, double wy)
{
    int i;
    double xconv(), yconv();
    double x = xconv(wx), y = yconv(wy);

    for (i = 0; i < maxgraph; i++) {
	if (g[i].active == ON) {
	    if ((g[i].v.xv1 <= x && g[i].v.xv2 >= x) && (g[i].v.yv1 <= y && g[i].v.yv2 >= y)) {
		return i;
	    }
	}
    }
    return gno;
}

int islogx(int gno)
{
    return (g[gno].type == LOGX || g[gno].type == LOGXY);
}

int islogy(int gno)
{
    return (g[gno].type == LOGY || g[gno].type == LOGXY);
}

char *graph_types(int it, int which)
{
    static char s[128];

    strcpy(s, "UNKNOWN");
    switch (it) {
    case XYFIXED:
	if (which) {
	    strcpy(s, "xyfixed");
	} else {
	    strcpy(s, "XYFIXED");
	}
	break;
    case XY:
	if (which) {
	    strcpy(s, "xy");
	} else {
	    strcpy(s, "XY");
	}
	break;
    case LOGX:
	if (which) {
	    strcpy(s, "logx");
	} else {
	    strcpy(s, "LOG-LINEAR");
	}
	break;
    case LOGY:
	if (which) {
	    strcpy(s, "logy");
	} else {
	    strcpy(s, "LINEAR-LOG");
	}
	break;
    case LOGXY:
	if (which) {
	    strcpy(s, "logxy");
	} else {
	    strcpy(s, "LOG-LOG");
	}
	break;
    case POLAR:
	strcpy(s, "polar");
	break;
    case BAR:
	if (which) {
	    strcpy(s, "bar");
	} else {
	    strcpy(s, "BAR");
	}
	break;
    case HBAR:
	if (which) {
	    strcpy(s, "hbar");
	} else {
	    strcpy(s, "HORIZONTAL BAR");
	}
	break;
    case PIE:
	strcpy(s, "pie");
	break;
    case STACKEDBAR:
	if (which) {
	    strcpy(s, "stackedbar");
	} else {
	    strcpy(s, "STACKED BAR");
	}
	break;
    case STACKEDHBAR:
	if (which) {
	    strcpy(s, "stackedhbar");
	} else {
	    strcpy(s, "STACKED HORIZONTAL BAR");
	}
	break;
    case STACKEDLINE:
	strcpy(s, "STACKED LINE");
	break;
    }
    return s;
}

int get_format_index(int f)
{
    int i = 0;

    while (f != format_types[i] && format_types[i] != 0)
	i++;
    return i;
}

char *get_format_types(int f)
{
    static char s[128];

    strcpy(s, "decimal");
    switch (f) {
    case DECIMAL:
	strcpy(s, "decimal");
	break;
    case EXPONENTIAL:
	strcpy(s, "exponential");
	break;
    case POWER:
	strcpy(s, "power");
	break;
    case GENERAL:
	strcpy(s, "general");
	break;
    case DDMMYY:
	strcpy(s, "ddmmyy");
	break;
    case MMDDYY:
	strcpy(s, "mmddyy");
	break;
    case MMYY:
	strcpy(s, "mmyy");
	break;
    case MMDD:
	strcpy(s, "mmdd");
	break;
    case MONTHDAY:
	strcpy(s, "monthday");
	break;
    case DAYMONTH:
	strcpy(s, "daymonth");
	break;
    case MONTHS:
	strcpy(s, "months");
	break;
    case MONTHL:
	strcpy(s, "monthl");
	break;
    case DAYOFWEEKS:
	strcpy(s, "dayofweeks");
	break;
    case DAYOFWEEKL:
	strcpy(s, "dayofweekl");
	break;
    case DAYOFYEAR:
	strcpy(s, "dayofyear");
	break;
    case HMS:
	strcpy(s, "hms");
	break;
    case MMDDHMS:
	strcpy(s, "mmddhms");
	break;
    case MMDDYYHMS:
	strcpy(s, "mmddyyhms");
	break;
    case DEGREESLON:
	strcpy(s, "degreeslon");
	break;
    case DEGREESMMLON:
	strcpy(s, "degreesmmlon");
	break;
    case DEGREESMMSSLON:
	strcpy(s, "degreesmmsslon");
	break;
    case MMSSLON:
	strcpy(s, "mmsslon");
	break;
    case DEGREESLAT:
	strcpy(s, "degreeslat");
	break;
    case DEGREESMMLAT:
	strcpy(s, "degreesmmlat");
	break;
    case DEGREESMMSSLAT:
	strcpy(s, "degreesmmsslat");
	break;
    case MMSSLAT:
	strcpy(s, "mmsslat");
	break;
    }
    return s;
}

void kill_graph(int gno)
{
    int i, j;

    if (gno == maxgraph) {
	for (i = 0; i < maxgraph; i++) {
	    set_default_graph(i);
	}
    } else {
	set_default_graph(gno);
    }
}

void copy_graph(int from, int to)
{
    int i, j;
    kill_graph(to);
    memcpy(&g[to], &g[from], sizeof(graph));
    set_graph_active(to);	/* TODO compare maxplots */
}

void swap_graph(int from, int to)
{
    graph gtmp;
    memcpy(&gtmp, &g[from], sizeof(graph));
    memcpy(&g[from], &g[to], sizeof(graph));
    memcpy(&g[to], &gtmp, sizeof(graph));
}

/*
 * for flipping
 */
void flipxy(int gno)
{
    int i, j;
    tickmarks t;
    double *x, *y;

    for (i = 0; i < 6; i += 2) {
	memcpy(&t, &g[gno].t[i], sizeof(tickmarks));
	memcpy(&g[gno].t[i], &g[gno].t[i + 1], sizeof(tickmarks));
	memcpy(&g[gno].t[i + 1], &t, sizeof(tickmarks));
	if (g[gno].t[i].t_op == RIGHT) {
	    g[gno].t[i].t_op = TOP;
	} else if (g[gno].t[i].t_op == LEFT) {
	    g[gno].t[i].t_op = BOTTOM;
	}
	if (g[gno].t[i].tl_op == RIGHT) {
	    g[gno].t[i].tl_op = TOP;
	} else if (g[gno].t[i].tl_op == LEFT) {
	    g[gno].t[i].tl_op = BOTTOM;
	}
	if (g[gno].t[i + 1].t_op == TOP) {
	    g[gno].t[i + 1].t_op = RIGHT;
	} else if (g[gno].t[i + 1].t_op == BOTTOM) {
	    g[gno].t[i + 1].t_op = LEFT;
	}
	if (g[gno].t[i + 1].tl_op == TOP) {
	    g[gno].t[i + 1].tl_op = RIGHT;
	} else if (g[gno].t[i + 1].tl_op == BOTTOM) {
	    g[gno].t[i + 1].tl_op = LEFT;
	}
    }
    if (g[gno].type == LOGX) {
	g[gno].type = LOGY;
    } else if (g[gno].type == LOGY) {
	g[gno].type = LOGX;
    }
    fswap(&g[gno].w.xg1, &g[gno].w.yg1);
    fswap(&g[gno].w.xg2, &g[gno].w.yg2);
    fswap(&g[gno].dsx, &g[gno].dsy);
    iswap(&g[gno].fx, &g[gno].fy);
    iswap(&g[gno].px, &g[gno].py);
    update_all(gno);
    drawgraph();
}

void invertx(int gno)
{
    int i, j;
    double *x;

    if (!islogx(gno)) {
	fswap(&g[gno].w.xg1, &g[gno].w.xg2);
	g[gno].w.xg1 = -g[gno].w.xg1;
	g[gno].w.xg2 = -g[gno].w.xg2;
	g[gno].dsx = -g[gno].dsx;
	for (i = 0; i < 6; i += 2) {
	    if (g[gno].t[i].tl_sign == NORMAL) {
		g[gno].t[i].tl_sign = NEGATE;
	    } else {
		g[gno].t[i].tl_sign = NORMAL;
	    }
	}
	update_all(gno);
	drawgraph();
    } else {
	errwin("Can't invert log axes");
    }
}

void inverty(int gno)
{
    int i, j;
    double *y;

    if (!islogy(gno)) {
	fswap(&g[gno].w.yg1, &g[gno].w.yg2);
	g[gno].w.yg1 = -g[gno].w.yg1;
	g[gno].w.yg2 = -g[gno].w.yg2;
	g[gno].dsy = -g[gno].dsy;
	for (i = 1; i < 6; i += 2) {
	    if (g[gno].t[i].tl_sign == NORMAL) {
		g[gno].t[i].tl_sign = NEGATE;
	    } else {
		g[gno].t[i].tl_sign = NORMAL;
	    }
	}
	update_all(gno);
	drawgraph();
    } else {
	errwin("Can't invert log axes");
    }
}

void do_flipxy(void)
{
    flipxy(cg);
}


void do_invertx(void)
{
    invertx(cg);
}

void do_inverty(void)
{
    inverty(cg);
}

/*
	 the following routines determine default scaling and tickmarks
*/
void defaultgraph(int gno)
{
    double x1, x2, y1, y2;
    double xgmax, xgmin, ygmax, ygmin;
    int i, first = 1;

    xgmax = xgmin = ygmax = ygmin = 0.0;
    for (i = 0; i < g[gno].maxplot; i++) {
	if (first) {
	    xgmin = x1;
	    xgmax = x2;
	    ygmin = y1;
	    ygmax = y2;
	    first = 0;
	} else {
	    xgmin = (x1 < xgmin) ? x1 : xgmin;
	    xgmax = (x2 > xgmax) ? x2 : xgmax;
	    ygmin = (y1 < ygmin) ? y1 : ygmin;
	    ygmax = (y2 > ygmax) ? y2 : ygmax;
	}
    }
    if (xgmin != xgmax) {
	g[gno].w.xg2 = xgmax;
	g[gno].w.xg1 = xgmin;
    } else {
	g[gno].w.xg1 = xgmin - 1.0;
	g[gno].w.xg2 = xgmin + 1.0;
    }
    if (ygmin != ygmax) {
	g[gno].w.yg2 = ygmax;
	g[gno].w.yg1 = ygmin;
    } else {
	g[gno].w.yg1 = ygmin - 1.0;
	g[gno].w.yg2 = ygmin + 1.0;
    }
    switch (g[gno].type) {
    case LOGX:
	if (g[gno].w.xg1 <= 0.0) {
	    errwin("can't set graph type to log-linear, X minimum = 0.0");
	    g[gno].type = XY;
	}
	break;
    case LOGY:
	if (g[gno].w.yg1 <= 0.0) {
	    errwin("can't set graph type to linear-log, Y minimum = 0.0");
	    g[gno].type = XY;
	}
	break;
    case LOGXY:
	if (g[gno].w.xg1 <= 0.0) {
	    errwin("can't set graph to log-log, X minimum <= 0.0");
	    g[gno].type = XY;
	} else if (g[gno].w.yg1 <= 0.0) {
	    errwin("can't set graph type to log-log, Y minimum <= 0.0");
	    g[gno].type = XY;
	}
	break;
    }
}

void defaultx(int gno, int setno)
{
    int i, first = 1;
    double xgmin, xgmax, xmax, xmin, tmp;

    xgmin = xgmax = 0.0;
    if (setno < 0) {
	for (i = 0; i < g[gno].maxplot; i++) {
	    if (first) {
		xgmin = xmin;
		xgmax = xmax;
		first = 0;
	    } else {
		xgmin = (xmin < xgmin) ? xmin : xgmin;
		xgmax = (xmax > xgmax) ? xmax : xgmax;
	    }
	}
    } else {
    }
    if (xgmin != xgmax) {
	g[gno].w.xg2 = xgmax;
	g[gno].w.xg1 = xgmin;
    } else {
	if ((xgmin == 0.0) && (xgmax == 0.0)) {
	    xgmin = 1.0;
	}
	g[gno].w.xg1 = xgmin - 0.1 * fabs(xgmin);
	g[gno].w.xg2 = xgmin + 0.1 * fabs(xgmin);
    }
    switch (g[gno].type) {
    case LOGX:
	if (g[gno].w.xg1 <= 0.0) {
	    errwin("can't set graph type to log-linear, X minimum = 0.0");
	    g[gno].type = XY;
	}
	break;
    case LOGY:
	if (g[gno].w.yg1 <= 0.0) {
	    errwin("can't set graph type to linear-log, Y minimum = 0.0");
	    g[gno].type = XY;
	}
	break;
    case LOGXY:
	if (g[gno].w.xg1 <= 0.0) {
	    errwin("can't set graph to log-log, X minimum <= 0.0");
	    g[gno].type = XY;
	} else if (g[gno].w.yg1 <= 0.0) {
	    errwin("can't set graph type to log-log, Y minimum <= 0.0");
	    g[gno].type = XY;
	}
	break;
    }
}

void defaulty(int gno, int setno)
{
    int i, first = 1;
    double ygmax, ygmin, ymin, ymax, tmp;

    ygmin = ygmax = 0.0;
    if (setno < 0) {
	for (i = 0; i < g[gno].maxplot; i++) {
	    if (first) {
		ygmin = ymin;
		ygmax = ymax;
		first = 0;
	    } else {
		ygmin = (ymin < ygmin) ? ymin : ygmin;
		ygmax = (ymax > ygmax) ? ymax : ygmax;
	    }
	}
    } else {
    }
    if (ygmin != ygmax) {
	g[gno].w.yg2 = ygmax;
	g[gno].w.yg1 = ygmin;
    } else {
	if ((ygmin == 0.0) && (ygmax == 0.0)) {
	    ygmin = 1.0;
	}
	g[gno].w.yg1 = ygmin - 0.1 * fabs(ygmin);
	g[gno].w.yg2 = ygmin + 0.1 * fabs(ygmin);
    }
    switch (g[gno].type) {
    case LOGX:
	if (g[gno].w.xg1 <= 0.0) {
	    errwin("can't set graph type to log-linear, X minimum = 0.0");
	    g[gno].type = XY;
	}
	break;
    case LOGY:
	if (g[gno].w.yg1 <= 0.0) {
	    errwin("can't set graph type to linear-log, Y minimum = 0.0");
	    g[gno].type = XY;
	}
	break;
    case LOGXY:
	if (g[gno].w.xg1 <= 0.0) {
	    errwin("can't set graph to log-log, X minimum <= 0.0");
	    g[gno].type = XY;
	} else if (g[gno].w.yg1 <= 0.0) {
	    errwin("can't set graph type to log-log, Y minimum <= 0.0");
	    g[gno].type = XY;
	}
	break;
    }
}

/*
 * label: test program to demonstrate nice graph axis labeling
 *
 * Paul Heckbert, 2 Dec 88
 */

double expt(double a, register int n), nicenum(double x, int round);

void default_ticks(int gno, int axis, double *gmin, double *gmax)
{
    tickmarks t;
    double range, d, tmpmax = *gmax, tmpmin = *gmin;

    get_graph_tickmarks(gno, &t, axis);
    if (axis % 2 && (g[gno].type == LOGY || g[gno].type == LOGXY)) {
	tmpmax = ceil(log10(tmpmax));
	tmpmin = floor(log10(tmpmin));
    } else if ((axis % 2 == 0) && (g[gno].type == LOGX || g[gno].type == LOGXY)) {
	tmpmax = ceil(log10(tmpmax));
	tmpmin = floor(log10(tmpmin));
    }
    range = nicenum(tmpmax - tmpmin, 0);
    d = nicenum(range / (t.t_num - 1), 1);
    tmpmin = floor(tmpmin / d) * d;
    tmpmax = ceil(tmpmax / d) * d;
    if (axis % 2 && (g[gno].type == LOGY || g[gno].type == LOGXY)) {
	*gmax = pow(10.0, tmpmax);
	*gmin = pow(10.0, tmpmin);
	t.tmajor = (int) d;
	if (t.tmajor == 0.0) {
	    t.tmajor = 1.0;
	}
	if ((int) t.tmajor < 2) {
	    t.tminor = 1.0;
	} else {
	    t.tminor = 0.0;
	}
	if (fabs(tmpmax) > 6.0 || fabs(tmpmin) > 6.0) {
	    t.tl_format = POWER;
	    t.tl_prec = 0;
	} else {
	    t.tl_format = DECIMAL;
	    t.tl_prec = 0;
	}
    } else if ((axis % 2 == 0) && (g[gno].type == LOGX || g[gno].type == LOGXY)) {
	*gmax = pow(10.0, tmpmax);
	*gmin = pow(10.0, tmpmin);
	t.tmajor = (int) d;
	if (t.tmajor == 0.0) {
	    t.tmajor = 1.0;
	}
	if (fabs(tmpmax) > 6.0 || fabs(tmpmin) > 6.0) {
	    t.tl_format = POWER;
	    t.tl_prec = 0;
	} else {
	    t.tl_format = DECIMAL;
	    t.tl_prec = 0;
	}
	if ((int) t.tmajor < 2) {
	    t.tminor = 1.0;
	} else {
	    t.tminor = 0.0;
	}
    } else {
	*gmax = tmpmax;
	*gmin = tmpmin;
	t.tmajor = d;
	t.tminor = d * 0.5;
    }
    set_graph_tickmarks(gno, &t, axis);
}

/*
 * nicenum: find a "nice" number approximately equal to x
 * round if round=1, ceil if round=0
 */

double nicenum(double x, int round)
{
    int exp;
    double f, y;

    exp = floor(log10(x));
    f = x / expt(10., exp);	/* fraction between 1 and 10 */
    if (round)
	if (f < 1.5)
	    y = 1.;
	else if (f < 3.)
	    y = 2.;
	else if (f < 7.)
	    y = 5.;
	else
	    y = 10.;
    else if (f <= 1.)
	y = 1.;
    else if (f <= 2.)
	y = 2.;
    else if (f <= 5.)
	y = 5.;
    else
	y = 10.;
    return y * expt(10., exp);
}

/*
 * expt(a,n)=a^n for integer n
 * roundoff errors in pow were causing problems, so I wrote my own
 */

double expt(double a, register int n)
{
    double x;

    x = 1.;
    if (n > 0)
	for (; n > 0; n--)
	    x *= a;
    else
	for (; n < 0; n++)
	    x /= a;
    return x;
}

void defaultsetgraph(int gno, int setno)
{
    double xmax, xmin, ymax, ymin;

/*
    getsetminmax(gno, setno, &xmin, &xmax, &ymin, &ymax);
*/
/* TODO need to hook up to data in xmvis */
    xmin = ymin = 0.0;
    xmax = ymax = 1.0;

    if (xmin != xmax) {
	g[gno].w.xg2 = xmax;
	g[gno].w.xg1 = xmin;
    } else {
	if ((xmin == 0.0) && (xmax == 0.0)) {
	    xmin = 1.0;
	}
	g[gno].w.xg1 = xmin - 0.1 * fabs(xmin);
	g[gno].w.xg2 = xmin + 0.1 * fabs(xmin);
    }
    if (ymin != ymax) {
	g[gno].w.yg2 = ymax;
	g[gno].w.yg1 = ymin;
    } else {
	if ((ymin == 0.0) && (ymax == 0.0)) {
	    ymin = 1.0;
	}
	g[gno].w.yg1 = ymin - 0.1 * fabs(ymin);
	g[gno].w.yg2 = ymin + 0.1 * fabs(ymin);
    }
}

void default_axis(int gno, int method, int axis)
{
    int cx;
    tickmarks t;
    world w;
    double llim, ulim;

    get_graph_tickmarks(gno, &t, axis);
    get_graph_world(gno, &w);
    if (axis % 2) {
	llim = w.yg1;
	ulim = w.yg2;
    } else {
	llim = w.xg1;
	ulim = w.xg2;
    }
    t.tmajor = (ulim - llim) / t.t_num;
    t.tminor = t.tmajor / 2.0;
    cx = (int) log10(t.tmajor);
    t.tl_prec = ((cx < 0) ? -cx + 1 : 1);
    if (t.tl_prec > 9) {
	t.tl_prec = 2;
	t.tl_format = EXPONENTIAL;
    }
    set_graph_tickmarks(gno, &t, axis);
    if (method == AUTO) {
	default_ticks(gno, axis, &llim, &ulim);
	if (axis % 2) {
	    w.yg1 = llim;
	    w.yg2 = ulim;
	} else {
	    w.xg1 = llim;
	    w.xg2 = ulim;
	}
	set_graph_world(gno, &w);
    }
    if (g[gno].type == XYFIXED) {
	setfixedscale(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2, &g[gno].w.xg1, &g[gno].w.yg1, &g[gno].w.xg2, &g[gno].w.yg2);
    }
}

void autoscale_graph(int gno, int axis)
{
    if (activeset(gno)) {
	switch (axis) {
	case -3:
	    defaultgraph(gno);
	    default_axis(gno, g[gno].auto_type, X_AXIS);
	    default_axis(gno, g[gno].auto_type, ZX_AXIS);
	    default_axis(gno, g[gno].auto_type, Y_AXIS);
	    default_axis(gno, g[gno].auto_type, ZY_AXIS);
	    break;
	case -2:
	    defaultx(gno, -1);
	    default_axis(gno, g[gno].auto_type, X_AXIS);
	    default_axis(gno, g[gno].auto_type, ZX_AXIS);
	    break;
	case -1:
	    defaulty(gno, -1);
	    default_axis(gno, g[gno].auto_type, Y_AXIS);
	    default_axis(gno, g[gno].auto_type, ZY_AXIS);
	    break;
	default:
	    if (axis % 2) {
		defaulty(gno, -1);
	    } else {
		defaultx(gno, -1);
	    }
	    default_axis(gno, g[gno].auto_type, axis);
	    break;
	}
	update_all(gno);
    }
}

void autoscale_proc(void)
{
    if (activeset(cg)) {
	defaultgraph(cg);
	default_axis(cg, g[cg].auto_type, X_AXIS);
	default_axis(cg, g[cg].auto_type, ZX_AXIS);
	default_axis(cg, g[cg].auto_type, Y_AXIS);
	default_axis(cg, g[cg].auto_type, ZY_AXIS);

	update_all(cg);
	drawgraph();
    } else {
	errwin("No active sets!");
    }
}

void autoticks_proc(void)
{
    default_axis(cg, g[cg].auto_type, X_AXIS);
    default_axis(cg, g[cg].auto_type, ZX_AXIS);
    default_axis(cg, g[cg].auto_type, Y_AXIS);
    default_axis(cg, g[cg].auto_type, ZY_AXIS);
    update_all(cg);
    drawgraph();
}

void autoscale_set(int gno, int setno, int axis)
{
}

void wipeout(int ask)
{
    if (ask && !yesno("Kill all graphs, sets, and annotation?", "", "YES", "NO")) {
	return;
    }
    kill_graph(maxgraph);
    do_clear_lines();
    do_clear_boxes();
    do_clear_text();
    cg = 0;
    drawgraph();
}

void update_all(int gno)
{
    set_display_pixmaps();	/* also calls update_display_items() */
    update_teanl_flow();
    update_adcirc_flow();
    update_draw();
    update_graph_items();
    update_world(gno);
    update_view(gno);
/*
    update_status(gno, cur_statusitem, -1);
    update_label_proc();
    update_editp_proc(gno, -1);
*/
    update_ticks(gno);
    update_autos(gno);
    update_locator_items(gno);
    update_draw();
    update_frame_items(gno);
    update_graph_items();
}

int object_isactive(int obj, int w, int w2)
{
    switch (obj) {
    case GRID:
	return grid[w].active == ON;
    case BOUNDARY:
	return bounds[w].active == ON;
    case TEANL:
	return flowf[w].active == ON;
    case ADCIRC:
	if (w == FLOW) {
	    return ((flowt[w2].active == ON) && (flowt[w2].f != NULL) && (flowt[w2].f[0].u != NULL));
	} else if (w == ELEV) {
	    return (flowt[w2].active == ON && flowt[w2].f != NULL && flowt[w2].f[0].e != NULL);
	}
    case TRANSECT:
	return trans[w].active == ON;
    case ADCIRC3DFLOW:
	return adc3d[w].active == ON;
    case HISTORY:
	if (w == FLOW) {
	    return flowh[w2].active == ON;
	} else {
	    return hist[w2].active == ON;
	}
    case ELA:
	return elaconc[w].active == ON;
    case STATION:
	if (nsta > 0 && sta != NULL && w < nsta) {
	    return sta[w].active == ON;
	} else {
	    return 0;
	}
    case TIDESTATION:
	if (tidestat == NULL || tidestat[w] == NULL) {
	    return 0;
	}
	return 1;
    case DROGUES:
	return drogues[w].active == ON;
    case TRACK:
	return track[w].active == ON;
    default:
	fprintf(stderr, "Object not found: %d\n", obj);
	return 0;
    }
}

int display_object(int gno, int obj, int w1, int w2, int w3)
{
    switch (obj) {
    case GRID:
	if (w1 == GRID) {
	    return g[gno].grid[w2].display != OFF;
	}
	if (w1 == BOUNDARY) {
	    return g[gno].grid[w2].display_boundary != OFF;
	}
	break;
    case BATH:
	return g[gno].grid[w1].display_bath != OFF;
    case BOUNDARY:
	return g[gno].bounds[w1].display != OFF;
    case TEANL:
	if (w1 == ELEV) {
	    return g[gno].flowf[w2].display_elev != OFF;
	}
	if (w1 == FLOW) {
	    return g[gno].flowf[w2].display != OFF;
	}
	if (w1 == PHASE) {
	    return g[gno].flowf[w2].display_phase != OFF;
	}
	if (w1 == AMP) {
	    return g[gno].flowf[w2].display_amp != OFF;
	}
	if (w1 == MAG) {
	    return g[gno].flowf[w2].display_mag != OFF;
	}
	if (w1 == MARKERS) {
	    return g[gno].flowf[w2].display_elevmarkers != OFF;
	}
	if (w1 == WETDRY) {
	    return g[gno].flowf[w2].display_inun != OFF;
	}
	break;
    case ADCIRC:
	if (w1 == ELEV && w2 == MAXP) {
	    return g[gno].flowt[w3].display_maxelev != OFF;
	}
	if (w1 == ELEV && w2 == NODES) {
	    return g[gno].flowt[w3].display_maxelevval != OFF;
	}
	if (w1 == ELEV) {
	    return g[gno].flowt[w2].display_elev != OFF;
	}
	if (w1 == FLOW) {
	    return g[gno].flowt[w2].display != OFF;
	}
	if (w1 == MAG) {
	    return g[gno].flowt[w2].display_mag != OFF;
	}
	if (w1 == MARKERS) {
	    return g[gno].flowt[w2].display_elevmarkers != OFF;
	}
	break;
    case ADCIRC3DFLOW:
	if (w1 == FLOW) {
	    return g[gno].flow3d[w2].display != OFF;
	}
	break;
    case HISTORY:
	if (w1 == ELEV) {
	    return g[gno].flowh[w2].display_elev != OFF;
	} else if (w1 == FLOW) {
	    return g[gno].flowh[w2].display != OFF;
	} else if (w1 == MARKERS) {
	    return g[gno].flowh[w2].display_elevmarkers != OFF;
	} else {
	    return g[gno].hbox[w2].display != OFF;
	}
    case TRANSECT:
	return g[gno].trans[w1].display != OFF;
    case ELA:
	return g[gno].elaconc[w1].display != OFF;
    case STATION:
	if (nsta > 0 && sta != NULL && w1 < nsta) {
	    return sta[w1].display != OFF;
	} else {
	    return 0;
	}
    case TIDESTATION:
	if (tidestat == NULL || tidestat[w1] == NULL || g[gno].tidestat == NULL) {
	    return 0;
	}
	return g[gno].tidestat[w1].display != OFF;
    case DROGUES:
	return g[gno].drogues[w1].display != OFF;
    case TRACK:
	return g[gno].track[w1].display != OFF;
    case SLICE:
	return g[gno].sbox[w1].display_slice != OFF || g[gno].sbox[w1].display_marker != OFF;
    case FLUX:
	return g[gno].fbox[w1].display_slice != OFF || g[gno].fbox[w1].display_marker != OFF;
    case ZOOM:
	return g[gno].zbox[w1].display != OFF || g[gno].zbox[w1].display_marker != OFF;
    case TIDALCLOCK:
	return g[gno].tidalclock.active == ON;
    case TIMEINFO:
	return g[gno].timeinfo.active == ON;
    case MAPSCALE:
	return g[gno].mapscale.active == ON;
    case VSCALE:
	return g[gno].vl.active == ON;
    case WSCALE:
	return g[gno].wl.active == ON;
    case FLOW:
	return g[gno].fl.active == ON;
    case TIMELINE:
	return g[gno].timeline.active == ON;
    default:
	break;
    }
    return 0;
}

void set_display_object(int gno, int obj, int w1, int w2, int w3)
{
/*
        g[gno].grid[i].display_boundary = OFF;
        g[gno].grid[i].display_nodes = OFF;
        g[gno].grid[i].display_elements = OFF;
        g[gno].grid[i].display_depths = OFF;
        g[gno].grid[i].display_courant = OFF;
        g[gno].grid[i].display_dimw = OFF;
        g[gno].grid[i].display_gridf = OFF;
        g[gno].grid[i].display_flags[EDIT_GRID_FILLED] = 0;
        g[gno].grid[i].display_flags[EDIT_GRID_NODE_NUMBERS] = 0;
        g[gno].grid[i].display_flags[EDIT_GRID_ELEMENT_NUMBERS] = 0;
        g[gno].grid[i].display_flags[EDIT_GRID_DEPTHS] = 0;
    }
    g[gno].grid[0].display_flags[EDIT_GRID_ISOLINES] = 1;
*/

    switch (obj) {
    case GRID:
	g[gno].grid[w1].display = w2;
	g[gno].grid[w1].display_flags[EDIT_GRID] = (w2 == ON);
	break;
    case BATH:
	g[gno].grid[w1].display_bath = w2;
	g[gno].grid[w1].display_flags[EDIT_GRID_ISOLINES] = (w2 == ON);
	break;
    case BOUNDARY:
	g[gno].grid[w1].display_boundary = w2;
	g[gno].grid[w1].display_flags[EDIT_BOUNDARY] = (w2 == ON);
	break;
    }
/*
    case TEANL:
	if (w1 == ELEV) {
	    return g[gno].flowf[w2].display_elev != OFF;
	}
	if (w1 == FLOW) {
	    return g[gno].flowf[w2].display != OFF;
	}
	if (w1 == PHASE) {
	    return g[gno].flowf[w2].display_phase != OFF;
	}
	if (w1 == AMP) {
	    return g[gno].flowf[w2].display_amp != OFF;
	}
	if (w1 == MARKERS) {
	    return g[gno].flowf[w2].display_elevmarkers != OFF;
	}
	break;
    case ADCIRC:
	if (w1 == ELEV && w2 == MAXP) {
	    return g[gno].flowt[w3].display_maxelev != OFF;
	}
	if (w1 == ELEV) {
	    return g[gno].flowt[w2].display_elev != OFF;
	}
	if (w1 == FLOW) {
	    return g[gno].flowt[w2].display != OFF;
	}
	if (w1 == MARKERS) {
	    return g[gno].flowt[w2].display_elevmarkers != OFF;
	}
	break;
    case HISTORY:
	if (w1 == ELEV) {
	    return g[gno].flowh[w2].display_elev != OFF;
	}
	if (w1 == FLOW) {
	    return g[gno].flowh[w2].display != OFF;
	}
	if (w1 == MARKERS) {
	    return g[gno].flowh[w2].display_elevmarkers != OFF;
	}
	break;
    case ELA:
	return g[gno].elaconc[w1].display != OFF;
	break;
    case STATION:
	return stations[w1].display != OFF;
	break;
    case FISH:
	return fishes[w1].display != OFF;
	break;
    case DROGUES:
	return g[gno].drogues[w1].display != OFF;
	break;
    case SLICE:
	return g[gno].sbox[w1].display_slice != OFF ||
		g[gno].sbox[w1].display_marker != OFF;
	break;
    case FLUX:
	return g[gno].fbox[w1].display_slice != OFF ||
		g[gno].fbox[w1].display_marker != OFF;
	break;
    case ZOOM:
	return g[gno].zbox[w1].display != OFF ||
		g[gno].zbox[w1].display_marker != OFF;
	break;
    case TIDALCLOCK:
	return g[gno].tidalclock.active == ON;
	break;
    case MAPSCALE:
	return g[gno].mapscale.active == ON;
	break;
    case VSCALE:
	return g[gno].vl.active == ON;
	break;
    case TIMELINE:
	return g[gno].timeline.active == ON;
	break;
    default:
	return 0;
	break;
    }
*/
}


void page(double p, int dir)
{
    int restart = 0;
    double dx = scrollper * (g[cg].w.xg2 - g[cg].w.xg1);
    double dy = scrollper * (g[cg].w.yg2 - g[cg].w.yg1);
    if (isrunning()) {
	restart = 1;
	setistop();
    }
    switch (dir) {
    case 0:
	g[cg].w.xg1 += dx;
	g[cg].w.xg2 += dx;
	break;
    case 1:
	g[cg].w.xg1 -= dx;
	g[cg].w.xg2 -= dx;
	break;
    case 2:
	g[cg].w.yg1 += dy;
	g[cg].w.yg2 += dy;
	break;
    case 3:
	g[cg].w.yg1 -= dy;
	g[cg].w.yg2 -= dy;
	break;
    case 4:
	dx = shexper * (g[cg].w.xg2 - g[cg].w.xg1) / 2.0;
	dy = shexper * (g[cg].w.yg2 - g[cg].w.yg1) / 2.0;
	g[cg].w.xg1 += dx;
	g[cg].w.xg2 -= dx;
	g[cg].w.yg1 += dy;
	g[cg].w.yg2 -= dy;
	break;
    case 5:
	dx = shexper * (g[cg].w.xg2 - g[cg].w.xg1);
	dy = shexper * (g[cg].w.yg2 - g[cg].w.yg1);
	g[cg].w.xg1 -= dx;
	g[cg].w.xg2 += dx;
	g[cg].w.yg1 -= dy;
	g[cg].w.yg2 += dy;
	break;
    }
    set_up_world(cg);
    if (!no_display) {
	setredraw_world();
    }
    if (restart) {
	setirun();
    }
}

int check_err = 0;

#define MAX_FONT 10
#define MAX_JUST 2
#define MAX_ARROW 3
#define MAX_PATTERN 30
#define MAX_PREC 10

int checkon(int prop, int old_val, int new_val)
{
    char buf[256];
    int retval = old_val;
    check_err = 0;
    switch (prop) {
    case LINEWIDTH:
	if (new_val >= 0 && new_val <= MAX_LINEWIDTH) {
	    retval = new_val;
	} else {
	    sprintf(buf, "LINEWIDTH out of bounds, should be from 0 to %d", MAX_LINEWIDTH);
	    check_err = 1;
	}
	break;
    case LINESTYLE:
	if (new_val >= 0 && new_val <= MAX_LINESTYLE) {
	    retval = new_val;
	} else {
	    sprintf(buf, "LINESTYLE out of bounds, should be from 0 to %d", MAX_LINESTYLE);
	    check_err = 1;
	}
	break;
    case COLOR:		/* TODO use MAX_COLOR */
	if (new_val >= 0 && new_val < 16) {
	    retval = new_val;
	} else {
	    sprintf(buf, "COLOR out of bounds, should be from 0 to %d", 16 - 1);
	    check_err = 1;
	}
	break;
    case JUST:
	if (new_val >= 0 && new_val <= MAX_JUST) {
	    retval = new_val;
	} else {
	    sprintf(buf, "JUST out of bounds, should be from 0 to %d", MAX_JUST);
	    check_err = 1;
	}
	break;
    case FONTP:
	if (new_val >= 0 && new_val < MAX_FONT) {
	    retval = new_val;
	} else {
	    sprintf(buf, "FONT out of bounds, should be from 0 to %d", MAX_FONT - 1);
	    check_err = 1;
	}
	break;
    case ARROW:
	if (new_val >= 0 && new_val <= MAX_ARROW) {
	    retval = new_val;
	} else {
	    sprintf(buf, "ARROW out of bounds, should be from 0 to %d", MAX_ARROW);
	    check_err = 1;
	}
	break;
    case PATTERN:
	if (new_val >= 0 && new_val < MAX_PATTERN) {
	    retval = new_val;
	} else {
	    sprintf(buf, "PATTERN out of bounds, should be from 0 to %d", MAX_PATTERN - 1);
	    check_err = 1;
	}
	break;
    case SYMBOL:
	if (new_val >= 0 && new_val < MAXSYM) {
	    retval = new_val;
	} else {
	    sprintf(buf, "SYMBOL out of bounds, should be from 0 to %d", MAXSYM - 1);
	    check_err = 1;
	}
	break;
    case PREC:
	if (new_val >= 0 && new_val < MAX_PREC) {
	    retval = new_val;
	} else {
	    sprintf(buf, "PREC out of bounds, should be from 0 to %d", MAX_PREC - 1);
	    check_err = 1;
	}
	break;
    }
    if (check_err) {
	errwin(buf);
    }
    return retval;
}
