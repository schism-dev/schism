/* $Id: regionutils.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * routines to allocate, manipulate, and return
 * information about regions.
 *
 */

#include <stdio.h>
#include <math.h>
#include "globals.h"

extern int regiontype, regionlinkto;

/*
 * see if (x,y) lies inside the plot
 */
int inbounds(int gno, double x, double y)
{
    return ((x >= g[gno].w.xg1 && x <= g[gno].w.xg2) && (y >= g[gno].w.yg1 && y <= g[gno].w.yg2));
}

int isactive_region(int regno)
{
    return (regno == MAXREGION || regno == MAXREGION + 1 || rg[regno].active == ON);
}

char *region_types(int it, int which)
{
    static char s[128];

    strcpy(s, "UNDEFINED");
    switch (it) {
    case LEFT:
	strcpy(s, "LEFT");
	break;
    case RIGHT:
	strcpy(s, "RIGHT");
	break;
    case ABOVE:
	strcpy(s, "ABOVE");
	break;
    case BELOW:
	strcpy(s, "BELOW");
	break;
    case POLYI:
	if (which) {
	    strcpy(s, "POLYI");
	} else {
	    strcpy(s, "INSIDE POLY");
	}
	break;
    case POLYO:
	if (which) {
	    strcpy(s, "POLYO");
	} else {
	    strcpy(s, "OUTSIDE POLY");
	}
	break;
    }
    return s;
}

void kill_region(int r)
{
    int i;

    if (rg[r].active == ON) {
	if (rg[r].x != NULL) {
	    free(rg[r].x);
	    rg[r].x = NULL;
	}
	if (rg[r].y != NULL) {
	    free(rg[r].y);
	    rg[r].y = NULL;
	}
    }
    rg[r].active = OFF;
    for (i = 0; i < maxgraph; i++) {
	rg[r].linkto[i] = FALSE;
    }
}

void activate_region(int r, int type)
{
    kill_region(r);
    rg[r].active = ON;
    rg[r].type = type;
}

void define_region(int nr, int regionlinkto, int rtype)
{
    kill_region(nr);
    switch (rtype) {
    case 0:
	regiontype = POLYI;
	do_select_region();
	break;
    case 1:
	regiontype = POLYO;
	do_select_region();
	break;
    case 2:
	regiontype = ABOVE;
	set_action(0);
	set_action(DEF_REGION1ST);
	break;
    case 3:
	regiontype = BELOW;
	set_action(0);
	set_action(DEF_REGION1ST);
	break;
    case 4:
	regiontype = LEFT;
	set_action(0);
	set_action(DEF_REGION1ST);
	break;
    case 5:
	regiontype = RIGHT;
	set_action(0);
	set_action(DEF_REGION1ST);
	break;
    }
}

void extract_region(int gno, int fromset, int toset, int regno)
{
    int i, j;
    double *x, *y;
    int startno, stopno;
    if (fromset == -1) {
	startno = 0;
	stopno = maxplot - 1;
    } else {
	startno = stopno = fromset;
    }
    if (regno >= MAXREGION) {
	for (j = startno; j <= stopno; j++) {
	    x = getx(cg, j);
	    y = gety(cg, j);
	    if (isactive(cg, j) && (toset != j || gno != cg)) {
		for (i = 0; i < getsetlength(cg, j); i++) {
		    if (regno == MAXREGION) {
			if (inbounds(cg, x[i], y[i])) {
			    add_point(gno, toset, x[i], y[i], 0.0, 0.0, XY);
			}
		    } else {
			if (!inbounds(cg, x[i], y[i])) {
			    add_point(gno, toset, x[i], y[i], 0.0, 0.0, XY);
			}
		    }
		}
	    }
	}
    } else {
	if (rg[regno].active == OFF) {
	    errwin("Region not active");
	    return;
	}
	if (rg[regno].linkto[cg] == FALSE) {
	    errwin("Region not linked to this graph");
	    return;
	}
	for (j = 0; j < g[cg].maxplot; j++) {
	    x = getx(cg, j);
	    y = gety(cg, j);
	    if (isactive(cg, j) && (toset != j || gno != cg)) {
		for (i = 0; i < getsetlength(cg, j); i++) {
		    if (inregion(regno, x[i], y[i])) {
			add_point(gno, toset, x[i], y[i], 0.0, 0.0, XY);
		    }
		}
	    }
	}
    }
    updatesetminmax(gno, toset);
    update_set_status(gno, toset);
    drawgraph();
}

void delete_byindex(int gno, int setno, int *ind)
{
    int i, j, cnt = 0;
    int ncols = getncols(gno, setno);
    for (i = 0; i < getsetlength(gno, setno); i++) {
	if (ind[i]) {
	    cnt++;
	}
    }
    if (cnt == getsetlength(gno, setno)) {
	killset(gno, setno);
	return;
    }
    cnt = 0;
    for (i = 0; i < getsetlength(gno, setno); i++) {
	if (ind[i] == 0) {
	    for (j = 0; j < ncols; j++) {
		g[gno].p[setno].ex[j][cnt] = g[gno].p[setno].ex[j][i];
	    }
	    cnt++;
	}
    }
    setlength(gno, setno, cnt);
}

void delete_region(int gno, int setno, int regno)
{
    int i, j, k, len, *ind = NULL;
    int gstart, gstop;
    int sstart, sstop;
    double *x, *y;

    if (regno < 0 || regno > MAXREGION + 1) {
	errwin("Invalid region");
	return;
    }
    if (regno < MAXREGION) {
	if (rg[regno].active == OFF) {
	    errwin("Region not active");
	    return;
	}
	if (gno >= 0) {
	    if (rg[regno].linkto[gno] == FALSE) {
		errwin("Region not linked to this graph");
		return;
	    }
	}
    }
    if (gno < 0) {		/* current graph */
	gstart = cg;
	gstop = cg;
    } else if (gno == maxgraph) {	/* all graphs */
	gstart = 0;
	gstop = maxgraph - 1;
    } else {
	gstart = gno;		/* particular graph */
	gstop = gno;
    }

    for (k = gstart; k <= gstop; k++) {
	if (isactive_graph(k)) {
            if (setno < 0) {
                sstart = 0;
                sstop = g[k].maxplot - 1;
            } else {
                sstart = setno;
                sstop = setno;
            }
	    for (j = sstart; j <= sstop; j++) {
		x = getx(k, j);
		y = gety(k, j);
		if (isactive(k, j)) {
		    len = getsetlength(k, j);
		    if (ind != NULL) {
			free(ind);
		    }
		    ind = (int *) malloc(len * sizeof(int));
		    if (ind == NULL) {
			errwin("Error mallocing memory in delete_region, operation cancelled");
			return;
		    }
		    for (i = 0; i < len; i++) {
			ind[i] = 0;
			if (regno == MAXREGION) {	/* inside world */
			    if (inbounds(k, x[i], y[i])) {
				ind[i] = 1;
			    }
			} else if (regno == MAXREGION + 1) {	/* outside world */
			    if (!inbounds(k, x[i], y[i])) {
				ind[i] = 1;
			    }
			} else {
			    if (inregion(regno, x[i], y[i])) {	/* inside region */
				ind[i] = 1;
			    }
			}
		    }
		    delete_byindex(k, j, ind);
		    updatesetminmax(k, j);
		    update_set_status(k, j);
		}
	    }
	}
    }
    drawgraph();
}

void evaluate_region(int regno, char *buf)
{
    double a, b, c, d;
    double *x, *y;
    int errpos;
    int i, j;
    extern double resx, resy;	/* result passed from the expression
				 * interpreter */

    if (regno >= MAXREGION) {
	for (j = 0; j < g[cg].maxplot; j++) {
	    x = getx(cg, j);
	    y = gety(cg, j);
	    if (isactive(cg, j)) {
		for (i = 0; i < getsetlength(cg, j); i++) {
		    if (regno == MAXREGION) {
			if (inbounds(cg, x[i], y[i])) {
			    scanner(buf, &x[i], &y[i], 1, &a, &b, &c, &d, 1, i, j, &errpos);
			    if (errpos) {
				updatesetminmax(cg, j);
				update_set_status(cg, j);
				return;
			    }
			}
		    } else {
			if (!inbounds(cg, x[i], y[i])) {
			    scanner(buf, &x[i], &y[i], 1, &a, &b, &c, &d, 1, i, j, &errpos);
			    if (errpos) {
				updatesetminmax(cg, j);
				update_set_status(cg, j);
				return;
			    }
			}
		    }
		}
		updatesetminmax(cg, j);
		update_set_status(cg, j);
	    }
	}
    } else {
	if (rg[regno].active == OFF) {
	    errwin("Region not active");
	    return;
	}
	if (rg[regno].linkto[cg] == FALSE) {
	    errwin("Region not linked to this graph");
	    return;
	}
	for (j = 0; j < g[cg].maxplot; j++) {
	    x = getx(cg, j);
	    y = gety(cg, j);
	    if (isactive(cg, j)) {
		for (i = 0; i < getsetlength(cg, j); i++) {
		    if (inregion(regno, x[i], y[i])) {
			scanner(buf, &x[i], &y[i], 1, &a, &b, &c, &d, 1, i, j, &errpos);
			if (errpos) {
			    updatesetminmax(cg, j);
			    update_set_status(cg, j);
			    return;
			}
		    }
		}
		updatesetminmax(cg, j);
		update_set_status(cg, j);
	    }
	}
    }
    drawgraph();
}

/*
 * extract sets from graph gfrom to graph gto
 * a set is in a region if any point of the set is in the region
 */
void extractsets_region(int gfrom, int gto, int rno)
{
    int i, j;
    double *x, *y;
    for (j = 0; j < g[gfrom].maxplot; j++) {
	x = getx(gfrom, j);
	y = gety(gfrom, j);
	if (isactive(gfrom, j)) {
	    for (i = 0; i < getsetlength(gfrom, j); i++) {
		if (inregion(rno, x[i], y[i])) {
		}
	    }
	}
    }
}

/*
 * delete sets from graph gno
 * a set is in a region if any point of the set is in the region
 */
void deletesets_region(int gno, int rno)
{
    int i, j;
    double *x, *y;
    for (j = 0; j < g[gno].maxplot; j++) {
	x = getx(gno, j);
	y = gety(gno, j);
	if (isactive(gno, j)) {
	    for (i = 0; i < getsetlength(gno, j); i++) {
		if (inregion(rno, x[i], y[i])) {
		    killset(gno, j);
		    break;	/* set no longer exists, so get out */
		}
	    }
	}
    }
}

void load_poly_region(int r, int n, double *x, double *y)
{
    int i;

    if (n > 2) {
	activate_region(r, regiontype);
	rg[r].n = n;
	rg[r].x = (double *) calloc(n, sizeof(double));
	rg[r].y = (double *) calloc(n, sizeof(double));
	for (i = 0; i < n; i++) {
	    rg[r].x[i] = x[i];
	    rg[r].y[i] = y[i];
	}
    }
}

void draw_region(int r)
{
    int i, c, s, w;
    double vx, vy, wx, wy;

    c = setcolor(rg[r].color);
    s = setlinestyle(rg[r].lines);
    w = setlinewidth(rg[r].linew);
    switch (rg[r].type) {
    case ABOVE:
	my_move2(rg[r].x1, rg[r].y1);
	my_draw2(rg[r].x2, rg[r].y2);
	world2view(rg[r].x1, rg[r].y1, &vx, &vy);
	view2world(vx, vy + 0.05, &wx, &wy);
	draw_arrow(rg[r].x1, rg[r].y1, rg[r].x1, wy, 2, 1.0, 0);
	world2view(rg[r].x2, rg[r].y2, &vx, &vy);
	view2world(vx, vy + 0.05, &wx, &wy);
	draw_arrow(rg[r].x2, rg[r].y2, rg[r].x2, wy, 2, 1.0, 0);
	break;
    case BELOW:
	my_move2(rg[r].x1, rg[r].y1);
	my_draw2(rg[r].x2, rg[r].y2);
	world2view(rg[r].x1, rg[r].y1, &vx, &vy);
	view2world(vx, vy - 0.05, &wx, &wy);
	draw_arrow(rg[r].x1, rg[r].y1, rg[r].x1, wy, 2, 1.0, 0);
	world2view(rg[r].x2, rg[r].y2, &vx, &vy);
	view2world(vx, vy - 0.05, &wx, &wy);
	draw_arrow(rg[r].x2, rg[r].y2, rg[r].x2, wy, 2, 1.0, 0);
	break;
    case LEFT:
	my_move2(rg[r].x1, rg[r].y1);
	my_draw2(rg[r].x2, rg[r].y2);
	world2view(rg[r].x1, rg[r].y1, &vx, &vy);
	view2world(vx - 0.05, vy, &wx, &wy);
	draw_arrow(rg[r].x1, rg[r].y1, wx, rg[r].y1, 2, 1.0, 0);
	world2view(rg[r].x2, rg[r].y2, &vx, &vy);
	view2world(vx - 0.05, vy, &wx, &wy);
	draw_arrow(rg[r].x2, rg[r].y2, wx, rg[r].y2, 2, 1.0, 0);
	break;
    case RIGHT:
	my_move2(rg[r].x1, rg[r].y1);
	my_draw2(rg[r].x2, rg[r].y2);
	world2view(rg[r].x1, rg[r].y1, &vx, &vy);
	view2world(vx + 0.05, vy, &wx, &wy);
	draw_arrow(rg[r].x1, rg[r].y1, wx, rg[r].y1, 2, 1.0, 0);
	world2view(rg[r].x2, rg[r].y2, &vx, &vy);
	view2world(vx + 0.05, vy, &wx, &wy);
	draw_arrow(rg[r].x2, rg[r].y2, wx, rg[r].y2, 2, 1.0, 0);
	break;
    case POLYI:
    case POLYO:
	if (rg[r].x != NULL && rg[r].n > 2) {
	    my_move2(rg[r].x[0], rg[r].y[0]);
	    for (i = 1; i < rg[r].n; i++) {
		my_draw2(rg[r].x[i], rg[r].y[i]);
	    }
	    my_draw2(rg[r].x[0], rg[r].y[0]);
	}
	break;
    }
    setcolor(c);
    setlinestyle(s);
    setlinewidth(w);
}

/*
 * routines to determine if a point lies in a polygon
*/
int intersect_to_left(double x, double y, double x1, double y1, double x2, double y2)
{
    double xtmp, m, b;

    /* ignore horizontal lines */
    if (y1 == y2) {
	return 0;
    }
    /* not contained vertically */
    if (((y < y1) && (y < y2)) || ((y > y1) && (y > y2))) {
	return 0;
    }
    /* none of the above, compute the intersection */
    if ((xtmp = x2 - x1) != 0.0) {
	m = (y2 - y1) / xtmp;
	b = y1 - m * x1;
	xtmp = (y - b) / m;
    } else {
	xtmp = x1;
    }
    if (xtmp <= x) {
	/* check for intersections at a vertex */
	/* if this is the max ordinate then accept */
	if (y == y1) {
	    if (y1 > y2) {
		return 1;
	    } else {
		return 0;
	    }
	}
	/* check for intersections at a vertex */
	if (y == y2) {
	    if (y2 > y1) {
		return 1;
	    } else {
		return 0;
	    }
	}
	/* no vertices intersected */
	return 1;
    }
    return 0;
}

/*
 * determine if (x,y) is in the polygon xlist[], ylist[]
 */
int inbound(double x, double y, double *xlist, double *ylist, int n)
{
    int i, l = 0, ll = 0;

    for (i = 0; i < n; i++) {
	l += intersect_to_left(x, y, xlist[i], ylist[i], xlist[(i + 1) % n], ylist[(i + 1) % n]);
    }
    return l % 2;
}

/*
 * routines to determine if a point lies to the left of an infinite line
*/
int isleft(double x, double y, double x1, double y1, double x2, double y2)
{
    double xtmp, m, b;

    /* horizontal lines */
    if (y1 == y2) {
	return 0;
    }
    /* none of the above, compute the intersection */
    if ((xtmp = x2 - x1) != 0.0) {
	m = (y2 - y1) / xtmp;
	b = y1 - m * x1;
	xtmp = (y - b) / m;
    } else {
	xtmp = x1;
    }
    if (xtmp >= x) {
	return 1;
    }
    return 0;
}


/*
 * routines to determine if a point lies to the left of an infinite line
*/
int isright(double x, double y, double x1, double y1, double x2, double y2)
{
    double xtmp, m, b;

    /* horizontal lines */
    if (y1 == y2) {
	return 0;
    }
    if ((xtmp = x2 - x1) != 0.0) {
	m = (y2 - y1) / xtmp;
	b = y1 - m * x1;
	xtmp = (y - b) / m;
    } else {
	xtmp = x1;
    }
    if (xtmp <= x) {
	return 1;
    }
    return 0;
}

/*
 * routines to determine if a point lies above an infinite line
*/
int isabove(double x, double y, double x1, double y1, double x2, double y2)
{
    double xtmp, ytmp, m, b;

    /* vertical lines */
    if (x1 == x2) {
	return 0;
    }
    if ((ytmp = y2 - y1) != 0.0) {
	m = ytmp / (x2 - x1);
	b = y1 - m * x1;
	ytmp = m * x + b;
    } else {
	ytmp = y1;
    }
    if (ytmp <= y) {
	return 1;
    }
    return 0;
}

/*
 * routines to determine if a point lies below an infinite line
*/
int isbelow(double x, double y, double x1, double y1, double x2, double y2)
{
    double xtmp, ytmp, m, b;

    /* vertical lines */
    if (x1 == x2) {
	return 0;
    }
    if ((ytmp = y2 - y1) != 0.0) {
	m = ytmp / (x2 - x1);
	b = y1 - m * x1;
	ytmp = m * x + b;
    } else {
	ytmp = y1;
    }
    if (ytmp >= y) {
	return 1;
    }
    return 0;
}

int inregion(int regno, double x, double y)
{
    int i;

    if (regno == MAXREGION) {
	return (inbounds(cg, x, y));
    }
    if (regno == MAXREGION + 1) {
	return (!inbounds(cg, x, y));
    }
    if (rg[regno].active == ON) {
	switch (rg[regno].type) {
	case POLYI:
	    if (inbound(x, y, rg[regno].x, rg[regno].y, rg[regno].n)) {
		return 1;
	    }
	    break;
	case POLYO:
	    if (!inbound(x, y, rg[regno].x, rg[regno].y, rg[regno].n)) {
		return 1;
	    }
	    break;
	case RIGHT:
	    if (isright(x, y, rg[regno].x1, rg[regno].y1, rg[regno].x2, rg[regno].y2)) {
		return 1;
	    }
	    break;
	case LEFT:
	    if (isleft(x, y, rg[regno].x1, rg[regno].y1, rg[regno].x2, rg[regno].y2)) {
		return 1;
	    }
	    break;
	case ABOVE:
	    if (isabove(x, y, rg[regno].x1, rg[regno].y1, rg[regno].x2, rg[regno].y2)) {
		return 1;
	    }
	    break;
	case BELOW:
	    if (isbelow(x, y, rg[regno].x1, rg[regno].y1, rg[regno].x2, rg[regno].y2)) {
		return 1;
	    }
	    break;
	}
    }
    return 0;
}
