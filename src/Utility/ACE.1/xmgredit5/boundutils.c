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
 * utilities for manipulating boundaries
 *
 */

#ifndef lint
static char RCSid[] = "$Id: boundutils.c,v 1.3 2006/07/31 23:42:52 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "defines.h"
#include "globals.h"

#define MAXBNODES 600000

void find_boundary_point(int gridno, double wx, double wy, int *bno, int *ind);
void check_gridboundary_intersect(int gridno);
void check_gridboundary_edges(int gridno);

void move_boundary_node(int gridno, int bno, int ind, double wx, double wy)
{
    boundary[bno].boundx[ind] = wx;
    boundary[bno].boundy[ind] = wy;
}

void delete_boundary_node(int gridno, double wx, double wy)
{
    int ind, bno, i;

    find_boundary_point(gridno, wx, wy, &bno, &ind);
    if (bno < 0)
	return;
    if (ind < 0)
	return;
    for (i = ind + 1; i < boundary[bno].nbpts; i++) {
	boundary[bno].boundx[i] = boundary[bno].boundx[i - 1];
	boundary[bno].boundy[i] = boundary[bno].boundy[i - 1];
    }
    Reallocate_boundary(bno, boundary[bno].nbpts - 1);
}

void add_boundary_node(int gridno, int bno, int b1, int b2, double wx, double wy)
{
    int i;
    Reallocate_boundary(bno, boundary[bno].nbpts + 1);
    for (i = boundary[bno].nbpts - 1; i > b2; i--) {
	boundary[bno].boundx[i] = boundary[bno].boundx[i - 1];
	boundary[bno].boundy[i] = boundary[bno].boundy[i - 1];
    }
    boundary[bno].boundx[b2] = wx;
    boundary[bno].boundy[b2] = wy;
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
	if ((y < ylist[i] && y < ylist[(i + 1) % n]) || (y > ylist[i] && y > ylist[(i + 1) % n])) {
	    continue;
	}
	if ((x < xlist[i] && x < xlist[(i + 1) % n])) {
	    continue;
	}
	l += intersect_to_left(x, y, xlist[i], ylist[i], xlist[(i + 1) % n], ylist[(i + 1) % n]);
    }
    return l % 2;
}

void clip_to_boundary(int gridno, double *bx, double *by, int *npts)
{
    int i, ngood = 0;

    for (i = 0; i < *npts; i++) {
	if (goodpoint(gridno, bx[i], by[i])) {
	    bx[ngood] = bx[i];
	    by[ngood] = by[i];
	    ngood++;
	}
    }
    *npts = ngood;
}

int goodpoint(int gridno, double x, double y)
{
    int i, ib, itmp;

    for (i = 0; i < grid[gridno].nbounds; i++) {
	ib = grid[gridno].boundaries[i];
	if (boundary[ib].bactive) {
	    itmp = inbound(x, y, boundary[ib].boundx, boundary[ib].boundy, boundary[ib].nbpts);
	    if (boundary[ib].boundtype == 0) {
		if (itmp == 0) {
		    return 0;
		}
	    } else if (itmp == 1) {
		return 0;
	    }
	}
    }
    return 1;
}

void find_boundary_point(int gridno, double wx, double wy, int *bno, int *ind)
{
    int i, j, ib;
    double tmp, radius = 1e307;

    if (grid[gridno].nbounds <= 0) {
	errwin("No boundaries");
	*bno = -1;
	*ind = -1;
	return;
    }
    for (i = 0; i < grid[gridno].nbounds; i++) {
	ib = grid[gridno].boundaries[i];
	for (j = 0; j < boundary[ib].nbpts; j++) {
	    tmp = hypot((wx - boundary[ib].boundx[j]), (wy - boundary[ib].boundy[j]));
	    if (tmp < radius) {
		radius = tmp;
		*ind = j;
		*bno = ib;
	    }
	}
    }
}

void find_external_boundary_point(int gridno, double wx, double wy, int *bno, int *ind)
{
    int j, ib;
    double tmp, radius = 1e307;

    *bno = -1;
    *ind = -1;
    if (grid[gridno].nbounds <= 0) {
	errwin("No external boundary");
	return;
    }
    ib = grid[gridno].boundaries[0];
    *bno = ib;
    for (j = 0; j < boundary[ib].nbpts; j++) {
	tmp = hypot((wx - boundary[ib].boundx[j]), (wy - boundary[ib].boundy[j]));
	if (tmp < radius) {
	    radius = tmp;
	    *ind = j;
	}
    }
}

void find_internal_boundary_point(int gridno, double wx, double wy, int *bno, int *ind)
{
    int i, j, ib;
    double tmp, radius = 1e307;

    if (grid[gridno].nbounds <= 1) {
	errwin("No internal boundaries");
	*bno = -1;
	*ind = -1;
	return;
    }
    for (i = 1; i < grid[gridno].nbounds; i++) {
	ib = grid[gridno].boundaries[i];
	for (j = 0; j < boundary[ib].nbpts; j++) {
	    tmp = hypot((wx - boundary[ib].boundx[j]), (wy - boundary[ib].boundy[j]));
	    if (tmp < radius) {
		radius = tmp;
		*ind = j;
		*bno = ib;
	    }
	}
    }
}

void associate_grid(int gridno, int bno)
{
    grid[gridno].boundaries[grid[gridno].nbounds] = bno;
    grid[gridno].nbounds++;
}

void disassociate_grid(int gridno, int bno)
{
    int i, k = 0, found = 0;

    if (grid[gridno].nbounds) {
	for (i = 0; i < grid[gridno].nbounds; i++) {
	    if (bno == grid[gridno].boundaries[i]) {
		if (found == 1) {
		    errwin("Call Paul, assumption incorrect in dis_bound()\n");
		}
		Free_boundary(bno);
		found = 1;
	    } else {
		grid[gridno].boundaries[k] = grid[gridno].boundaries[i];
		k++;
	    }
	}
    }
    if (found) {
	grid[gridno].nbounds--;
    } else {
	errwin("No boundary found, internal error in dis_bound()\n");
    }
}

static int nb;

int findadjelem(int gridno, int el, int n1, int n2)
{
    int i, j, k;

    for (i = 0; i < grid[gridno].nmel; i++) {
	if (!(grid[gridno].ellist[i] == 3)) {
	    if (i != el) {
		for (j = 0; j < 3; j++) {
		    if (n1 == grid[gridno].icon[i].nl[j]) {
			for (k = 0; k < 3; k++) {
			    if (n2 == grid[gridno].icon[i].nl[k]) {
				return i;
			    }
			}
		    }
		}
	    }
	}
    }
    return -1;
}

void getboundary(int gridno, int (*clist)[2], int *nb)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[40], v, nmnp;
    int nn, mnode;
    double xmin;
    struct edgelist {
	int *list;
    } *e;

    nmnp = grid[gridno].nmnp;
    e = (struct edgelist *) malloc(nmnp * sizeof(struct edgelist));
    if (e == NULL) {
	errwin("Can't malloc edgelist");
	return;
    }
    for (i = 0; i < nmnp; i++) {
	e[i].list = (int *) calloc(2, sizeof(int));
	e[i].list[0] = -1;
	e[i].list[1] = 0;
    }
/* find the starting node */
    xmin = grid[gridno].xord[0];
    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (grid[gridno].xord[i] <= xmin) {
	    xmin = grid[gridno].xord[i];
	    mnode = i;
	}
    }

    for (i = 0; i < grid[gridno].nmel; i++) {
/*
	if (grid[gridno].icon[i].type != grid[gridno].icon[i].ngeom) {
*/
	switch (grid[gridno].icon[i].type) {
	case 0:
	    for (j = 0; j < 3; j++) {
		el[j] = grid[gridno].icon[i].nl[j];
	    }
	    break;
	case 1:
	    for (j = 0; j < 3; j++) {
		el[j] = grid[gridno].icon[i].nl[j];
	    }
	    break;
	case 2:
	    for (j = 0; j < 2; j++) {
		el[j] = grid[gridno].icon[i].nl[j];
	    }
	    break;
	case 3:
	    for (j = 0; j < 3; j++) {
		el[j] = grid[gridno].icon[i].nl[j];
	    }
	    break;
	case 4:
	    for (j = 0; j < 4; j++) {
		el[j] = grid[gridno].icon[i].nl[j];
	    }
	    break;
	case 5:
	    for (j = 0; j < 5; j++) {
		el[j] = grid[gridno].icon[i].nl[j];
	    }
	    break;
	case 6:
	    for (j = 0; j < 3; j++) {
		el[2 * j] = grid[gridno].icon[i].nl[j];
	    }
	    for (j = 0; j < 3; j++) {
		el[2 * j + 1] = grid[gridno].icon[i].nl[j + 3];
	    }
	    break;
	case 8:
	    for (j = 0; j < 4; j++) {
		el[2 * j] = grid[gridno].icon[i].nl[j];
	    }
	    for (j = 0; j < 4; j++) {
		el[2 * j + 1] = grid[gridno].icon[i].nl[j + 4];
	    }
	    break;
	default:
	    for (j = 0; j < 10; j++) {
		el[j] = -1;
	    }
	    break;
	}
	/*
	 * if (grid[gridno].icon[i].type != 1 && grid[gridno].icon[i].type !=
	 * 5) {
	 */
	switch (grid[gridno].icon[i].type) {
	case 0:
	    nn = 0;
	    break;
	case 1:
	    nn = 1;
	    break;
	case 2:
	    nn = 2;
	    break;
	case 5:
	    nn = 2;
	    break;
	default:
	    nn = grid[gridno].icon[i].type;
	    break;
	}
	if (nn > 2) {
	    for (j = 0; j < nn; j++) {
/* search for edge */
/*
		if (nn == 2) {
		    if (j == 0) {
			n1 = el[1];
			n2 = el[0];
		    } else if (j == 1) {
			n1 = el[2];
			n2 = el[1];
		    }
		} else {
		    n1 = el[j];
		    n2 = el[(j + 1) % nn];
		}
*/
		n1 = el[j];
		n2 = el[(j + 1) % nn];
		cnt1 = 0;
		cnt2 = 0;
		f1 = 0;
		f2 = 0;

		while (e[n1].list[2 * cnt1] != -1 && e[n1].list[2 * cnt1] != n2) {
		    cnt1++;
		}
		while (e[n2].list[2 * cnt2] != -1 && e[n2].list[2 * cnt2] != n1) {
		    cnt2++;
		}

		if (e[n1].list[2 * cnt1] == -1) {
		    f1 = 0;
		} else {
		    f1 = 1;
		}
		if (e[n2].list[2 * cnt2] == -1) {
		    f2 = 0;
		} else {
		    f2 = 1;
		}

/* check cnt1 and cnt2, should only appear once */
		if (f1 && f2) {
		    printf("Error: duplicate edge, send mail to ace_bugs, f1 = %d f2 = %d element number = %d element type = %d\n", f1, f2, i + 1, grid[gridno].icon[i].type);
		    exit(1);
		}
/* if cnt1 != 0 then increment element counter */
		if (f1 && !f2) {
		    e[n1].list[2 * cnt1 + 1]++;
		}
/* if cnt2 != 0 then increment element counter */
		else if (f2 && !f1) {
		    e[n2].list[2 * cnt2 + 1]++;
		}
/* if cnt1 == 0 && cnt2 == 0 then add edge n1-n2 and increment counter */
		else if (!f1 && !f2) {
		    e[n1].list = (int *) realloc(e[n1].list, 2 * (cnt1 + 2) * sizeof(int));
		    e[n1].list[2 * cnt1] = n2;
		    e[n1].list[2 * cnt1 + 1]++;
		    cnt1++;
		    e[n1].list[2 * cnt1] = -1;
		    e[n1].list[2 * cnt1 + 1] = 0;
		} else {
		    printf("Error: duplicate edge (phase 2), send mail to ace_bugs\n");
		    exit(1);
		}
	    }
	}
    }

    *nb = 0;
    for (i = 0; i < nmnp; i++) {
	cnt1 = 0;
	while (e[(i + mnode) % nmnp].list[2 * cnt1] != -1) {
	    if (e[(i + mnode) % nmnp].list[2 * cnt1 + 1] == 1) {
		clist[*nb][0] = (i + mnode) % nmnp;
		clist[*nb][1] = e[(i + mnode) % nmnp].list[2 * cnt1];
		*nb += 1;
	    }
	    cnt1++;
	}
    }
    for (i = 0; i < nmnp; i++) {
	if (e[i].list != NULL) {
	    free(e[i].list);
	}
    }
    free(e);
}

/*
 * Update the boundary menu on the main panel
 */
void update_boundary_menu(int gridno)
{
    if (grid[gridno].nbounds > 1) {
	ebound_defined = 1;
	ibound_defined = 1;
	update_fuzz_items();
    }
    if (grid[gridno].nbounds == 1) {
	ibound_defined = 0;
	ebound_defined = 1;
	update_fuzz_items();
    }
}

void load_boundary(int gridno, int (*clist)[2])
{
    int start, bno, done, n = 0, *outb;
    int i, j, n1, n2, n3, found, cnt, next, tcnt = 0, skipcnt = 0, *skipn;
    double tmpx, tmpy;
    Boundary b;
    outb = (int *) malloc(MAXBNODES * sizeof(int));
    if (outb == NULL) {
	return;
    }
    skipn = (int *) malloc(MAXBNODES * sizeof(int));
    if (skipn == NULL) {
	free(outb);
	return;
    }
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = 0;
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].type == 5) {
	    skipn[skipcnt++] = grid[gridno].icon[i].nl[1];
	}
    }
    skipcnt = 0;
    if (skipcnt) {
	for (i = 0; i < nb; i++) {
	    for (j = 0; j < skipcnt; j++) {
		if (clist[i][0] == skipn[j]) {
		    clist[i][0] = -1;
		    /* clist[i][1] = -1; */
		}
		if (clist[i][1] == skipn[j]) {
		    clist[i][1] = -1;
		    /* clist[i][0] = -1; */
		}
	    }
	}
    }
    if (skipcnt) {
	for (i = 0; i < nb; i++) {
	    if (clist[i][1] == -1 && clist[i][0] != -1) {
		next = clist[i][0];
		outb[0] = clist[i][0];
		clist[i][0] = -1;
		break;
	    }
/*
else if (clist[i][1] == -1 && clist[i][0] != -1) {
		next = clist[i][0];
		outb[0] = clist[i][0];
		clist[i][0] = -1;
		break;
	    }
*/
	}
    } else {
	outb[0] = clist[0][0];
	next = clist[0][1];
	clist[0][1] = -1;
    }
    cnt = 1;
    done = 0;
    while (!done) {
	for (i = 0; i < nb; i++) {
	    for (j = 0; j < nb; j++) {
		if (next == clist[j][0] && next != -1) {
		    outb[cnt++] = next;
		    next = clist[j][1];
		    clist[j][0] = -1;
		    clist[j][1] = -1;
		} else if (next == clist[j][1] && next != -1) {
		    outb[cnt++] = next;
		    next = clist[j][0];
		    clist[j][0] = -1;
		    clist[j][1] = -1;
		}
		if (next == -1)
		    goto getout;
	    }
	}
getout:;
	if (cnt > 0) {
	    bno = nextboundary();
	    if (bno == -1) {
		errwin("Maximum number of boundaries (1000) exceeded");
		goto giveup;
	    }
	    grid[gridno].boundaries[n] = bno;
	    n++;
	    if (skipcnt == 0) {
		if (n == 1) {
		    Allocate_boundary(bno, cnt - 1, 0);
		} else if (n > 1) {
		    Allocate_boundary(bno, cnt - 1, 1);
		}
		for (i = 0; i < cnt - 1; i++) {
		    boundary[bno].boundx[i] = grid[gridno].xord[outb[i]];
		    boundary[bno].boundy[i] = grid[gridno].yord[outb[i]];
		    grid[gridno].nlist[outb[i]] = 1;
		    boundary[bno].nodes[i] = outb[i];
		}
	    } else {
		Allocate_boundary(bno, cnt, 3);
		for (i = 0; i < cnt; i++) {
		    boundary[bno].boundx[i] = grid[gridno].xord[outb[i]];
		    boundary[bno].boundy[i] = grid[gridno].yord[outb[i]];
		    grid[gridno].nlist[outb[i]] = 1;
		    boundary[bno].nodes[i] = outb[i];
		}
	    }
	    cnt = 0;
	    done = 1;
	    for (i = 0; i < nb; i++) {
		if (!skipcnt) {
		    if (clist[i][0] != -1 && clist[i][1] != -1) {
			outb[cnt++] = clist[i][0];
			next = clist[i][1];
			clist[i][1] = -1;
			done = 0;
			goto getout2;
		    }
		} else {
		    if (clist[i][0] == -1 && clist[i][1] != -1) {
			next = clist[i][1];
			outb[cnt++] = clist[i][1];
			clist[i][1] = -1;
			done = 0;
			goto getout2;
		    } else if (clist[i][1] == -1 && clist[i][0] != -1) {
			next = clist[i][0];
			outb[cnt++] = clist[i][0];
			clist[i][0] = -1;
			done = 0;
			goto getout2;
		    }
		}
	    }
	} else {
	    done = 1;
	}
getout2:;
    }
giveup:;
    cnt = 0;
    for (i = 0; i < cnt; i++) {
	found = outb[i];
	for (j = 0; j < cnt; j++) {
	    if (i != j) {
		if (outb[i] == outb[j]) {
		    found = -1;
		}
	    }
	}
	if (found >= 0) {
	    bno = nextboundary();
	    if (bno == -1) {
		errwin("Maximum number of boundaries (1000) exceeded");
		goto giveup2;
	    }
	    grid[gridno].boundaries[n] = bno;
	    n++;
	    Allocate_boundary(bno, 1, 3);
	    boundary[bno].boundx[0] = grid[gridno].xord[found];
	    boundary[bno].boundy[0] = grid[gridno].yord[found];
	    boundary[bno].nodes[0] = found;
	}
    }
giveup2:;
    grid[gridno].nbounds = n;
    free(outb);
    free(skipn);
    if (n > 1) {
	ebound_defined = 1;
	ibound_defined = 1;
	update_fuzz_items();
    }
    if (n == 1) {
	ibound_defined = 0;
	ebound_defined = 1;
	update_fuzz_items();
    }
}

static int clist[MAXBNODES][2];
void compute_boundary(int gridno)
{
    int i;

    nb = 0;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	Free_boundary(grid[gridno].boundaries[i]);
    }
    grid[gridno].nbounds = 0;
    getboundary(gridno, clist, &nb);
    load_boundary(gridno, clist);
}

extern int dispbound_flag;

void define_extbound(void)
{
    set_action(0);
    write_mode_str("Define external boundary - use left button to mark, middle button to register");
    set_action(DEF_EXT_BOUND);
}

void move_extbound(void)
{
    set_action(0);
    write_mode_str("Moving external boundary points - use left button to mark and again to move");
    set_action(MOVE_EXT_BOUND1ST);
}

void insert_boundpt(void)
{
    set_action(0);
    write_mode_str("Insert boundary point - mark first node");
    set_action(ADD_BOUND1);
}

void define_intbound(void)
{
    set_action(0);
    write_mode_str("Define an internal boundary - use left button to mark, middle button to register");
    set_action(DEF_INT_BOUND);
}

void move_intbound(void)
{
    set_action(0);
    write_mode_str("Moving internal boundary points - use left button to mark and again to move");
    set_action(MOVE_INT_BOUND1ST);
}

void del_intbound(void)
{
    set_action(0);
    write_mode_str("Delete an internal boundary - use left button to mark");
    set_action(DEL_INT_BOUND);
}

void compute_bound(void)
{
    set_action(0);
    write_mode_str("Computing boundary from editable grid, please wait...");
    compute_boundary(curgrid);
    draw_boundary(curgrid);
    write_mode_str(NULL);
}

void clear_extbound(void)
{
    write_mode_str("Clearing external boundary definition");
    if (grid[curgrid].nbounds > 0) {
	Free_boundary(grid[curgrid].boundaries[0]);
	grid[curgrid].boundaries[0] = -1;
	grid[curgrid].nbounds--;
	ebound_defined = 0;
	update_fuzz_items();
    }
    do_drawgrid();
}

void clear_intbound(void)
{
    int i;
    write_mode_str("Clearing internal boundary definition(s)");
    for (i = 1; i < grid[curgrid].nbounds; i++) {
	Free_boundary(grid[curgrid].boundaries[i]);
	grid[curgrid].boundaries[i] = -1;
	grid[curgrid].nbounds--;
    }
    if (grid[curgrid].nbounds < 0) {
	grid[curgrid].nbounds = 0;
    }
    ibound_defined = 0;
    update_fuzz_items();
    do_drawgrid();
}

void line_bound(void)
{
    write_mode_str("Draw boundary as line");
    dispbound_flag = 0;
    do_drawgrid();
}

void points_bound(void)
{
    write_mode_str("Draw boundary as points");
    dispbound_flag = 1;
    do_drawgrid();
}

void numbers_bound(void)
{
    write_mode_str("Display boundary numbers");
    dispbound_flag = 2;
    do_drawgrid();
}

void markboundaries(int gridno)
{
    int i, j, ib, itmp;

    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = 0;
    }
    for (i = 0; i < grid[gridno].nbounds; i++) {
	ib = grid[gridno].boundaries[i];
	if (boundary[ib].bactive) {
	    for (j = 0; j < boundary[ib].nbpts; j++) {
		itmp = boundary[ib].nodes[j];
		grid[gridno].nlist[itmp] = 1;
	    }
	}
    }
}

int locate_bound_pts(int n, double *bs, double dist, int *i1, int *i2)
{
    int j, ib;
    j = 0;
    while (j < n && bs[j] < dist)
	j++;
    *i1 = j - 1;
    *i2 = j % n;
    return 1;
}

int compute_bound_dist(int gridno, int bno, double x, double y, double dist, int *i1, int *i2)
{
    int i, j, ib;
    double x1, y1, x2, y2;
    if (grid[gridno].nbounds == 0) {
	errwin("No boundary\n");
	return 0;
    }
    if (boundary[bno].bactive) {
	j = 0;
	while (j < boundary[bno].nbpts < dist) {
	    j++;
	}
	*i1 = j - 1;
	*i2 = j % boundary[bno].nbpts;
    }
    return 1;
}

/*
 * return TRUE if area is greater than 0.0
 * FALSE otherwise
 */
int check_boundary(int ib)
{
    double a;
    double comp_area();
    a = comp_area(boundary[ib].nbpts, boundary[ib].boundx, boundary[ib].boundy);
    return a > 0.0;
}

/*
 * reverse the order of the boundary
 */
void reverse_boundary(int ib)
{
    int i, j;
    double *x = boundary[ib].boundx;
    double *y = boundary[ib].boundy;
    int n = boundary[ib].nbpts;
    for (i = 0; i < n / 2; i++) {
	j = (n - 1) - i;
	fswap(&x[i], &x[j]);
	fswap(&y[i], &y[j]);
    }
}

/*
 * check the boundaries of grid gridno for proper
 * orientation
 */
void check_boundaries(int gridno)
{
    int i, ib;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	ib = grid[gridno].boundaries[i];
	if (boundary[ib].bactive) {
	    if (i == 0) {	/* check external boundary */
		if (!check_boundary(ib)) {
		    /* not properly formed */
		    reverse_boundary(ib);
		    errwin("External boundary not properly oriented\nShould be counterclockwise - reversed");
		}
	    } else {
		if (check_boundary(ib)) {
		    /* not properly formed */
		    reverse_boundary(ib);
		    errwin("Internal boundary not properly oriented\nShould be clockwise - reversed");
		}
	    }
	}
    }
}

/*
 * Compute the cumulative distance along a boundary
 */
double *compute_dist_array(int ib)
{
    double dist, *bs;
    double x1, y1, x2, y2;
    int j;
    dist = 0.0;
    bs = (double *) calloc(boundary[ib].nbpts + 2, sizeof(double));
    for (j = 0; j < boundary[ib].nbpts; j++) {
	x1 = boundary[ib].boundx[j];
	y1 = boundary[ib].boundy[j];
	x2 = boundary[ib].boundx[(j + 1) % boundary[ib].nbpts];
	y2 = boundary[ib].boundy[(j + 1) % boundary[ib].nbpts];
	dist = dist + hypot(x1 - x2, y1 - y2);
	bs[j] = dist;
    }
    return bs;
}

double *compute_dist_array2(int ib, int start, double x, double y, double d)
{
    double *bs;
    double dist;
    double x1, y1, x2, y2;
    int j;
    bs = (double *) calloc(boundary[ib].nbpts + 2, sizeof(double));
    for (j = start; j < boundary[ib].nbpts; j++) {
	x1 = boundary[ib].boundx[j];
	y1 = boundary[ib].boundy[j];
	bs[j - start] = hypot(x - x1, y - y1);
    }
    return bs;
}

/*
 * Walk the coastal boundary and place points in the edit grid boundary
 * such that the selected criteria is satisfied.
 *
 * p1 = maxcour
 * p2 = dt
 * p3 = tol
 *
 */
int do_newbound(int gridno, int crit, double p1, double p2, double p3, int dgrid)
{
    int accept, iter, i1, i2, j1, j2, i, j, ib, itmp1, itmp2, cur, cnt, bno;
    double x1, y1, x2, y2, distx, disty, dist, d1, d2, dx, sumdx;
    double dt, cu, A0, s0, tmp, tmp1, tmp2, c1, c2, x0, y0, a, b;
    double t1, t2, xn, yn;
    extern double cboundmin, celevmin;
    char buf[256];
    int errcnt = 0;
    time_t tp, st, et;		/* if timing needed */
/*
    FILE *fp = fopen("debug.d", "w");
*/

    double *bxf, *byf, *bs, *bs2;
    int *bnd;
    double xf, yf, *xs, dx1, dy1, px[15000], py[15000], theta, left;
    int *nd;
    int ind, starti = 0, ip[15000], tcnt = 0;
    double get_depth_element(), comp_area();
    st = time(0);
    switch (crit) {
    case 0:
	break;
    case 1:			/* courant number */
	cu = p1;
	dt = p2;
	break;
    case 2:			/* dimensionless wavel */
	cu = p1;
	dt = p2;
	break;
    }
    if (grid[dgrid].nbounds == 0) {
	errwin("No coastal outline from which to generate boundary\n");
	return 0;
    }
    for (i = 0; i < grid[gridno].nbounds; i++) {
	Free_boundary(grid[gridno].boundaries[i]);
    }
    grid[gridno].nbounds = 0;
    write_mode_str("Processing coastal outlines, please wait...");
    create_timer_popup();
    for (i = 0; i < grid[dgrid].nbounds; i++) {
	ib = grid[dgrid].boundaries[i];
	if (boundary[ib].bactive) {

	    dist = 0.0;
	    sumdx = 0.0;
	    starti = 0;
	    cnt = 1;
	    cur = 1;
/* compute the distance along the boundary (follows the boundary) NOT USED*/
/*	    bs = compute_dist_array(ib); */
/* init for first point */
	    xn = px[0] = boundary[ib].boundx[0];
	    yn = py[0] = boundary[ib].boundy[0];
	    ip[0] = tcnt++;
	    i1 = 0;
	    i2 = 1;
	    errcnt = 0;

	    et = time(0) - st;
	    sprintf(buf, "Processing coastal outlines\nStart\nElapsed time %ld seconds", et);
	    set_timer_item(buf, 0.0, (double) (boundary[ib].nbpts), 0.0);
	    while (cur < boundary[ib].nbpts) {

		et = time(0) - st;
		sprintf(buf, "Processing coastal outlines\nAdjusted %d elevations\nElapsed time %ld seconds", errcnt, et);
		set_timer_item(buf, 0.0, (double) (boundary[ib].nbpts), (double) (cur));
		iter = 0;
		while (iter < 4) {
		    if (debuglevel == 5) {
			printf("current xn, yn = %lf %lf\n", xn, yn);
		    }
		    find_element(dgrid, xn, yn, &ind);
		    if (ind == -1) {
			d1 = cboundmin;
			printf("Element not found at start of loop\n");
		    } else {
			d1 = get_depth_element(dgrid, ind, xn, yn);
			if (d1 < cboundmin) {
			    d1 = cboundmin;
			    errcnt++;
			}
		    }
/* init for first iteration */
		    if (iter == 0) {
			A0 = (M_PI * 9.8 * d1 * (dt * dt)) / (4 * cu * cu);
			s0 = sqrt(2 * A0 / (sin(M_PI / 3.0)));
			if (s0 < 2.0)
			    s0 = 2.0;
			if (debuglevel == 5) {
			    printf("initial A0, d1, and s0 = %lf %lf %lf\n", A0, d1, s0);
			}
		    }
/* compute the distances from xn, yn to succeeding boundary points */
		    bs2 = compute_dist_array2(ib, cur, xn, yn, s0);
/* find the 2 boundary points containing s0 */
		    j = cur;
		    while (j < boundary[ib].nbpts && s0 > bs2[j - cur]) {
			j++;
		    }
		    free(bs2);
		    j1 = j - 1;
		    j2 = j % boundary[ib].nbpts;
		    if (j1 < 0 || j2 < cur) {
			cur = boundary[ib].nbpts;
			if (debuglevel == 5) {
			    printf("bail out at j1 < 0 || j2 < cur, %d %d %d\n", j1, j2, cur);
			}
			goto out1;
		    }
		    if (debuglevel == 6) {
			printf("s0 in %d %d = %lf %lf %lf\n", j1, j2, s0, A0, d1);
		    }
		    if (debuglevel == 5) {
			printf("s0 in %d %d = %lf %lf %lf\n", j1, j2, s0, A0, d1);
		    }
/* compute the intersection of s0 with the boundary */
		    x1 = boundary[ib].boundx[j1];
		    y1 = boundary[ib].boundy[j1];
		    x2 = boundary[ib].boundx[j2];
		    y2 = boundary[ib].boundy[j2];
		    if (debuglevel == 5) {
			printf("boundary points = %lf %lf %lf %lf\n", x1, y1, x2, y2);
		    }
		    c1 = x1 - xn;
		    c2 = y1 - yn;
		    a = x2 - x1;
		    b = y2 - y1;
		    tmp = pow((c1 * a + c2 * b), 2.0) - (a * a + b * b) * (c1 * c1 + c2 * c2 - s0 * s0);
		    if (debuglevel == 5) {
			printf("discriminant = %lf\n", tmp);
		    }
/* using the parametric form of the line, compute the intersection */
		    if (tmp >= 0.0) {
			t1 = (-(c1 * a + c2 * b) + sqrt(tmp)) / (a * a + b * b);
			t2 = (-(c1 * a + c2 * b) - sqrt(tmp)) / (a * a + b * b);
			if (t1 >= 0.0 && t1 <= 1.0) {
			    xf = x1 + a * t1;
			    yf = y1 + b * t1;
			    if (debuglevel == 5) {
				printf("%lf %lf %lf %lf\n", xf, yf, t1, t2);
			    }
			} else if (t2 >= 0.0 && t2 <= 1.0) {
			    xf = x1 + a * t1;
			    yf = y1 + b * t1;
			    if (debuglevel == 5) {
				printf("%lf %lf %lf %lf\n", xf, yf, t1, t2);
			    }
			} else {
			    printf("No intersection\n");
			}
		    } else {
			printf("Discriminant < 0\n");
		    }
/* compute the location of the point on the interior of the boundary
   that forms a equilateral triangle with xn, yn and the intersection
   of s0 with the boundary
   px[], py[] are the locations of the points along the bounday so far.
   xf, yf is the location of the current candidate for inclusion into
   px[], py[]. The distance is computed from the previous point (px[cnt - 1],
   py[cnt - 1]).
*/
		    {
			double dx = xf - px[cnt - 1], dy = yf - py[cnt - 1];
			double dist = hypot(dx, dy);
			double tpx, tpy, ang, ax[3], ay[3];
			ang = atan2(dy, dx) + M_PI / 3.0;
			tpx = dist * cos(ang);
			tpy = dist * sin(ang);
			if (debuglevel == 5) {
			    printf("s0 and dist = %lf %lf \n", s0, dist);
			}
			ax[0] = px[cnt - 1];
			ax[1] = xf;
			ax[2] = px[cnt - 1] + tpx;
			ay[0] = py[cnt - 1];
			ay[1] = yf;
			ay[2] = py[cnt - 1] + tpy;
			d1 = 0.0;
/*
   locate the element(s) that contain the corner nodes and compute the average
   depth in the element
*/
			for (j = 0; j < 3; j++) {
			    find_element(dgrid, ax[j], ay[j], &ind);
			    if (ind == -1) {
				if (debuglevel == 5) {
				    printf("Element not found in loc point\n");
				}
				d1 = d1 + celevmin;
			    } else {
				double dtmp = get_depth_element(dgrid, ind, ax[j], ay[j]);
				if (dtmp < celevmin) {
				    dtmp = celevmin;
				    errcnt++;
				}
				d1 = d1 + dtmp;
			    }
			}
			d1 *= 0.333333333333;
/*
   if the depth is less than the specified minimum, set the average depth
   to the minimum
*/
			if (d1 < cboundmin) {
			    d1 = cboundmin;
			}
/*
  check the criteria
*/
			A0 = comp_area(3, ax, ay);
			if (debuglevel == 5) {
			    printf("before adding point A0 = %lf d1 = %lf\n", A0, d1);
			}
			switch (crit) {
			case 0:
			    break;
			case 1:
			    accept = (A0 <= (M_PI * 9.8 * d1 * (dt * dt)) / (4 * cu * cu));
			    break;
			case 2:
			    accept = (A0 >= (M_PI * 9.8 * d1 * (dt * dt)) / (4 * cu * cu));
			    break;
			}
/*
   if accepted, then add the point to px[], py[]
*/
			if (accept) {
			    if (debuglevel == 5) {
				printf("adding point %lf %lf\n", xf, yf);
			    }
			    px[cnt] = xf;
			    py[cnt] = yf;
			    ip[cnt] = tcnt++;
			    xn = xf;
			    yn = yf;
			    setcolor(cnt % 4 + 2);
			    box(px[cnt - 1] + tpx, py[cnt - 1] + tpy);
			    my_move2(ax[0], ay[0]);
			    my_draw2(ax[1], ay[1]);
			    my_draw2(ax[2], ay[2]);
			    setcolor(1);
			    cnt++;
			    cur = j2;
			    goto out1;
			} else {
/*
   point rejected, shorten or lengthen the distance (depending on the
   selected criterion) and iterate.
*/
			    if (debuglevel == 5) {
				printf("Rejected point %lf %lf\n", xf, yf);
				setcolor(cnt % 4 + 2);
				solidbox(xf, yf);
			    }
			    s0 = sqrt(2 * A0 / (sin(M_PI / 3.0)));
			    switch (crit) {
			    case 0:
				break;
			    case 1:
				s0 = s0 - s0 * 0.1;
				break;
			    case 2:
				s0 = s0 + s0 * 0.1;
				break;
			    }
			    if (debuglevel == 5) {
				printf("new distance = %lf %lf %d\n", s0, A0, iter);
			    }
			    iter++;
			}
		    }		/* compute 3rd point */
		}		/* while iter */
/* exceeded iter */
		px[cnt] = xf;
		py[cnt] = yf;
		ip[cnt] = tcnt++;
		xn = xf;
		yn = yf;
		cnt++;
		cur = j2;
	out1:	;
	    }			/* while cur */
	}			/* if */
	if (cnt > 2) {
	    bno = nextboundary();
	    if (bno >= 0) {
		Allocate_boundary(bno, cnt, i == 0 ? 0 : 1);
		grid[gridno].boundaries[grid[gridno].nbounds] = bno;
		grid[gridno].nbounds++;
		for (j = 0; j < cnt; j++) {
		    boundary[bno].boundx[j] = px[j];
		    boundary[bno].boundy[j] = py[j];
		    boundary[bno].nodes[j] = ip[j];
		}
	    } else {
		errwin("Maximum number of boundaries (1000) exceeded");
	    }
	} else {
	    printf("Insufficient points for boundary\n");
	}
    }				/* for */
    if (grid[gridno].nbounds > 1) {
	ebound_defined = 1;
	ibound_defined = 1;
	update_fuzz_items();
    }
    if (grid[gridno].nbounds == 1) {
	ibound_defined = 0;
	ebound_defined = 1;
	update_fuzz_items();
    }
/*
    for (i = 0; i < grid[gridno].nbounds; i++) {
	bno = grid[gridno].boundaries[i];
	if (boundary[bno].bactive) {
	    fprintf(fp, "#cnt = %d\n", boundary[bno].nbpts);
	    for (j = 0; j < boundary[bno].nbpts; j++) {
		fprintf(fp, "%lf %lf\n", boundary[bno].boundx[j], boundary[bno].boundy[j]);
	    }

	    fprintf(fp, "&\n");
	}
    }
    fclose(fp);
*/
    draw_boundary(gridno);
    if (debuglevel == 5) {
	printf("returning %d\n", grid[gridno].nbounds);
    }
    return grid[gridno].nbounds;
}

void find_nearest_boundary_point(int *boundaries, int nb, double wx, double wy, int *bno, int *ind)
{
    int i, j, ib;
    double tmp, radius = 1e307;

    if (nb <= 0) {
	errwin("No boundaries");
	*bno = -1;
	*ind = -1;
	return;
    }
    for (i = 0; i < nb; i++) {
	ib = boundaries[i];
	for (j = 0; j < boundary[ib].nbpts; j++) {
	    tmp = hypot((wx - boundary[ib].boundx[j]), (wy - boundary[ib].boundy[j]));
	    if (tmp < radius) {
		radius = tmp;
		*ind = j;
		*bno = ib;
	    }
	}
    }
}

void find_nearest_boundary_point2(int bno, double wx, double wy, int *ind)
{
    int i, j, ib;
    double tmp, radius = 1e307;

    for (j = 0; j < boundary[bno].nbpts; j++) {
	tmp = hypot((wx - boundary[bno].boundx[j]), (wy - boundary[bno].boundy[j]));
	if (tmp < radius) {
	    radius = tmp;
	    *ind = j;
	}
    }
}

void do_check_gridboundary_intersect(void)
{
    char buf[256];
    sprintf(buf, "Performing edge and intersection checks:\n");
    stufftext(buf, 1);
    check_gridboundary_intersect(0);
    check_gridboundary_edges(0);
    stufftext("Done boundary checks\n", 0);
}

void check_gridboundary_edges(int gridno)
{
    char buf[256];
    int i, j, k, l, bno, n1, n2;
    double bx1, by1, bx2, by2, x1, x2, y1, y2;
    sprintf(buf, "Checking that every boundary edge is also an element edge:\n");
    stufftext(buf, 0);
    for (i = 0; i < grid[gridno].nbounds; i++) {
	bno = grid[gridno].boundaries[i];
	for (j = 0; j < boundary[bno].nbpts; j++) {
	    bx1 = boundary[bno].boundx[j];
	    by1 = boundary[bno].boundy[j];
	    bx2 = boundary[bno].boundx[(j + 1) % boundary[bno].nbpts];
	    by2 = boundary[bno].boundy[(j + 1) % boundary[bno].nbpts];
	    for (k = 0; k < grid[gridno].nmel; k++) {
		for (l = 0; l < 3; l++) {
		    n1 = grid[gridno].icon[k].nl[l];
		    n2 = grid[gridno].icon[k].nl[(l + 1) % 3];
		    x1 = grid[gridno].xord[n1];
		    y1 = grid[gridno].yord[n1];
		    x2 = grid[gridno].xord[n2];
		    y2 = grid[gridno].yord[n2];
		    if ((bx1 == x1 && by1 == y1) && (bx2 == x2 && by2 == y2)) {
			goto out;
		    }
		}
	    }
	    sprintf(buf, "In boundary %d, edge %d (%d, %d) is not an element edge\n",
		    i + 1, j + 1, j + 1, (j + 1) % boundary[bno].nbpts + 1);
	    stufftext(buf, 0);
    out:    ;
	}
    }
    stufftext("Done edge checks\n", 0);
}

void check_gridboundary_intersect(int gridno)
{
    char buf[256];
    int i, j, b1, b2;
    sprintf(buf, "Checking for intersecting boundary edges:\n");
    stufftext(buf, 0);
    for (i = 0; i < grid[gridno].nbounds; i++) {
	b1 = grid[gridno].boundaries[i];
	for (j = i; j < grid[gridno].nbounds; j++) {
	    b2 = grid[gridno].boundaries[j];
	    if (check_boundary_intersect(b1, b2)) {
		sprintf(buf, "boundaries %d and %d do not intersect\n", i + 1, j + 1);
		stufftext(buf, 0);
	    }
	}
    }
    stufftext("Done intersection checks\n", 0);
}

int check_boundary_intersect(int b1, int b2)
{
    char buf[256];
    int i, j, ib, iprev, inext;
    double tmp, radius = 1e307;
    double xk, yk, xl, yl, xm, ym, xn, yn, x, y;

    for (j = 0; j < boundary[b1].nbpts; j++) {
	xk = boundary[b1].boundx[j];
	yk = boundary[b1].boundy[j];
	xl = boundary[b1].boundx[(j + 1) % boundary[b1].nbpts];
	yl = boundary[b1].boundy[(j + 1) % boundary[b1].nbpts];
	if (b1 == b2) {
	    if (j == 0) {
		iprev = boundary[b2].nbpts - 1;
	    } else {
		iprev = (j - 1);
	    }
	    if (j == boundary[b1].nbpts - 1) {
		inext = 0;
	    } else {
		inext = (j + 1);
	    }
	}
	for (i = 0; i < boundary[b2].nbpts; i++) {
	    if (b1 == b2) {
		if (i == iprev || i == j || i == inext) {
		    continue;
		}
	    }
	    xm = boundary[b2].boundx[i];
	    ym = boundary[b2].boundy[i];
	    xn = boundary[b2].boundx[(i + 1) % boundary[b2].nbpts];
	    yn = boundary[b2].boundy[(i + 1) % boundary[b2].nbpts];
	    if (intersect_lines(xk, yk, xl, yl, xm, ym, xn, yn, &x, &y)) {
		sprintf(buf, "Boundary segments intersect %d:%d %d:%d at %lf %lf\n", b1 + 1, i + 1, b2 + 1, j + 1, x, y);
		stufftext(buf, 1);
		return 0;
	    }
	}
    }
    return 1;
}

int intersect_lines(double xk, double yk, double xl, double yl,
	   double xm, double ym, double xn, double yn, double *x, double *y)
{
    double xlk, ylk, xnm, ynm, xmk, ymk, s, t, det, detinv;
    xlk = xl - xk;
    ylk = yl - yk;
    xnm = xn - xm;
    ynm = yn - ym;
    xmk = xm - xk;
    ymk = ym - yk;
    det = xnm * ylk - ynm * xlk;
    if (fabs(det) < 1.0e-11) {
	return 0;
    }
    detinv = 1.0 / det;
    s = (xnm * ymk - ynm * xmk) * detinv;
    t = (xlk * ymk - ylk * xmk) * detinv;
    if ((s < 0.0 || s > 1.0) || (t < 0.0 || t > 1.0)) {
	return 0;
    }
    *x = xk + xlk * s;
    *y = yk + ylk * s;
    return 1;
}
