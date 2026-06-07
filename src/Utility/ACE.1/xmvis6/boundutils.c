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
 * utilities for manipulating boundaries
 *
 */

#ifndef lint
static char RCSid[] = "$Id: boundutils.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";

#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

static int nb;

#define MAXB 100000

void getboundary(int gridno, int (*clist)[2], int *nb)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[10], v, nmnp;
    int mnode;
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
	if (grid[gridno].icon[i].type != 0 && grid[gridno].icon[i].type != 2 && grid[gridno].icon[i].type != 5 && grid[gridno].icon[i].type != 1) {
	    switch (grid[gridno].icon[i].type) {
	    case 3:
		for (j = 0; j < 3; j++) {
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
	    case 4:
		for (j = 0; j < 4; j++) {
		    el[j] = grid[gridno].icon[i].nl[j];
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
	    }
	    for (j = 0; j < grid[gridno].icon[i].nn; j++) {
/* search for edge */
		cnt1 = 0;
		cnt2 = 0;
		f1 = 0;
		f2 = 0;
		n1 = el[j];
		n2 = el[(j + 1) % grid[gridno].icon[i].nn];

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
		    printf("Error: duplicate edge, send mail to ace_bugs\n");
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

void load_boundary(int gridno, int (*clist)[2])
{
    int start, bno, done, n = 0, *outb;
    int i, j, cnt, next;
    double uniform(), tmpx, tmpy;
    FILE *fp;
    Boundary b;

    outb = (int *) malloc(MAXB * sizeof(int));

    outb[0] = clist[0][0];
    next = clist[0][1];
    clist[0][1] = -1;
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
	bno = n;
	n++;
	Allocate_boundaries(gridno, n);
	if (n == 1) {
	    Allocate_boundary(gridno, bno, cnt - 1, 0);
	} else if (n > 1) {
	    Allocate_boundary(gridno, bno, cnt - 1, 1);
	}
	for (i = 0; i < cnt - 1; i++) {
	    grid[gridno].boundaries[bno].boundx[i] = grid[gridno].xord[outb[i]];
	    grid[gridno].boundaries[bno].boundy[i] = grid[gridno].yord[outb[i]];
	    grid[gridno].boundaries[bno].nodes[i] = outb[i];
	}
	cnt = 0;
	done = 1;
	for (i = 0; i < nb; i++) {
	    if (clist[i][0] != -1 && clist[i][1] != -1) {
		outb[cnt++] = clist[i][0];
		next = clist[i][1];
		clist[i][1] = -1;
		done = 0;
		goto getout2;
	    }
	}
      getout2:;
    }
  giveup:;
    grid[gridno].nbounds = n;
    free(outb);
}

void compute_boundary(int gridno)
{
    int i, clist[MAXB][2];

    nb = 0;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	Free_boundary(gridno, i);
    }
    grid[gridno].nbounds = 0;
    getboundary(gridno, clist, &nb);
    load_boundary(gridno, clist);
}

void compute_nearest_point(int gridno, double xn, double yn, double *xi, double *yi)
{
    int i, j, ib;
    double t, mdist, x1, y1, x2, y2, tmp, a, b, c1, c2, dist;
    Boundary boundary;

    for (i = 0; i < grid[gridno].nbounds; i++) {
	boundary = grid[gridno].boundaries[i];
	for (j = 0; j < boundary.nbpts; j++) {
	    x1 = boundary.boundx[j];
	    y1 = boundary.boundy[j];
	    x2 = boundary.boundx[(j + 1) % boundary.nbpts];
	    y2 = boundary.boundy[(j + 1) % boundary.nbpts];
	    c1 = x1 - xn;
	    c2 = y1 - yn;
	    a = x2 - x1;
	    b = y2 - y1;
	    tmp = a * a + b * b;
	    t = -((a * c1) + (b * c2)) / tmp;
	    dist = hypot((x1 - xn + t * a), (y1 - yn + t * b));
	    if (j == 0) {
		mdist = dist;
	    } else {
		if (t >= 0.0 && t <= 1.0) {
		} else {
		}
	    }
	    if (t >= 0.0 && t <= 1.0) {
		*xi = x1 + t * a;
		*yi = y1 + t * a;
	    }
	}
    }
}

int compute_intersection(int gridno, double xl, double yl, double xk, double yk, double xm, double ym, double xn, double yn, double *xi, double *yi)
{
    double t1, t2, det, detinv, xlk, ylk, xnm, ynm, xmk, ymk, dist;

    xlk = xl - xk;
    ylk = yl - yk;
    xnm = xn - xm;
    ynm = yn - ym;
    xmk = xm - xk;
    ymk = ym - yk;
    det = xnm * ylk - ynm * xlk;
    if (fabs(det) < 1.0e-14) {
	return 0;
    }
    detinv = 1.0 / det;
    t1 = (xnm * ymk - ynm * xmk) * detinv;
    t2 = (xlk * ymk - ylk * xmk) * detinv;
    if (t1 < 0.0 || t1 > 1.0 || t2 < 0.0 || t2 > 1.0) {
	return 0;
    }
    *xi = xk + xlk * t1;
    *yi = yk + ylk * t1;
    return 1;
}
