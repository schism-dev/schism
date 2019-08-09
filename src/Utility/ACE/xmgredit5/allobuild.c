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
 * allobuild.c - allocate/delete/add to a set of build points
 *
 */

#ifndef lint
static char RCSid[] = "$Id: allobuild.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "globals.h"

/* allobuild.c */
int Allocate_build(int bno, int npts);
void Reallocate_build(int bno, int npts);
void add_build(int bno, double wx, double wy, double d);
void add_builds(int bno, int n, double *wx, double *wy, double *d);
void Free_build(int bno);
void Free_builds(void);
void Init_builds(void);
int nextbuild(void);
void init_qadd_build(int buildno);
void flush_qadd_build(int buildno);
void qadd_build(int buildno, double x, double y, double d);

/*
 * Allocate storage for build points struct bno of length npts
 */
int Allocate_build(int bno, int npts)
{
    build[bno].nbuild = npts;
    build[bno].buildactive = 1;
    build[bno].bx = (double *) calloc(npts, sizeof(double));
    build[bno].by = (double *) calloc(npts, sizeof(double));
    build[bno].db = (double *) calloc(npts, sizeof(double));
    if (build[bno].nbuild) {
	build_defined = 1;
	update_fuzz_items();
    } else {
	build_defined = 0;
	update_fuzz_items();
    }
    return 1;
}

/*
 * Resize a build bpoints structure
 */
void Reallocate_build(int bno, int npts)
{
    build[bno].nbuild = npts;
    build[bno].bx = (double *) realloc(build[bno].bx, npts * sizeof(double));
    build[bno].by = (double *) realloc(build[bno].by, npts * sizeof(double));
    build[bno].db = (double *) realloc(build[bno].db, npts * sizeof(double));
    if (build[bno].nbuild) {
	build_defined = 1;
	update_fuzz_items();
    } else {
	build_defined = 0;
	update_fuzz_items();
    }
}

/*
 * Add a single build point
 */
void add_build(int bno, double wx, double wy, double d)
{
    int npts;

    if (build[bno].nbuild == 0) {
	Free_build(bno);
	Allocate_build(bno, 1);
	build[bno].bx[0] = wx;
	build[bno].by[0] = wy;
	build[bno].db[0] = d;
    } else {
	build[bno].nbuild++;
	npts = build[bno].nbuild;
	build[bno].bx = (double *) realloc(build[bno].bx, npts * sizeof(double));
	build[bno].by = (double *) realloc(build[bno].by, npts * sizeof(double));
	build[bno].db = (double *) realloc(build[bno].db, npts * sizeof(double));
	build[bno].bx[npts - 1] = wx;
	build[bno].by[npts - 1] = wy;
	build[bno].db[npts - 1] = d;
    }
    if (build[bno].nbuild) {
	build_defined = 1;
	update_fuzz_items();
    } else {
	build_defined = 0;
	update_fuzz_items();
    }
}

/*
 * Add an array of build points
 */
void add_builds(int bno, int n, double *wx, double *wy, double *d)
{
    int i;
    int npts = build[bno].nbuild;

    if (npts) {
	build[bno].nbuild += n;
	build[bno].bx = (double *) realloc(build[bno].bx, (npts + n) * sizeof(double));
	build[bno].by = (double *) realloc(build[bno].by, (npts + n) * sizeof(double));
	build[bno].db = (double *) realloc(build[bno].db, (npts + n) * sizeof(double));
	for (i = npts; i < npts + n; i++) {
	    build[bno].bx[i] = wx[i - npts];
	    build[bno].by[i] = wy[i - npts];
	    build[bno].db[i] = d[i - npts];
	}
    } else {
	build[bno].nbuild = n;
	build[bno].bx = (double *) calloc((npts + n), sizeof(double));
	build[bno].by = (double *) calloc((npts + n), sizeof(double));
	build[bno].db = (double *) calloc((npts + n), sizeof(double));
	for (i = 0; i < n; i++) {
	    build[bno].bx[i] = wx[i];
	    build[bno].by[i] = wy[i];
	    build[bno].db[i] = d[i];
	}
    }
    if (build[bno].nbuild) {
	build_defined = 1;
	update_fuzz_items();
    } else {
	build_defined = 0;
	update_fuzz_items();
    }
}

/*
 * Free a build points structure
 */
void Free_build(int bno)
{
    build[bno].nbuild = 0;
    build[bno].buildactive = 0;
    if (build[bno].bx != NULL) {
	free(build[bno].bx);
    }
    if (build[bno].by != NULL) {
	free(build[bno].by);
    }
    if (build[bno].db != NULL) {
	free(build[bno].db);
    }
    build[bno].bx = build[bno].by = build[bno].db = NULL;
    build_defined = 0;
    update_fuzz_items();
}

/*
 * Free all build points structures
 */
void Free_builds(void)
{
    int i;

    for (i = 0; i < MAXBUILD; i++) {
	Free_build(i);
    }
}

/*
 * Initialize all build points structures
 */
void Init_builds(void)
{
    int i;

    for (i = 0; i < MAXBUILD; i++) {
	build[i].nbuild = 0;
	build[i].buildactive = 0;
	build[i].bx = NULL;
	build[i].by = NULL;
	build[i].db = NULL;
    }
    build_defined = 0;
    update_fuzz_items();
}

/*
 * Get the next available build points structure
 */
int nextbuild(void)
{
    int i;

    for (i = 0; i < MAXBUILD; i++) {
	if (build[i].buildactive == 0) {
	    return i;
	}
    }
    return -1;
}

/*
 * The follow was implemented for performance reasons, the purpose
 * is to queue build points to be added, then at BUFSIZE, call the
 * bulk build points adder routine.
 */
#define BUFSIZE 1000

static double *qx, *qy, *qd;
static int qn;

void init_qadd_build(int buildno)
{
    if (qx) {
	free(qx);
	free(qy);
	free(qd);
    }
    qn = 0;
    qx = (double *) calloc(BUFSIZE, sizeof(double));
    qy = (double *) calloc(BUFSIZE, sizeof(double));
    qd = (double *) calloc(BUFSIZE, sizeof(double));
}

void flush_qadd_build(int buildno)
{
    add_builds(buildno, qn, qx, qy, qd);
    if (qx) {
	free(qx);
	free(qy);
	free(qd);
	qx = NULL;
    }
}

void qadd_build(int buildno, double x, double y, double d)
{
    qx[qn] = x;
    qy[qn] = y;
    qd[qn] = d;
    qn++;
    if (qn % BUFSIZE == 0) {
	qx = (double *) realloc(qx, (qn + BUFSIZE) * sizeof(double));
	qy = (double *) realloc(qy, (qn + BUFSIZE) * sizeof(double));
	qd = (double *) realloc(qd, (qn + BUFSIZE) * sizeof(double));
    }
}
