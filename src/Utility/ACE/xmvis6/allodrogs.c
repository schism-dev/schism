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
 * allodrogs.c - allocate a single time step for drogues
 *
 */

#ifndef lint
static char RCSid[] = "$Id: allodrogs.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

void Free_drogues(int drno);

int Allocate_drogstep(int drno, int nt, double time, int npts)
{
    drogues[drno].p[nt].npts = npts;
    drogues[drno].p[nt].x = (double *) malloc(npts * sizeof(double));
    drogues[drno].p[nt].y = (double *) malloc(npts * sizeof(double));
    drogues[drno].p[nt].z = (double *) malloc(npts * sizeof(double));
    drogues[drno].p[nt].drnum = (int *) malloc(npts * sizeof(int));
    drogues[drno].p[nt].time = time;
    return 1;
}

int Allocate_drogues(int drno, int nsteps)
{
    Free_drogues(drno);
    drogues[drno].active = ON;
    drogues[drno].nsteps = nsteps;
    drogues[drno].p = (Particles *) malloc(nsteps * sizeof(Particles));
    return 1;
}

void Free_drogues(int drno)
{
    int i;
    for (i = 0; i < drogues[drno].nsteps; i++) {
	if (drogues[drno].p[i].x != NULL) {
	    drogues[drno].p[i].x = NULL;
	}
	if (drogues[drno].p[i].y != NULL) {
	    free(drogues[drno].p[i].y);
	}
	if (drogues[drno].p[i].z != NULL) {
	    free(drogues[drno].p[i].z);
	}
	if (drogues[drno].p[i].drnum != NULL) {
	    free(drogues[drno].p[i].drnum);
	}
    }
    if (drogues[drno].p != NULL) {
	free(drogues[drno].p);
    }
    drogues[drno].p = NULL;
    drogues[drno].nsteps = 0;
    drogues[drno].active = OFF;
}

void elim_drogs(int drno)
{
    int n = drogues[drno].nsteps;
    int final_pts = drogues[drno].p[n - 1].npts;
    if (final_pts == 0) {
	errwin("No survivors! operation canceled");
	return;
    }
    if (final_pts == drogues[drno].p[0].npts) {
	errwin("All drogues survived");
	return;
    } else {
	errwin("Can't allocate enough memory, operation canceled");
	return;
    }
}

void copy_drogs(int drfrom, int drto)
{
}
