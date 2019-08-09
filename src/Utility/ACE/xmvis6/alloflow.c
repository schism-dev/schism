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
 * alloflow.c - allocate a flow field
 *
 */

#ifndef lint
static char RCSid[] = "$Id: alloflow.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

int Allocate_flowf(int flowno, int nmnp, int nfreq, int use_nodes)
{
    int i;

    flowf[flowno].active = ON;
    flowf[flowno].nfreq = nfreq;
    if (!use_nodes) {
	flowf[flowno].x = (double *) calloc(nmnp, sizeof(double));
	flowf[flowno].y = (double *) calloc(nmnp, sizeof(double));
    } else {
	flowf[flowno].x = NULL;
	flowf[flowno].y = NULL;
    }
    return 1;
}

void Free_flowf(int flowfno)
{
    int i;

    if (flowf[flowfno].x != NULL) {
	free(flowf[flowfno].x);
    }
    if (flowf[flowfno].y != NULL) {
	free(flowf[flowfno].y);
    }
    for (i = 0; i < flowf[flowfno].nfreq; i++) {
	if (flowf[flowfno].elamp[i] != NULL) {
	    free(flowf[flowfno].elamp[i]);
	}
	if (flowf[flowfno].elphase[i] != NULL) {
	    free(flowf[flowfno].elphase[i]);
	}
	if (flowf[flowfno].ampx[i] != NULL) {
	    free(flowf[flowfno].ampx[i]);
	}
	if (flowf[flowfno].ampy[i] != NULL) {
	    free(flowf[flowfno].ampy[i]);
	}
	if (flowf[flowfno].phax[i] != NULL) {
	    free(flowf[flowfno].phax[i]);
	}
	if (flowf[flowfno].phay[i] != NULL) {
	    free(flowf[flowfno].phay[i]);
	}
    }
    flowf[flowfno].nfreq = 0;
    flowf[flowfno].active = OFF;
}

/*
 * ADCIRC
 */
int Allocate_flowt(int flowno, int nmnp, int nfreq, int use_nodes)
{
    int i;
    flowt[flowno].active = 1;
    if (!use_nodes) {
	flowt[flowno].x = (double *) calloc(nmnp, sizeof(double));
	flowt[flowno].y = (double *) calloc(nmnp, sizeof(double));
    } else {
	flowt[flowno].x = NULL;
	flowt[flowno].y = NULL;
    }
    return 1;
}

void Free_flowt(int flowno)
{
    int i;
    if (flowt[flowno].x != NULL) {
	free(flowt[flowno].x);
    }
    if (flowt[flowno].y != NULL) {
	free(flowt[flowno].y);
    }
    for (i = 0; i < flowt[flowno].nsteps; i++) {
	if (flowt[flowno].f) {
	    cxfree(flowt[flowno].f[i].u);
	    cxfree(flowt[flowno].f[i].v);
	    cxfree(flowt[flowno].f[i].e);
	}
    }
    if (flowt[flowno].f) {
	cxfree(flowt[flowno].f);
	flowt[flowno].f = NULL;
    }
    flowt[flowno].active = OFF;
}

/*
 * Time histories of flow
*/
int Allocate_flowh(int flowno, int nmnp, int use_nodes)
{
    int i;

    flowh[flowno].active = 1;
    if (!use_nodes) {
	flowh[flowno].x = (double *) calloc(nmnp, sizeof(double));
	flowh[flowno].y = (double *) calloc(nmnp, sizeof(double));
    } else {
	flowh[flowno].x = NULL;
	flowh[flowno].y = NULL;
    }
    return 1;
}

void Free_flowh(int flowno)
{
    int i;
    if (flowh[flowno].x != NULL) {
	free(flowh[flowno].x);
    }
    if (flowh[flowno].y != NULL) {
	free(flowh[flowno].y);
    }
    for (i = 0; i < flowh[flowno].nsteps; i++) {
	if (flowh[flowno].f) {
	    cxfree(flowh[flowno].f[i].u);
	    cxfree(flowh[flowno].f[i].v);
	    cxfree(flowh[flowno].f[i].e);
	}
    }
    if (flowh[flowno].f) {
	cxfree(flowh[flowno].f);
	flowh[flowno].f = NULL;
    }
    flowh[flowno].active = 0;
}
