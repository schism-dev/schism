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
 *  flowio.c - read/write flow fields
 *
 */

#ifndef lint
static char RCSid[] = "$Id: flowio.c,v 1.21 2008/04/10 21:39:57 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"
#include "elio.h"

/*
 * Read a frequency domain flow file
 */
int readflow(int flowno, int gridno, char *fname)
{
    FILE *fp;
    int i, j, slen, ifl, ni, nfreq, nread = 0, nmnp = grid[gridno].nmnp;
    int imin, imax, touflag = 0;
    double omega, pcor;
    double e, emin = BIG, emax = MBIG;
    double temin = 0.0, temax = 0.0;
    char buf[257], s[5];
    double missing;
    int itmp;

    if ((fp = fopen(fname, "r")) == NULL) {
	return 0;
    }

    slen = strlen(fname);
    s[0] = fname[slen - 3];
    s[1] = fname[slen - 2];
    s[2] = fname[slen - 1];
    s[3] = 0;
    fgets(buf, 256, fp);
    fgets(buf, 256, fp);
    sscanf(buf, "%d", &nfreq);
    if (nfreq == 0) {
	errwin("Number of frequencies = 0\n");
	return 0;
    }
    flowf[flowno].nread = nfreq;
    Allocate_flowf(flowno, grid[gridno].nmnp, nfreq, 1, 0);
    flowf[flowno].grid = gridno;
    for (j = 0; j < flowf[flowno].nfreq; j++) {
	ifl = 0;
	if (fgets(buf, 256, fp) == NULL) {
	    goto branch;
	}
	if (strcmp(s, "tou")) {	/* a tct file */
	    sscanf(buf, "%lf %d", &omega, &ifl);
	    fgets(buf, 256, fp);
	    buf[strlen(buf) - 1] = 0;
	    flowf[flowno].sign_of_phase = 1;
	    strcpy(flowf[flowno].freqname[nread], buf);
	    missing = BIG;
	} else {		/* a tou file */
	    buf[strlen(buf) - 1] = 0;
	    strcpy(flowf[flowno].freqname[nread], buf);
	    fgets(buf, 256, fp);
	    sscanf(buf, "%lf", &omega);
	    ifl = 1;
	    touflag = 1;
	    flowf[flowno].sign_of_phase = -1;
	    missing = 8888.0;
	}
	if (ifl) {		/* got a frequency */
	    flowf[flowno].omega[nread] = omega;
	    buf[strlen(buf) - 1] = 0;
	    flowf[flowno].elamp[nread] = (double *) calloc(nmnp, sizeof(double));
	    flowf[flowno].elphase[nread] = (double *) calloc(nmnp, sizeof(double));
	    flowf[flowno].ampx[nread] = (double *) calloc(nmnp, sizeof(double));
	    flowf[flowno].ampy[nread] = (double *) calloc(nmnp, sizeof(double));
	    flowf[flowno].phax[nread] = (double *) calloc(nmnp, sizeof(double));
	    flowf[flowno].phay[nread] = (double *) calloc(nmnp, sizeof(double));

/* read elevations */
	    for (i = 0; i < grid[gridno].nmnp; i++) {
		fgets(buf, 256, fp);
		sscanf(buf, "%d %lf %lf", &ni, &flowf[flowno].elamp[nread][i], &flowf[flowno].elphase[nread][i]);
	    }
	    minmax2(flowf[flowno].elamp[nread], nmnp, missing, &flowf[flowno].ampmin[nread], &flowf[flowno].ampmax[nread], &imin, &imax);
	    minmax2(flowf[flowno].elphase[nread], nmnp, missing, &flowf[flowno].phamin[nread], &flowf[flowno].phamax[nread], &imin, &imax);
	    if (temax < flowf[flowno].ampmax[nread]) {
		temax = flowf[flowno].ampmax[nread];
	    }
/* read flow */
	    for (i = 0; i < grid[gridno].nmnp; i++) {
		fgets(buf, 256, fp);
		sscanf(buf, "%d %lf %lf %lf %lf", &ni, &flowf[flowno].ampx[nread][i], &flowf[flowno].phax[nread][i], &flowf[flowno].ampy[nread][i], &flowf[flowno].phay[nread][i]);
	    }
	    nread++;
	}
    }
  branch:;
    flowf[flowno].emin = -temax;
    flowf[flowno].emax = temax;
    fclose(fp);
    flowf[flowno].nfreq = nread;
    flowf[flowno].nread = nfreq;
    flowf[flowno].active = ON;
    return (1);
}

/*
 * filet = 0 - both elevations and flows
 * filet = 1 - elevations
 * filet = 2 - flows
 */
int readflowt2d(int flowno, int filet, char *gfile, char *fname, int keep)
{
    FILE *fp;
    char *s;
    int scnt = 0, i, j, k, itmp, ni, nread = 0, nmnp;
    double omega, e, emin, emax;
    double temin, temax;
    char buf[1024];
    Grid *g;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readflow, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    for (j = 0; j < flowt[flowno].nsteps; j++) {
	if (flowt[flowno].f) {
	    if (flowt[flowno].f[j].e) {
		cxfree(flowt[flowno].f[j].e);
		flowt[flowno].f[j].e = NULL;
	    }
	    if (flowt[flowno].f[j].u) {
		cxfree(flowt[flowno].f[j].u);
		flowt[flowno].f[j].u = NULL;
	    }
	    if (flowt[flowno].f[j].v) {
		cxfree(flowt[flowno].f[j].v);
		flowt[flowno].f[j].v = NULL;
	    }
	}
    }


    g = &flowt[flowno].g;
    ReadGrid(gfile, g);

    nmnp = g->nmnp;

    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);

    flowt[flowno].nsteps = 0;
    if (flowt[flowno].f) {
	cxfree(flowt[flowno].f);
	flowt[flowno].f = NULL;
    }
    fgets(buf, 256, fp);
    fgets(buf, 256, fp);
    sscanf(buf, "%d", &flowt[flowno].nsteps);
    flowt[flowno].npts = flowt[flowno].g.nmnp;
    /*flowt[flowno].f = (Flow_gen *) calloc(flowt[flowno].nsteps, sizeof(Flow_gen)); */
/*
 * in case the number of time steps in the header do not match the number of
 * time steps in the file, just allocate what is needed
 */
    flowt[flowno].f = (Flow_gen *) malloc(sizeof(Flow_gen));
    for (j = 0; j < flowt[flowno].nsteps; j++) {
	switch (filet) {
	case 0:
	    flowt[flowno].f[j].e = (float *) calloc(nmnp, sizeof(float));
	    flowt[flowno].f[j].u = NULL;
	    flowt[flowno].f[j].v = NULL;
	    break;
	case 1:
	    flowt[flowno].f[j].u = (float *) calloc(nmnp, sizeof(float));
	    flowt[flowno].f[j].v = (float *) calloc(nmnp, sizeof(float));
	    flowt[flowno].f[j].e = NULL;
	    break;
	case 2:
	    flowt[flowno].f[j].u = (float *) calloc(nmnp, sizeof(float));
	    flowt[flowno].f[j].v = (float *) calloc(nmnp, sizeof(float));
	    flowt[flowno].f[j].e = (float *) calloc(nmnp, sizeof(float));
	    break;
	}
	if (fgets(buf, 1024, fp) == NULL) {
	    goto bustout;
	}
	sscanf(buf, "%lf", &flowt[flowno].f[j].time);
/*
   printf("ADCIRC step # %d time = %lf\n", j, flowt[flowno].f[j].time);
 */
	for (i = 0; i < flowt[flowno].npts; i++) {
	    if (fgets(buf, 1024, fp) == NULL) {
		goto bustout;
	    }
	    switch (filet) {
	    case 0:
		sscanf(buf, "%*d %f", &flowt[flowno].f[j].e[i]);
		break;
	    case 1:
		sscanf(buf, "%*d %f %f", &flowt[flowno].f[j].u[i], &flowt[flowno].f[j].v[i]);
		break;
	    case 2:
		sscanf(buf, "%*d %f %f %f", &flowt[flowno].f[j].e[i], &flowt[flowno].f[j].u[i], &flowt[flowno].f[j].v[i]);
		break;
	    }
	    if (filet == 0 || filet == 2) {
		e = flowt[flowno].f[j].e[i];
		if (i == 0) {
		    emin = emax = e;
		    if (j == 0) {
			temin = temax = e;
		    }
		}
		if (e < temin) {
		    temin = e;
		}
		if (e > temax) {
		    temax = e;
		}
		if (e < emin) {
		    emin = e;
		}
		if (e > emax) {
		    emax = e;
		}
		flowt[flowno].f[j].emin = emin;
		flowt[flowno].f[j].emax = emax;
	    }
	}
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (j + 2) * sizeof(Flow_gen));
	scnt++;
    }
  bustout:;			/* failed to read nsteps timesteps */
    if (scnt != flowt[flowno].nsteps) {
	sprintf(buf, "File claims %d timesteps, read only %d,\n resetting number of steps\n", flowt[flowno].nsteps, scnt);
	errwin(buf);
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, scnt * sizeof(Flow_gen));
	flowt[flowno].nsteps = scnt;
    }
    if (filet == 0 || filet == 2) {
	flowt[flowno].emin = temin;
	flowt[flowno].emax = temax;
    }
    flowt[flowno].type = filet;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    if (flowt[flowno].nsteps > 0) {
	flowt[flowno].step = flowt[flowno].f[1].time - flowt[flowno].f[0].time;
    } else {
	flowt[flowno].step = 0.0;
    }
    flowt[flowno].active = ON;
    fclose(fp);
    return (1);
}

int readflowh2d(int flowno, int gridno, char *fname)
{
    FILE *fp;
    char *s;
    int npts, i, j, k, itmp, ifl, ni, nfreq, nread = 0, nmnp = flowt[flowno].g.nmnp;
    double omega, e, emin = 1e307, emax = -1e307;
    double temin = 1e307, temax = -1e307;
    char buf[257];

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readflow, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    for (j = 0; j < flowh[flowno].nsteps; j++) {
	if (flowh[flowno].f) {
	    cxfree(flowh[flowno].f[j].u);
	    cxfree(flowh[flowno].f[j].v);
	    cxfree(flowh[flowno].f[j].e);
	}
    }
    if (flowh[flowno].f) {
	cxfree(flowh[flowno].f);
	flowh[flowno].f = NULL;
    }
    fgets(buf, 256, fp);
    fgets(buf, 256, fp);
    sscanf(buf, "%d %d", &flowh[flowno].npts, &flowh[flowno].nsteps);
    flowh[flowno].f = (Flow_gen *) calloc(flowh[flowno].nsteps, sizeof(Flow_gen));
    flowh[flowno].x = (double *) calloc(flowh[flowno].npts, sizeof(double));
    flowh[flowno].y = (double *) calloc(flowh[flowno].npts, sizeof(double));
    for (j = 0; j < flowh[flowno].nsteps; j++) {
	flowh[flowno].f[j].e = (float *) calloc(nmnp, sizeof(float));
	flowh[flowno].f[j].u = (float *) calloc(nmnp, sizeof(float));
	flowh[flowno].f[j].v = (float *) calloc(nmnp, sizeof(float));
    }
    for (j = 0; j < flowh[flowno].npts; j++) {
	fgets(buf, 1024, fp);
	sscanf(buf, "%lf %lf", &flowh[flowno].x[j], &flowh[flowno].y[j]);
    }
    for (j = 0; j < flowh[flowno].nsteps; j++) {
	fgets(buf, 1024, fp);
	sscanf(buf, "%lf", &flowh[flowno].f[j].time);
	for (i = 0; i < flowh[flowno].npts; i++) {
	    fgets(buf, 1024, fp);
	    sscanf(buf, "%*d %lf %lf %lf", &flowh[flowno].f[j].e[i], &flowh[flowno].f[j].u[i], &flowh[flowno].f[j].v[i]);
	}
    }
    flowh[flowno].start = flowh[flowno].f[0].time;
    flowh[flowno].stop = flowh[flowno].f[flowh[flowno].nsteps - 1].time;
    flowh[flowno].active = ON;
    fclose(fp);
    return (1);
}

int ReadElcirc(int flowno, char *fname, int level, int start, int stop, int skip, int missing, double mval, int append)
{
    FILE *fp;
    int ftype;
    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcirc(): Unable to open file %s\n", fname);
	return 0;
    }
    ftype = ElioGetFileType(fp);
    fclose(fp);
    if (ftype == -1) {
	errwin("ReadElcirc(): Unable to read old format files with this version of ACE/vis");
	return 1;
    } else {
	return ReadElcircNew(flowno, fname, level, start, stop, skip, missing, mval, append);
    }
    return 0;
}

int ReadElcircNew(int flowno, char *fname, int level, int start, int stop, int skip, int missing, double mval, int append)
{
    FILE *fp;
    int *readdata;
    int cnt, scnt = 0, i, j, k, itmp;
    double emin, emax;
    double temin, temax;
    double umin, umax;
    double tumin, tumax;
    double vmin, vmax;
    double tvmin, tvmax;
    float *tmpe, *data;
    char buf[1024];
    int first, nsteps, err, ftype;
    long offset;
    ElcircHeader h;
    ElcircTimeStep t;
    Grid *g;

    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to open file %s\n", fname);
	return 0;
    }
    if (!append) {
	Free_flowt(flowno);
    }

    ftype = ElioGetFileType(fp);
/* Get the header */
    if (err = ElioGetHeader(fname, &h)) {
	fprintf(stderr, "ReadElcircNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
/*
    ElioPrintHeader(h);
*/

/*
   set the time.
   12/02/2002 00:00:00 PST
   YYYYMMDDHHMMSS
*/
    if (!append) {
	settime_info_start(cg, h.start_time);
	g = &flowt[flowno].g;
	GetElcircGrid(&h, g);
    }

    nsteps = ElioGetNStepsInFile(fname, &h);
/*
    printf("ReadElcircNew(): Number of steps in header = %d, really = %d\n", h.nsteps, scnt);
*/
    if (start < 0) {
	start = 0;
	stop = nsteps - 1;
	skip = 1;
    } else {
	if (start >= nsteps) {
	}
	if (stop >= nsteps) {
	    stop = nsteps - 1;
	}
    }
    scnt = 0;
    readdata = (int *) malloc(nsteps * sizeof(int));
    for (i = 0; i < nsteps; i++) {
	readdata[i] = 0;
    }
    for (i = start; i <= stop; i += skip) {
	readdata[i] = 1;
	scnt++;
    }
    if (scnt == 0) {
	errwin("ReadElcircNew(): No complete time step available, cancelling read");
	flowt[flowno].active = OFF;
	free(readdata);
	return 0;
    }

    if (append) {
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (flowt[flowno].nsteps + scnt) * sizeof(Flow_gen));
	scnt = flowt[flowno].nsteps;
    } else {
	flowt[flowno].npts = h.np;
	flowt[flowno].nsteps = scnt;
	flowt[flowno].f = (Flow_gen *) malloc(scnt * sizeof(Flow_gen));
	scnt = 0;
    }

    if (debuglevel == 13) {
	printf("Trying to read %d steps:\n", flowt[flowno].nsteps);
    }

    if (flowt[flowno].f == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to allocate memory for reading timesteps\n");
	free(readdata);
	return 0;
    }

/* conversion array */
    data = (float *) malloc(h.nvrt * h.np * h.ivs * 4);
    if (data == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to allocate memory for reading data, %d %d %d\n", h.nvrt, h.np, h.ivs);
	Free_flowt(flowno);
	free(readdata);
	return 0;
    }

/* pre-allocate memory for timestep struct */
    ElioAllocateTimeStep(&h, &t);

    first = 1;
    for (j = 0; j < nsteps; j++) {
	if (readdata[j] == 0) {
	    continue;
	}
	if (h.ivs == 2) {	/* Vector data */
	    flowt[flowno].f[scnt].e = NULL;
	    flowt[flowno].f[scnt].u = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].v = (float *) malloc(h.np * sizeof(float));
	} else {		/* scalars */
	    tmpe = flowt[flowno].f[scnt].e = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].u = NULL;
	    flowt[flowno].f[scnt].v = NULL;
	}

	if (ElioGetTimeStep(fp, j, &h, &t)) {
	    break;
	}
	if (h.ivs == 2) {
	    ElioMakeVectorsOld(&h, t.d, data);
	} else {
	    ElioMakeScalarsOld(&h, t.d, data);
	}

/* start of timestep */
	flowt[flowno].f[scnt].time = t.t;
/*
        printf("Time = %.1f, scnt = %d\n", t.t, scnt);
*/

	cnt = 0;
	for (i = level; i < (h.np * h.nvrt); i += h.nvrt) {
	    if (h.ivs == 2) {
		flowt[flowno].f[scnt].u[cnt] = data[2 * i];
		flowt[flowno].f[scnt].v[cnt] = data[2 * i + 1];
	    } else {
// hack for selfe elevations
//	if (h.nvrt == 1 && (data[i] + h.d[i]) < 0) {
//printf("fixup %d %f %f\n", h.nvrt, data[i], h.d[i]);
//		tmpe[cnt] = -999.0;
//	} else { 
		tmpe[cnt] = data[i];
//	}
	    }
	    cnt++;
	}

	if (h.ivs == 2) {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].u, mval, &umin, &umax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].u, &umin, &umax);
	    }
	    if (first) {
		tumin = umin;
		tumax = umax;
		first = 0;
	    }
	    if (umin < tumin) {
		tumin = umin;
	    }
	    if (umax > tumax) {
		tumax = umax;
	    }
	    flowt[flowno].f[scnt].umin = umin;
	    flowt[flowno].f[scnt].umax = umax;
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].v, mval, &vmin, &vmax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].v, &vmin, &vmax);
	    }
	    if (scnt == 0) {
		tvmin = vmin;
		tvmax = vmax;
	    }
	    if (vmin < tvmin) {
		tvmin = vmin;
	    }
	    if (vmax > tvmax) {
		tvmax = vmax;
	    }
	    flowt[flowno].f[scnt].vmin = vmin;
	    flowt[flowno].f[scnt].vmax = vmax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = u(%.3lf %.3lf) - v(%.3lf %.3lf)\n", j, umin, umax, vmin, vmax);
	    }
	} else {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].e, mval, &emin, &emax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].e, &emin, &emax);
	    }
	    if (first) {
		temin = emin;
		temax = emax;
		first = 0;
	    }
	    if (emin < temin) {
		temin = emin;
	    }
	    if (emax > temax) {
		temax = emax;
	    }
	    flowt[flowno].f[scnt].emin = emin;
	    flowt[flowno].f[scnt].emax = emax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = %.3lf %.3lf\n", j, emin, emax);
	    }
	}
	scnt++;
    }

  out:;
    flowt[flowno].nsteps = scnt;
/*
    if (nsteps != ndsetse) {
	sprintf(buf, "File is incomplete, scanned %d timesteps of %d\n", nsteps, ndsetse);
	errwin(buf);
    }
*/
    free(readdata);
    ElioFreeHeader(&h);
    ElioFreeTimeStep(&t);

    if (h.ivs == 2) {
	flowt[flowno].umin = tumin;
	flowt[flowno].umax = tumax;
	flowt[flowno].vmin = tvmin;
	flowt[flowno].vmax = tvmax;
    } else {
	flowt[flowno].emin = temin;
	flowt[flowno].emax = temax;
    }
    flowt[flowno].type = h.ivs;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;
    fclose(fp);
    return 1;
}

int ReadElcircSurfNew(int flowno, char *fname, int surf, int start, int stop, int skip, int missing, double mval, int append);
int ReadElcircSurf(int flowno, char *fname, int surf, int start, int stop, int skip, int missing, double mval, int append)
{
    FILE *fp;
    int ftype;
    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "In readbin_elcirc, unable to open file %s\n", fname);
	return 0;
    }
    ftype = ElioGetFileType(fp);
    fclose(fp);
    if (ftype == -1) {
	errwin("ReadElcircSurf(): Unable to read old format files with this version of ACE/vis");
	return 1;
    } else {
	return ReadElcircSurfNew(flowno, fname, surf, start, stop, skip, missing, mval, append);
    }
    return 0;
}

int ReadElcircSurfNew(int flowno, char *fname, int surf, int start, int stop, int skip, int missing, double mval, int append)
{
    FILE *fp;
    int *readdata;
    int cnt, scnt = 0, i, j, k, itmp;
    double emin, emax;
    double temin, temax;
    double umin, umax;
    double tumin, tumax;
    double vmin, vmax;
    double tvmin, tvmax;
    float *tmpe;
    char buf[1024];
    int first, ind, nsteps, err, ftype;
    ElcircHeader h;
    ElcircTimeStep t;
    Grid *g;
    int elevflag = 0;

    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcircSurfNew(): Unable to open file %s\n", fname);
	return 0;
    }
    if (!append) {
	Free_flowt(flowno);
    }

    ftype = ElioGetFileType(fp);

/* Get the header */
    if (err = ElioGetHeader(fname, &h)) {
	fprintf(stderr, "ReadElcircSurfNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
    //ElioPrintHeader(&h);
    if (strstr(buf, "elevation") != NULL) {
	elevflag = 1;
    } else {
	elevflag = 0;
    }

    if (!append) {
	settime_info_start(cg, h.start_time);
	g = &flowt[flowno].g;
	AllocateGrid(g, h.ne, h.np);
	GetElcircGrid(&h, g);
    }

    nsteps = ElioGetNStepsInFile(fname, &h);
/*
    printf("ReadElcircSurfNew(): Number of steps in header = %d, really = %d\n", h.nsteps, scnt);
*/

    if (start < 0) {
	start = 0;
	stop = nsteps - 1;
	skip = 1;
    } else {
	if (start >= nsteps) {
	}
	if (stop >= nsteps) {
	    stop = nsteps - 1;
	}
    }
    scnt = 0;
    readdata = (int *) malloc(nsteps * sizeof(int));
    for (i = 0; i < nsteps; i++) {
	readdata[i] = 0;
    }
    for (i = start; i <= stop; i += skip) {
	readdata[i] = 1;
	scnt++;
    }
    if (scnt == 0) {
	errwin("ReadElcircSurfNew(): No complete time step available, cancelling read");
	flowt[flowno].active = OFF;
	free(readdata);
	return 0;
    }

    if (append) {
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (flowt[flowno].nsteps + scnt) * sizeof(Flow_gen));
	scnt = flowt[flowno].nsteps;
    } else {
	flowt[flowno].npts = h.np;
	flowt[flowno].nsteps = scnt;
	flowt[flowno].f = (Flow_gen *) malloc(scnt * sizeof(Flow_gen));
	scnt = 0;
    }
    if (flowt[flowno].f == NULL) {
	fprintf(stderr, "ReadElcircSurfNew(): Unable to allocate memory for reading timesteps\n");
	free(readdata);
	return 0;
    }

/* pre-allocate memory for timestep struct */
    ElioAllocateTimeStep(&h, &t);
    first = 1;
    for (j = 0; j < nsteps; j++) {
	if (readdata[j] == 0) {
	    continue;
	}
	if (h.ivs == 2) {
	    flowt[flowno].f[scnt].e = NULL;
	    flowt[flowno].f[scnt].u = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].v = (float *) malloc(h.np * sizeof(float));
	} else {		/* scalars */
	    tmpe = flowt[flowno].f[scnt].e = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].u = NULL;
	    flowt[flowno].f[scnt].v = NULL;
	}

	if (ElioGetTimeStep(fp, j, &h, &t)) {
	    break;
	}
/* start of timestep */
	flowt[flowno].f[scnt].time = t.t;
/*
        printf("Time = %.1f, it = %d, ", tt, it);
*/

	cnt = 0;
	for (i = 0; i < h.np; i++) {
	    if (surf == 0) {	/* get the surface */
/* 
   h.bi[] starts from 0, surfind[] starts from one, no[i] points to the start of the
   level data
*/
		ind = (t.surfind[i] == 0) ? -1 : h.no[i] + (t.surfind[i] - h.bi[i] - 1) * h.ivs;
	    } else {		/* get the bottom */
		ind = h.no[i];	/* offset to ith node */
	    }
	    if (h.ivs == 2) {
		flowt[flowno].f[scnt].u[cnt] = (ind == -1) ? 0.0 : t.d[ind];
		flowt[flowno].f[scnt].v[cnt] = (ind == -1) ? 0.0 : t.d[ind + 1];
	    } else {
// hack for elevations
if (elevflag && ind >= 0 && (t.d[ind] + h.d[ind]) < 0.0) {
    tmpe[cnt] = -999.0;
//printf("ind == %d %f\n", ind);
} else {
		tmpe[cnt] = (ind == -1) ? -99.0 : t.d[ind];
}
	    }
	    cnt++;
	}
	if (h.ivs == 2) {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].u, mval, &umin, &umax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].u, &umin, &umax);
	    }
	    if (first) {
		tumin = umin;
		tumax = umax;
		first = 0;
	    }
	    if (umin < tumin) {
		tumin = umin;
	    }
	    if (umax > tumax) {
		tumax = umax;
	    }
	    flowt[flowno].f[scnt].umin = umin;
	    flowt[flowno].f[scnt].umax = umax;
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].v, mval, &vmin, &vmax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].v, &vmin, &vmax);
	    }
	    if (first) {
		tvmin = vmin;
		tvmax = vmax;
		first = 0;
	    }
	    if (vmin < tvmin) {
		tvmin = vmin;
	    }
	    if (vmax > tvmax) {
		tvmax = vmax;
	    }
	    flowt[flowno].f[scnt].vmin = vmin;
	    flowt[flowno].f[scnt].vmax = vmax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = u(%.3lf %.3lf) - v(%.3lf %.3lf)\n", j, umin, umax, vmin, vmax);
	    }
	} else {
/*
int ij;
for (ij=0;ij<h.np;ij += 1000) {
	printf("%d %lf\n", ij, flowt[flowno].f[scnt].e[ij]);
}
*/
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].e, mval, &emin, &emax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].e, &emin, &emax);
	    }
	    if (first) {
		temin = emin;
		temax = emax;
		first = 0;
	    }
	    if (emin < temin) {
		temin = emin;
	    }
	    if (emax > temax) {
		temax = emax;
	    }
	    flowt[flowno].f[scnt].emin = emin;
	    flowt[flowno].f[scnt].emax = emax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = %.3lf %.3lf\n", j, emin, emax);
	    }
	}

	scnt++;
    }

  out:;
    flowt[flowno].nsteps = scnt;
    free(readdata);
    ElioFreeHeader(&h);
    ElioFreeTimeStep(&t);

    if (h.ivs == 2) {
	flowt[flowno].umin = tumin;
	flowt[flowno].umax = tumax;
	flowt[flowno].vmin = tvmin;
	flowt[flowno].vmax = tvmax;
    } else {
	flowt[flowno].emin = temin;
	flowt[flowno].emax = temax;
    }
    flowt[flowno].type = 0;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;
// printf("t1, t2, nsteps: %lf %lf %d\n", flowt[flowno].start, flowt[flowno].stop, flowt[flowno].nsteps);
    fclose(fp);
    return 1;
}

int ReadElcircRegionFile(int flowno, char *fname)
{
    char buf[2048];
    int i, n, npts, *isin, n1, n2, n3;
    double *x, *y;
    FILE *fp;
    Grid *g = &flowt[flowno].g;
    g->ellist = NULL;
    if ((fp = fopen(fname, "rb")) == NULL) {
	sprintf(buf, "Unable to open region file %s", fname);
	errwin(buf);
	return 1;
    }
    if (fgets(buf, 255, fp) == NULL) {
	fprintf(stderr, "Error reading 1st line of region file %s\n", fname);
	return (1);
    }
    if (fgets(buf, 255, fp) == NULL) {
	fprintf(stderr, "Error reading 2nd line of region file %s\n", fname);
	return (1);
    }
    if (fgets(buf, 255, fp) == NULL) {
	fprintf(stderr, "Error reading 3rd line of region file %s\n", fname);
	return (1);
    }
    sscanf(buf, "%d", &npts);
    x = (double *) malloc(npts * sizeof(double));
    y = (double *) malloc(npts * sizeof(double));
    for (i = 0; i < npts; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    fprintf(stderr, "Error reading line %d of region file %s\n", i + 4, fname);
	    free(x);
	    free(y);
	    return (1);
	}
	sscanf(buf, "%lf %lf", &x[i], &y[i]);
    }
    fclose(fp);
    g->ellist = (int *) malloc(g->nmel * sizeof(int));
    isin = (int *) malloc(g->nmnp * sizeof(int));
    n = 1;			/* start from 1 for truth testing (0 is not in) */
    for (i = 0; i < g->nmnp; i++) {
	if (ElioInPolygon(g->xord[i], g->yord[i], npts, x, y)) {
	    isin[i] = n;
	    n++;
	} else {
	    isin[i] = 0;
	}
    }
    for (i = 0; i < g->nmel; i++) {
	n1 = g->icon[i].nl[0];
	n2 = g->icon[i].nl[1];
	n3 = g->icon[i].nl[2];
	if (isin[n1] && isin[n2] && isin[n3]) {
	    g->ellist[i] = 1;
	} else {
	    g->ellist[i] = 0;
	}
    }
    free(isin);
    free(x);
    free(y);
}

int ReplaceElcircGrid(int flowno, char *fname)
{
    char buf[2048];
    int nmel, nmnp;
    FILE *fp;
    Grid *g = &flowt[flowno].g;
    int i, j, itmp, type;
    if ((fp = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "ReplaceElcircGrid(): unable to open file %s\n", fname);
	return 1;
    }
    fgets(buf, 255, fp);	/* first line is ignored */
    fgets(buf, 255, fp);
    convertchar(buf);
    sscanf(buf, "%d %d", &nmel, &nmnp);
/* Number of nodes must be the same, but the table of elements can be different */
    if (g->nmnp != nmnp) {
	fprintf(stderr, "ReplaceElcircGrid(): Grid nodes incorrect, %d %d\n", nmnp, g->nmnp);
	return 1;
    }
    fclose(fp);
    FreeGrid(g);
    AllocateGrid(g, nmel, nmnp);
    ReadGrid(fname, g);
    return 0;
}

/*
 * read fort.63 file
 */
int readbin_adcirc_elev(int flowno, char *gfile, char *fname)
{
    FILE *fp;
    char *s;
    int cnt, scnt = 0, i, j, k, itmp, nmnp;
    double emin, emax;
    double temin, temax;
    float ttime, e, *tmpe, *tmpu, *tmpv;
    char buf[257];
    char rundes[33], runid[25], agrid[25];
    int ndsetse, np, nspoolge, irtype;
    float dtnspoolge;
    int it;
    Flow_gen *tflow;
    float *tfloat;
    Grid *g;

    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "readbin_adcirc_elev(): Unable to open file %s\n", fname);
	return 0;
    }
    for (j = 0; j < flowt[flowno].nsteps; j++) {
	if (flowt[flowno].f) {
	    cxfree(flowt[flowno].f[j].u);
	    cxfree(flowt[flowno].f[j].v);
	    cxfree(flowt[flowno].f[j].e);
	}
    }
    flowt[flowno].nsteps = 0;
    if (flowt[flowno].f) {
	cxfree(flowt[flowno].f);
	flowt[flowno].f = NULL;
    }

    g = &flowt[flowno].g;
    ReadGrid(gfile, g);

    nmnp = g->nmnp;
    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);

    fread(rundes, sizeof(char), 32, fp);
    fread(runid, sizeof(char), 24, fp);
    fread(agrid, sizeof(char), 24, fp);

    read_int(&ndsetse, 1, fp);
    read_int(&np, 1, fp);
    read_float(&dtnspoolge, 1, fp);
    read_int(&nspoolge, 1, fp);
    read_int(&irtype, 1, fp);
    rundes[31] = 0;
    runid[23] = 0;
    agrid[23] = 0;
    if (debuglevel == 13) {
	printf("Reading fort.63 data from file: %s\n", fname);
	printf("%s\n", rundes);
	printf("%s\n", runid);
	printf("%s\n", agrid);
	printf("ndsetse = %d  np = %d\n", ndsetse, np);
	printf("irtype = %d\n", irtype);
	printf("dt*nspoolge = %f  nspoolge = %d\n", dtnspoolge, nspoolge);
    }
    flowt[flowno].nsteps = ndsetse;

    if (nmnp != np) {
	errwin("Number of nodes in data file does not match the number of nodes in the grid");
	return 0;
    }
    flowt[flowno].npts = nmnp;

    flowt[flowno].f = (Flow_gen *) calloc(1, sizeof(Flow_gen));
    scnt = 0;
    for (j = 0; j < flowt[flowno].nsteps; j++) {
	if (j == 0) {
	    tmpe = (float *) calloc(nmnp, sizeof(float));
	}
	tfloat = (float *) calloc(nmnp, sizeof(float));
	if (tfloat == (float *) NULL) {
	    errwin("Insufficient memory for all steps.");
	    break;
	} else {
	    flowt[flowno].f[j].e = tfloat;
	}
	flowt[flowno].f[j].u = NULL;
	flowt[flowno].f[j].v = NULL;
	if (read_float(&ttime, 1, fp) != 1) {
	    /* depart for-loop if out of data in the file */
	    break;
	}
	if (read_int(&it, 1, fp) != 1) {
	    /* depart for-loop if out of data in the file */
	    break;
	}
	if (debuglevel == 13) {
	    printf("time = %f it = %d\n", ttime, it);
	}
	if (read_float(tmpe, np, fp) != np) {
	    /* depart for-loop if out of data in the file */
	    break;
	}
	flowt[flowno].f[j].time = ttime;

	for (i = 0; i < nmnp; i++) {
	    e = flowt[flowno].f[j].e[i] = tmpe[i];
	    if (i == 0) {
		emin = emax = e;
		if (j == 0) {
		    temin = temax = e;
		}
	    }
	    if (e < temin) {
		temin = e;
	    }
	    if (e > temax) {
		temax = e;
	    }
	    if (e < emin) {
		emin = e;
	    }
	    if (e > emax) {
		emax = e;
	    }
	}
	flowt[flowno].f[j].emin = emin;
	flowt[flowno].f[j].emax = emax;
	scnt++;
	tflow = (Flow_gen *) realloc(flowt[flowno].f, (j + 2) * sizeof(Flow_gen));
	if (tflow == (Flow_gen *) NULL) {
	    errwin("Insufficient memory for all steps.");
	    break;
	} else {
	    flowt[flowno].f = tflow;
	}
    }

    if (scnt == 0) {
	errwin("No complete time step read, cancelling read");
	flowt[flowno].active = OFF;
	free(flowt[flowno].f);
	return 0;
    }
    if (scnt != flowt[flowno].nsteps) {
	sprintf(buf, "File claims %d timesteps, read only %d,\n resetting number of steps\n", flowt[flowno].nsteps, scnt);
	errwin(buf);
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, scnt * sizeof(Flow_gen));
	flowt[flowno].nsteps = scnt;
    }
    free(tmpe);
    flowt[flowno].emin = temin;
    flowt[flowno].emax = temax;
    flowt[flowno].type = 0;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;

    fclose(fp);
    return (1);
}

/*
 * flowno: this flow
 * gridno: the domain to use
 * fname: the file name
 * keep: boolean to determine if this is a standalone file or connected
 *       to a fort.63 file - keep == 1 => yes, both fort.63 and fort.64 files
 *       are to be read to this flow.
 */
int readbin_adcirc_flow(int flowno, char *gfile, char *fname, int keep)
{
    FILE *fp;
    char *s;
    int cnt, scnt = 0, i, j, k, itmp, nmnp;
    double umin, umax;
    double tumin, tumax;
    double vmin, vmax;
    double tvmin, tvmax;
    float ttime, u, v, *tmpu, *tmpv;
    char buf[257];
    char rundes[33], runid[25], agrid[25];
    int ndsetse, np, nspoolge, irtype;
    float dtnspoolge;
    int it, err;
    Grid *g;

    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "In readflow, unable to open file %s\n", fname);
	return 0;
    }
    if (!keep) {
	for (j = 0; j < flowt[flowno].nsteps; j++) {
	    if (flowt[flowno].f) {
		cxfree(flowt[flowno].f[j].u);
		cxfree(flowt[flowno].f[j].v);
		cxfree(flowt[flowno].f[j].e);
	    }
	}
	flowt[flowno].nsteps = 0;
	if (flowt[flowno].f) {
	    cxfree(flowt[flowno].f);
	    flowt[flowno].f = NULL;
	}
    }

    g = &flowt[flowno].g;
    ReadGrid(gfile, g);

    nmnp = g->nmnp;
    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);

    fread(rundes, sizeof(char), 32, fp);
    fread(runid, sizeof(char), 24, fp);
    fread(agrid, sizeof(char), 24, fp);

    read_int(&ndsetse, 1, fp);
    read_int(&np, 1, fp);
    read_float(&dtnspoolge, 1, fp);
    read_int(&nspoolge, 1, fp);
    read_int(&irtype, 1, fp);
    rundes[31] = 0;
    runid[23] = 0;
    agrid[23] = 0;
    if (debuglevel == 13) {
	printf("%s\n", rundes);
	printf("%s\n", runid);
	printf("%s\n", agrid);
	printf("ndsetse = %d  np = %d\n", ndsetse, np);
	printf("dt*nspoolge = %f  nspoolge = %d\n", dtnspoolge, nspoolge);
	printf("irtype = %d\n", irtype);
    }
    flowt[flowno].nsteps = ndsetse;

    if (nmnp != np) {
	errwin("Number of nodes in data file does not match the number of nodes in the grid");
	return 0;
    }
    flowt[flowno].npts = nmnp;

    flowt[flowno].f = (Flow_gen *) malloc(sizeof(Flow_gen));
    scnt = 0;
/*
 * note that U and V are interlaced, so both U and V are read in at one time
 * then placed into the proper arrays
 */
    for (j = 0; j < flowt[flowno].nsteps; j++) {
	if (j == 0) {
	    tmpu = (float *) malloc(2 * nmnp * sizeof(float));	/* twice the size, for U
								 * and V */
	}
	if (!keep) {
	    flowt[flowno].f[j].e = NULL;
	}
	flowt[flowno].f[j].u = (float *) malloc(nmnp * sizeof(float));
	flowt[flowno].f[j].v = (float *) malloc(nmnp * sizeof(float));
	if (read_float(&ttime, 1, fp) != 1) {
	    /* depart for-loop if out of data in the file */
	    if (debuglevel == 13) {
		printf("Ended at reading time\n");
	    }
	    break;
	}
	if (read_int(&it, 1, fp) != 1) {
	    /* depart for-loop if out of data in the file */
	    if (debuglevel == 13) {
		printf("Ended at reading step number\n");
	    }
	    break;
	}
	if (debuglevel == 13) {
	    printf("time = %f it = %d\n", ttime, it);
	}
	if ((err = read_float(tmpu, 2 * np, fp)) != (2 * np)) {
	    /* depart for-loop if out of data in the file */
	    if (debuglevel == 13) {
		printf("Ended at reading uv, %d %d\n", err, 2 * np);
	    }
	    break;
	}
	flowt[flowno].f[j].time = ttime;
/*
 * straighten out U and V
 */
	for (i = 0; i < nmnp; i++) {
	    flowt[flowno].f[j].u[i] = tmpu[2 * i];
	    flowt[flowno].f[j].v[i] = tmpu[2 * i + 1];
	}

	for (i = 0; i < nmnp; i++) {
	    u = flowt[flowno].f[j].u[i];
	    v = flowt[flowno].f[j].v[i];
	    if (i == 0) {
		umin = umax = u;
		vmin = vmax = v;
		if (j == 0) {
		    tumin = tumax = u;
		    tvmin = tvmax = v;
		}
	    }
	    if (u < tumin) {
		tumin = u;
	    }
	    if (u > tumax) {
		tumax = u;
	    }
	    if (u < umin) {
		umin = u;
	    }
	    if (u > umax) {
		umax = u;
	    }
	    if (v < tvmin) {
		tvmin = v;
	    }
	    if (v > tvmax) {
		tvmax = v;
	    }
	    if (v < vmin) {
		vmin = v;
	    }
	    if (v > vmax) {
		vmax = v;
	    }
	}
	flowt[flowno].f[j].umin = umin;
	flowt[flowno].f[j].umax = umax;
	flowt[flowno].f[j].vmin = vmin;
	flowt[flowno].f[j].vmax = vmax;
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (j + 2) * sizeof(Flow_gen));
	scnt++;
    }

    if (scnt == 0) {
	errwin("No complete time step read, cancelling read");
	flowt[flowno].active = OFF;
	free(flowt[flowno].f);
	return 0;
    }
    if (scnt != flowt[flowno].nsteps) {
	sprintf(buf, "File claims %d timesteps, read only %d,\n resetting number of steps\n", flowt[flowno].nsteps, scnt);
	errwin(buf);
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, scnt * sizeof(Flow_gen));
	flowt[flowno].nsteps = scnt;
    }
    free(tmpu);
    flowt[flowno].umin = tumin;
    flowt[flowno].umax = tumax;
    flowt[flowno].vmin = tvmin;
    flowt[flowno].vmax = tvmax;
    flowt[flowno].type = 0;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;

    fclose(fp);
    return (1);
}

/*
 * Read values at a given depth, interpolating for each node, for sigma grids.
 *
 * fname is data files, zname is zcor file.
 *
 */
int ReadElcircDepth(int flowno, char *fname, char *zname, double depth, int start, int stop, int skip, int missing, double mval, int append)
{
    FILE *fp, *fpz;
    int *readdata;
    int ind, cnt, scnt = 0, i, j, k, itmp;
    float *data, *dataz;
    double emin, emax;
    double temin, temax;
    double umin, umax;
    double tumin, tumax;
    double vmin, vmax;
    double tvmin, tvmax;
    double u[200], v[200], s[200], z[200];
    char buf[1024];
    int first, nsteps, err, ftype;
    long offset;
    ElcircHeader h;
    ElcircTimeStep t;
    ElcircHeader hz;
    ElcircTimeStep tz;
    Grid *g;
    char *ptr = (char *) NULL;

    if (zname == (char *) NULL) {
	int nc = strlen(fname);
	zname = (char *) malloc(strlen(fname) + 1);
	strcpy(zname, fname);
	if (nc > 4) {
	    ptr = strstr(zname, ".6");
	    if (ptr != NULL) {
		ptr -= 4;
		strcpy(ptr, "zcor.63");
	    }
	}
    }
//printf("%d %s %s %lf %d %d %d\n", flowno, fname, zname, depth, start, stop, skip);
    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to open file %s\n", fname);
	return 0;
    }
    if ((fpz = fopen(zname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to open file %s\n", zname);
	return 0;
    }
    if (!append) {
	Free_flowt(flowno);
    }

    ftype = ElioGetFileType(fp);
/* Get the header */
    if (err = ElioGetHeader(fname, &h)) {
	fprintf(stderr, "ReadElcircNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
    if (err = ElioGetHeader(zname, &hz)) {
	fprintf(stderr, "ReadElcircNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
/*
    ElioPrintHeader(h);
*/

/*
   set the time.
   12/02/2002 00:00:00 PST
   YYYYMMDDHHMMSS
*/
    if (!append) {
	settime_info_start(cg, h.start_time);
	g = &flowt[flowno].g;
	GetElcircGrid(&h, g);
    }

    nsteps = ElioGetNStepsInFile(fname, &h);
/*
    printf("ReadElcircNew(): Number of steps in header = %d, really = %d\n", h.nsteps, scnt);
*/
    if (start < 0) {
	start = 0;
	stop = nsteps - 1;
	skip = 1;
    } else {
	if (start >= nsteps) {
	}
	if (stop >= nsteps) {
	    stop = nsteps - 1;
	}
    }
    scnt = 0;
    readdata = (int *) malloc(nsteps * sizeof(int));
    for (i = 0; i < nsteps; i++) {
	readdata[i] = 0;
    }
    for (i = start; i <= stop; i += skip) {
	readdata[i] = 1;
	scnt++;
    }
    if (scnt == 0) {
	errwin("ReadElcircNew(): No complete time step available, cancelling read");
	flowt[flowno].active = OFF;
	free(readdata);
	return 0;
    }

    if (append) {
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (flowt[flowno].nsteps + scnt) * sizeof(Flow_gen));
	scnt = flowt[flowno].nsteps;
    } else {
	flowt[flowno].npts = h.np;
	flowt[flowno].nsteps = scnt;
	flowt[flowno].f = (Flow_gen *) malloc(scnt * sizeof(Flow_gen));
	scnt = 0;
    }

    if (debuglevel == 13) {
	printf("Trying to read %d steps:\n", flowt[flowno].nsteps);
    }

    if (flowt[flowno].f == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to allocate memory for reading timesteps\n");
	free(readdata);
	return 0;
    }

/* pre-allocate memory for timestep struct */
    ElioAllocateTimeStep(&h, &t);
/* zcor */
    ElioAllocateTimeStep(&hz, &tz);

    first = 1;
    for (j = 0; j < nsteps; j++) {
	if (readdata[j] == 0) {
	    continue;
	}
	if (h.ivs == 2) {	/* Vector data */
	    flowt[flowno].f[scnt].e = NULL;
	    flowt[flowno].f[scnt].u = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].v = (float *) malloc(h.np * sizeof(float));
	} else {		/* scalars */
	    flowt[flowno].f[scnt].e = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].u = NULL;
	    flowt[flowno].f[scnt].v = NULL;
	}

	if (ElioGetTimeStep(fp, j, &h, &t)) {
	    break;
	}
	if (ElioGetTimeStep(fpz, j, &hz, &tz)) {
	    break;
	}
	if (h.v == 5) {
	    data = (float *) malloc(h.nvrt * h.np * h.ivs * 4);
	    dataz = (float *) malloc(hz.nvrt * hz.np * hz.ivs * 4);
	    if (h.ivs == 2) {
		ElioMakeVectorsOld(&h, t.d, data);
	    } else {
		ElioMakeScalarsOld(&h, t.d, data);
	    }
	    t.d = data;
	    ElioMakeDepthsOld(&hz, tz.d, dataz);
	    tz.d = dataz;
	}
/* start of timestep */
	flowt[flowno].f[scnt].time = t.t;
	cnt = 0;
	for (i = 0; i < h.np; i++) {
/* do interpolation here */
	    if (h.ivs == 2) {
		for (k = 0; k < h.nvrt; k++) {
		    z[k] = tz.d[i * h.nvrt + k];
		    u[k] = t.d[2 * i * h.nvrt + 2 * k];
		    v[k] = t.d[2 * i * h.nvrt + 2 * k + 1];
		}
		ind = ElioFindIndex(hz.nvrt, z, -depth);
		if (ind == hz.nvrt || ind == -1) {
		    flowt[flowno].f[scnt].u[cnt] = 0.0;
		    flowt[flowno].f[scnt].v[cnt] = 0.0;
		} else {
		    flowt[flowno].f[scnt].u[cnt] = ElioInterpolateAtIndex(h.nvrt, z, u, ind, -depth);
		    flowt[flowno].f[scnt].v[cnt] = ElioInterpolateAtIndex(h.nvrt, z, v, ind, -depth);
		}
	    } else {
		for (k = 0; k < h.nvrt; k++) {
		    z[k] = tz.d[i * h.nvrt + k];
		    s[k] = t.d[i * h.nvrt + k];
		}
		ind = ElioFindIndex(hz.nvrt, z, -depth);
/* out of range, skip */
		if (ind >= hz.nvrt || ind <= -1) {
		    flowt[flowno].f[scnt].e[cnt] = -9999.0;
		} else {
//printf("Index = %d\n", ind);
		    flowt[flowno].f[scnt].e[cnt] = ElioInterpolateAtIndex(h.nvrt, z, s, ind, -depth);
		}
	    }
	    cnt++;
	}

	if (h.ivs == 2) {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].u, mval, &umin, &umax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].u, &umin, &umax);
	    }
	    if (first) {
		tumin = umin;
		tumax = umax;
		first = 0;
	    }
	    if (umin < tumin) {
		tumin = umin;
	    }
	    if (umax > tumax) {
		tumax = umax;
	    }
	    flowt[flowno].f[scnt].umin = umin;
	    flowt[flowno].f[scnt].umax = umax;
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].v, mval, &vmin, &vmax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].v, &vmin, &vmax);
	    }
	    if (scnt == 0) {
		tvmin = vmin;
		tvmax = vmax;
	    }
	    if (vmin < tvmin) {
		tvmin = vmin;
	    }
	    if (vmax > tvmax) {
		tvmax = vmax;
	    }
	    flowt[flowno].f[scnt].vmin = vmin;
	    flowt[flowno].f[scnt].vmax = vmax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = u(%.3lf %.3lf) - v(%.3lf %.3lf)\n", j, umin, umax, vmin, vmax);
	    }
	} else {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].e, mval, &emin, &emax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].e, &emin, &emax);
	    }
	    if (first) {
		temin = emin;
		temax = emax;
		first = 0;
	    }
	    if (emin < temin) {
		temin = emin;
	    }
	    if (emax > temax) {
		temax = emax;
	    }
	    flowt[flowno].f[scnt].emin = emin;
	    flowt[flowno].f[scnt].emax = emax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = %.3lf %.3lf\n", j, emin, emax);
	    }
	}
	scnt++;
    }

  out:;
    flowt[flowno].nsteps = scnt;
/*
    if (nsteps != ndsetse) {
	sprintf(buf, "File is incomplete, scanned %d timesteps of %d\n", nsteps, ndsetse);
	errwin(buf);
    }
*/
    free(readdata);
    ElioFreeHeader(&h);
    ElioFreeTimeStep(&t);
    ElioFreeHeader(&hz);
    ElioFreeTimeStep(&tz);

    if (h.ivs == 2) {
	flowt[flowno].umin = tumin;
	flowt[flowno].umax = tumax;
	flowt[flowno].vmin = tvmin;
	flowt[flowno].vmax = tvmax;
    } else {
	flowt[flowno].emin = temin;
	flowt[flowno].emax = temax;
    }
    flowt[flowno].type = h.ivs;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;
    fclose(fp);
    fclose(fpz);
    return 1;
}

/*
 * Read values at a given depth below the free surface, interpolating for each node, for sigma grids.
 *
 * fname is data files, zname is zcor file.
 *
 */
int ReadElcircDepthFromFreeSurface(int flowno, char *fname, char *zname, double depth, int start, int stop, int skip, int missing, double mval, int append)
{
    FILE *fp, *fpz;
    int *readdata;
    int ind, cnt, scnt = 0, i, j, k, itmp;
    float *data, *dataz;
    double emin, emax;
    double temin, temax;
    double umin, umax;
    double tumin, tumax;
    double vmin, vmax;
    double tvmin, tvmax;
    double u[200], v[200], s[200], z[200];
    char buf[1024];
    int first, nsteps, err, ftype;
    long offset;
    ElcircHeader h;
    ElcircTimeStep t;
    ElcircHeader hz;
    ElcircTimeStep tz;
    Grid *g;
    char *ptr = (char *) NULL;

    if (zname == (char *) NULL) {
	int nc = strlen(fname);
	zname = (char *) malloc(strlen(fname) + 1);
	strcpy(zname, fname);
	if (nc > 4) {
	    ptr = strstr(zname, ".6");
	    if (ptr != NULL) {
		ptr -= 4;
		strcpy(ptr, "zcor.63");
	    }
	}
    }
//printf("%d %s %s %lf %d %d %d\n", flowno, fname, zname, depth, start, stop, skip);
    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to open file %s\n", fname);
	return 0;
    }
    if ((fpz = fopen(zname, "rb")) == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to open file %s\n", zname);
	return 0;
    }
    if (!append) {
	Free_flowt(flowno);
    }

    ftype = ElioGetFileType(fp);
/* Get the header */
    if (err = ElioGetHeader(fname, &h)) {
	fprintf(stderr, "ReadElcircNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
    if (err = ElioGetHeader(zname, &hz)) {
	fprintf(stderr, "ReadElcircNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
/*
    ElioPrintHeader(h);
*/

/*
   set the time.
   12/02/2002 00:00:00 PST
   YYYYMMDDHHMMSS
*/
    if (!append) {
	settime_info_start(cg, h.start_time);
	g = &flowt[flowno].g;
	GetElcircGrid(&h, g);
    }

    nsteps = ElioGetNStepsInFile(fname, &h);
/*
    printf("ReadElcircNew(): Number of steps in header = %d, really = %d\n", h.nsteps, scnt);
*/
    if (start < 0) {
	start = 0;
	stop = nsteps - 1;
	skip = 1;
    } else {
	if (start >= nsteps) {
	}
	if (stop >= nsteps) {
	    stop = nsteps - 1;
	}
    }
    scnt = 0;
    readdata = (int *) malloc(nsteps * sizeof(int));
    for (i = 0; i < nsteps; i++) {
	readdata[i] = 0;
    }
    for (i = start; i <= stop; i += skip) {
	readdata[i] = 1;
	scnt++;
    }
    if (scnt == 0) {
	errwin("ReadElcircNew(): No complete time step available, cancelling read");
	flowt[flowno].active = OFF;
	free(readdata);
	return 0;
    }

    if (append) {
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (flowt[flowno].nsteps + scnt) * sizeof(Flow_gen));
	scnt = flowt[flowno].nsteps;
    } else {
	flowt[flowno].npts = h.np;
	flowt[flowno].nsteps = scnt;
	flowt[flowno].f = (Flow_gen *) malloc(scnt * sizeof(Flow_gen));
	scnt = 0;
    }

    if (debuglevel == 13) {
	printf("Trying to read %d steps:\n", flowt[flowno].nsteps);
    }

    if (flowt[flowno].f == NULL) {
	fprintf(stderr, "ReadElcircNew(): Unable to allocate memory for reading timesteps\n");
	free(readdata);
	return 0;
    }

/* pre-allocate memory for timestep struct */
    ElioAllocateTimeStep(&h, &t);
/* zcor */
    ElioAllocateTimeStep(&hz, &tz);

    first = 1;
    for (j = 0; j < nsteps; j++) {
	if (readdata[j] == 0) {
	    continue;
	}
	if (h.ivs == 2) {	/* Vector data */
	    flowt[flowno].f[scnt].e = NULL;
	    flowt[flowno].f[scnt].u = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].v = (float *) malloc(h.np * sizeof(float));
	} else {		/* scalars */
	    flowt[flowno].f[scnt].e = (float *) malloc(h.np * sizeof(float));
	    flowt[flowno].f[scnt].u = NULL;
	    flowt[flowno].f[scnt].v = NULL;
	}

	if (ElioGetTimeStep(fp, j, &h, &t)) {
	    break;
	}
	if (ElioGetTimeStep(fpz, j, &hz, &tz)) {
	    break;
	}
	if (h.v == 5) {
	    data = (float *) malloc(h.nvrt * h.np * h.ivs * 4);
	    dataz = (float *) malloc(hz.nvrt * hz.np * hz.ivs * 4);
	    if (h.ivs == 2) {
		ElioMakeVectorsOld(&h, t.d, data);
	    } else {
		ElioMakeScalarsOld(&h, t.d, data);
	    }
	    t.d = data;
	    ElioMakeDepthsOld(&hz, tz.d, dataz);
	    tz.d = dataz;
	}
/* start of timestep */
	flowt[flowno].f[scnt].time = t.t;
	cnt = 0;
	for (i = 0; i < h.np; i++) {
/* do interpolation here */
	    if (h.ivs == 2) {
		for (k = 0; k < h.nvrt; k++) {
		    z[k] = tz.d[i * h.nvrt + k];
		    u[k] = t.d[2 * i * h.nvrt + 2 * k];
		    v[k] = t.d[2 * i * h.nvrt + 2 * k + 1];
		}
		ind = ElioFindIndex(hz.nvrt, z, -depth + t.e[i]);
		if (ind == hz.nvrt || ind == -1) {
		    flowt[flowno].f[scnt].u[cnt] = 0.0;
		    flowt[flowno].f[scnt].v[cnt] = 0.0;
		} else {
		    flowt[flowno].f[scnt].u[cnt] = ElioInterpolateAtIndex(h.nvrt, z, u, ind, -depth + t.e[i]);
		    flowt[flowno].f[scnt].v[cnt] = ElioInterpolateAtIndex(h.nvrt, z, v, ind, -depth + t.e[i]);
		}
	    } else {
		double ss = 0.0;
		for (k = 0; k < h.nvrt; k++) {
		    z[k] = tz.d[i * h.nvrt + k];
		    s[k] = t.d[i * h.nvrt + k];
		}
		ind = ElioFindIndex(hz.nvrt, z, -depth + t.e[i]);
/* out of range, skip */
		if (ind >= hz.nvrt || ind <= -1) {
		    flowt[flowno].f[scnt].e[cnt] = -9999.0;
		} else {
		    ss = flowt[flowno].f[scnt].e[cnt] = ElioInterpolateAtIndex(h.nvrt, z, s, ind, -depth + t.e[i]);
//printf("Index = %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", ind, flowno, -depth, t.e[i], z[ind], -depth + t.e[i], z[ind+1], s[ind], ss, s[ind+1]);
		}
	    }
	    cnt++;
	}

	if (h.ivs == 2) {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].u, mval, &umin, &umax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].u, &umin, &umax);
	    }
	    if (first) {
		tumin = umin;
		tumax = umax;
		first = 0;
	    }
	    if (umin < tumin) {
		tumin = umin;
	    }
	    if (umax > tumax) {
		tumax = umax;
	    }
	    flowt[flowno].f[scnt].umin = umin;
	    flowt[flowno].f[scnt].umax = umax;
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].v, mval, &vmin, &vmax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].v, &vmin, &vmax);
	    }
	    if (scnt == 0) {
		tvmin = vmin;
		tvmax = vmax;
	    }
	    if (vmin < tvmin) {
		tvmin = vmin;
	    }
	    if (vmax > tvmax) {
		tvmax = vmax;
	    }
	    flowt[flowno].f[scnt].vmin = vmin;
	    flowt[flowno].f[scnt].vmax = vmax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = u(%.3lf %.3lf) - v(%.3lf %.3lf)\n", j, umin, umax, vmin, vmax);
	    }
	} else {
	    if (missing) {
		minmaxm(h.np, flowt[flowno].f[scnt].e, mval, &emin, &emax);
	    } else {
		minmax(h.np, flowt[flowno].f[scnt].e, &emin, &emax);
	    }
	    if (first) {
		temin = emin;
		temax = emax;
		first = 0;
	    }
	    if (emin < temin) {
		temin = emin;
	    }
	    if (emax > temax) {
		temax = emax;
	    }
	    flowt[flowno].f[scnt].emin = emin;
	    flowt[flowno].f[scnt].emax = emax;
	    if (debuglevel == 13) {
		printf("index = %d, min/max = %.3lf %.3lf\n", j, emin, emax);
	    }
	}
	scnt++;
    }

  out:;
    flowt[flowno].nsteps = scnt;
/*
    if (nsteps != ndsetse) {
	sprintf(buf, "File is incomplete, scanned %d timesteps of %d\n", nsteps, ndsetse);
	errwin(buf);
    }
*/
    free(readdata);
    ElioFreeHeader(&h);
    ElioFreeTimeStep(&t);
    ElioFreeHeader(&hz);
    ElioFreeTimeStep(&tz);

    if (h.ivs == 2) {
	flowt[flowno].umin = tumin;
	flowt[flowno].umax = tumax;
	flowt[flowno].vmin = tvmin;
	flowt[flowno].vmax = tvmax;
    } else {
	flowt[flowno].emin = temin;
	flowt[flowno].emax = temax;
    }
    flowt[flowno].type = h.ivs;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;
    fclose(fp);
    fclose(fpz);
    return 1;
}
