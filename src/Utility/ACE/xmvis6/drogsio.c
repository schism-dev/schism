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
 *  drogsio.c - read drogues
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drogsio.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

int readdrogues(int drno, char *fname, int start, int stop, int skip)
{
    char alphid[257], buf[257];
    FILE *fp;
    int i, k, npts, nsteps, ntotal, scnt = 0, *readdata;
    double xtmp, ytmp, time;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "Error in procedure readdrogues(): can't open %s\n", fname);
	errwin(buf);
	return 0;
    }

    Free_drogues(drno);

    if (fgets(alphid, 256, fp) != NULL) {
	alphid[strlen(alphid) - 1] = 0;
	strcpy(drogues[drno].fname, alphid);
    } else {
	sprintf(buf, "Error in readdrogues(): can't read alphid from %s\n", fname);
	errwin(buf);
	fclose(fp);
	return 0;
    }

    if (fgets(alphid, 256, fp) != NULL) {
	sscanf(alphid, "%d", &nsteps);
    } else {
	sprintf(buf, "Error in readdrogues(): can't read #steps from %s\n", fname);
	errwin(buf);
	fclose(fp);
	return 0;
    }

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
	errwin("No complete time step available, cancelling read");
	free(readdata);
	fclose(fp);
	return 0;
    }

    Allocate_drogues(drno, scnt);

    scnt = 0;

    for (k = 0; k < nsteps; k++) {
	if (fgets(buf, 256, fp) != NULL) {
	    /*convertchar(buf); */
	    sscanf(buf, "%lf %d", &time, &npts);
	    if (readdata[k]) {
		Allocate_drogstep(drno, scnt, time, npts, 0);
	    }
	    for (i = 0; i < npts; i++) {
		if (fgets(buf, 256, fp) != NULL) {
		    convertchar(buf);
		    if (readdata[k]) {
			sscanf(buf, "%*s %lf %lf %lf", &drogues[drno].p[scnt].x[i], &drogues[drno].p[scnt].y[i], &drogues[drno].p[scnt].z[i]);
			drogues[drno].p[scnt].drnum[i] = i + 1;
		    }
		} else {
		    goto out;
		}
	    }
	    if (readdata[k]) {
		scnt++;
	    }
	} else {
	    ntotal = scnt;
	    break;
	}
    }
  out:;
    ntotal = scnt;
    drogues[drno].nsteps = ntotal;
    drogues[drno].start = drogues[drno].p[0].time;
    drogues[drno].stop = drogues[drno].p[drogues[drno].nsteps - 1].time;
    fclose(fp);
    return 1;
}
