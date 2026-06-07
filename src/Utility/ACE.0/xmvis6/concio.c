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
 *  concio.c - read/write conc fields
 *
 */

#ifndef lint
static char RCSid[] = "$Id: concio.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

int readelaconc(int concno, int gridno, char *fname)
{
    FILE *fp;
    char *s;
    int i, j, k, itmp, n, nmnp = grid[gridno].nmnp;
    double c, cmin, cmax;
    double tcmin, tcmax;
    char buf[257];

    if ((fp = fopen(fname, "r")) == NULL) {
	return 0;
    }
    if (elaconc[concno].data != NULL) {
	for (j = 0; j < elaconc[concno].nsteps; j++) {
	    free(elaconc[concno].data[j].c);
	}
	free(elaconc[concno].data);
	elaconc[concno].data = NULL;
    }
    n = 0;
    while (fgets(buf, 256, fp) != NULL) {
	if (elaconc[concno].data == NULL) {
	    elaconc[concno].data = (Conc_gen2d *) calloc(n + 1, sizeof(Conc_gen2d));
	} else {
	    elaconc[concno].data = (Conc_gen2d *) realloc(elaconc[concno].data, (n + 1) * sizeof(Conc_gen2d));
	}
	sscanf(buf, "%lf", &elaconc[concno].data[n].time);
	elaconc[concno].data[n].c = (double *) calloc(nmnp, sizeof(double));

	for (i = 0; i < nmnp; i++) {
	    if (fgets(buf, 256, fp) == NULL) {
		goto bustout;
	    }
	    sscanf(buf, "%*d %lf", &elaconc[concno].data[n].c[i]);
	    c = elaconc[concno].data[n].c[i];
	    if (i == 0) {
		cmin = cmax = c;
		if (j == 0) {
		    tcmin = tcmax = c;
		}
	    }
	    if (c < cmin) {
		cmin = c;
	    }
	    if (c > cmax) {
		cmax = c;
	    }
	    if (c < tcmin) {
		tcmin = c;
	    }
	    if (c > tcmax) {
		tcmax = c;
	    }
	}
	n++;
    }
  bustout:;
    elaconc[concno].active = ON;
    elaconc[concno].grid = gridno;
    elaconc[concno].cmin = cmin;
    elaconc[concno].cmax = cmax;
    elaconc[concno].nsteps = n;
    elaconc[concno].npts = nmnp;
    elaconc[concno].start = elaconc[concno].data[0].time;
    elaconc[concno].stop = elaconc[concno].data[elaconc[concno].nsteps - 1].time;
    if (elaconc[concno].nsteps > 0) {
	elaconc[concno].step = elaconc[concno].data[1].time - elaconc[concno].data[0].time;
    } else {
	elaconc[concno].step = 0.0;
    }
    fclose(fp);
    return (1);
}

int readelaconcbin(int concno, int gridno, char *fname)
{
    FILE *fp;
    char *s;
    int ftype, i, j, k, n, itmp, nmnp = grid[gridno].nmnp;
    double c, cmin, cmax;
    double tcmin, tcmax;
    float *cr, ttmp;
    char buf[257];

    if ((fp = fopen(fname, "r")) == NULL) {
	return 0;
    }
    if (elaconc[concno].data != NULL) {
	for (j = 0; j < elaconc[concno].nsteps; j++) {
	    free(elaconc[concno].data[j].c);
	}
	free(elaconc[concno].data);
	elaconc[concno].data = NULL;
    }
    fread(&itmp, sizeof(int), 1, fp);
    fread(&ftype, sizeof(int), 1, fp);
    fread(&itmp, sizeof(int), 1, fp);

    if (!(ftype == 22 || ftype == 21)) {
	fclose(fp);
	return 0;		/* not a proper binary file */
    }
    fread(&itmp, sizeof(int), 1, fp);
    fread(&nmnp, sizeof(int), 1, fp);
    fread(&itmp, sizeof(int), 1, fp);

    fread(&itmp, sizeof(int), 1, fp);
    fread(&n, sizeof(int), 1, fp);
    fread(&itmp, sizeof(int), 1, fp);
    elaconc[concno].nsteps = n;

    if (ftype == 22) {
	cr = (float *) calloc(nmnp, sizeof(float));
    }
    elaconc[concno].data = (Conc_gen2d *) calloc(n, sizeof(Conc_gen2d));
    for (j = 0; j < n; j++) {
	fread(&itmp, sizeof(int), 1, fp);
	fread(&ttmp, sizeof(float), 1, fp);
	fread(&itmp, sizeof(int), 1, fp);
	elaconc[concno].data[j].time = ttmp;

	elaconc[concno].data[j].c = (double *) calloc(nmnp, sizeof(double));

	fread(&itmp, sizeof(int), 1, fp);
	if (ftype == 21) {
	    fread(elaconc[concno].data[j].c, nmnp, sizeof(double), fp);
	} else if (ftype == 22) {
	    fread(cr, sizeof(float), nmnp, fp);
	    for (k = 0; k < nmnp; k++) {
		c = elaconc[concno].data[j].c[k] = cr[k];
		if (k == 0) {
		    cmin = cmax = c;
		    if (j == 0) {
			tcmin = tcmax = c;
		    }
		}
		if (c < cmin) {
		    cmin = c;
		}
		if (c > cmax) {
		    cmax = c;
		}
		if (c < tcmin) {
		    tcmin = c;
		}
		if (c > tcmax) {
		    tcmax = c;
		}
	    }
	}
	fread(&itmp, sizeof(int), 1, fp);
	elaconc[concno].data[j].cmin = cmin;
	elaconc[concno].data[j].cmax = cmax;
    }
    free(cr);
    elaconc[concno].active = ON;
    elaconc[concno].grid = gridno;
    elaconc[concno].cmin = tcmin;
    elaconc[concno].cmax = tcmax;
    elaconc[concno].nsteps = n;
    elaconc[concno].npts = nmnp;
    elaconc[concno].start = elaconc[concno].data[0].time;
    elaconc[concno].stop = elaconc[concno].data[elaconc[concno].nsteps - 1].time;
    if (elaconc[concno].nsteps > 0) {
	elaconc[concno].step = elaconc[concno].data[1].time - elaconc[concno].data[0].time;
    } else {
	elaconc[concno].step = 0.0;
    }
    fclose(fp);
    return (1);
}
