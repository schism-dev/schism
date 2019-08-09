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
 * utilities for tide stations
 */

#include "symdefs.h"
#include "defines.h"
#include "globals.h"

/*
 * a few tidal constituents
 */
TidalConstituent tc[] = {
    "S0", 0.0000000000E+00, 1e307,
    "Mm", 0.2639202648E-05, 661.309,
    "MSf", 0.4925201210E-05, 354.367,
    "Mf", 0.5323414371E-05, 327.859,
    "O1", 0.6759774260E-04, 25.819,
    "P1", 0.7252294745E-04, 24.066,
    "K1", 0.7292116061E-04, 23.934,
    "N2", 0.1378797024E-03, 12.658,
    "M2", 0.1405188959E-03, 12.421,
    "L2", 0.1431581040E-03, 12.192,
    "S2", 0.1454441081E-03, 12.000,
    "K2", 0.1458423212E-03, 11.967,
    "2MS2", 0.1355936984E-03, 12.872,
    "MN4", 0.2783985983E-03, 6.269,
    "MS4", 0.2859630040E-03, 6.103,
    "2MN6", 0.4189175088E-03, 4.166,
    "2MS6", 0.4264819145E-03, 4.092,
    "2SM6", 0.4314070975E-03, 4.046,
    "M4", 0.2810377919E-03, 6.210,
    "S4", 0.2908882161E-03, 6.000,
    "M6", 0.4215567023E-03, 4.140,
    "S6", 0.4363323096E-03, 4.000,
    "M8", 0.5620755837E-03, 3.105,
    "Sa", 0.1990969736E-06, 8766.227,
    "Ssa", 0.3982128760E-06, 4382.905,
    "2Q1", 0.6231934094E-04, 28.006,
    "SIGMA1", 0.6267254503E-04, 27.848,
    "Q1", 0.6495854177E-04, 26.868,
    "RO1", 0.6531174586E-04, 26.723,
    "M1", 0.7028195250E-04, 24.833,
    "CSI1", 0.7063515659E-04, 24.709,
    "PI1", 0.7232384814E-04, 24.132,
    "S1", 0.7272205403E-04, 24.000,
    "PSI1", 0.7312025264E-04, 23.869,
    "FI1", 0.7331937377E-04, 23.804,
    "TETA1", 0.7520715735E-04, 23.207,
    "J1", 0.7556036144E-04, 23.098,
    "2N2", 0.1352404943E-03, 12.905,
    "NIU2", 0.1382329065E-03, 12.626,
    "(HORN1)", 0.1403197966E-03, 12.438,
    "(HORN2)", 0.1407180098E-03, 12.403,
    "LAMBDA2", 0.1428049000E-03, 12.222,
    "T2", 0.1452450088E-03, 12.016,
    "R2", 0.1456432074E-03, 11.984,
    "ZETA2", 0.1481283107E-03, 11.783,
    "MP1", 0.6799595576E-04, 25.668,
    "SO1", 0.7784635818E-04, 22.420,
    "OO1", 0.7824457134E-04, 22.306,
    "OQ2", 0.1325562916E-03, 13.167,
    "MNS2", 0.1329545048E-03, 13.127,
    "OP2", 0.1401206828E-03, 12.456,
    "MKS2", 0.1409171091E-03, 12.386,
    "MSN2", 0.1480833016E-03, 11.786,
    "KJ2", 0.1484815148E-03, 11.755,
    "2SM2", 0.1503693056E-03, 11.607,
    "MO3", 0.2081166458E-03, 8.386,
    "SO3", 0.2130418434E-03, 8.192,
    "MK3", 0.2134400565E-03, 8.177,
    "SK3", 0.2183652687E-03, 7.993,
    "SN4", 0.2833238104E-03, 6.160,
    "MK4", 0.2863612317E-03, 6.095,
    "SK4", 0.2912864147E-03, 5.992,
    "MSN6", 0.4238427209E-03, 4.118,
    "2MK6", 0.4268801131E-03, 4.089,
    "MSK6", 0.4318053252E-03, 4.042,
    "M3", 0.2107783512E-03, 8.280,
    "S8", 0.5817764322E-03, 3.000,
    NULL, 0.0, 0.0
};

TideStation *ReadTideStation(FILE * fp);

double GetFrequency(char *name)
{
    TidalConstituent *t = tc;
    while (t->name != NULL) {
	if (!strcmp(t->name, name)) {
	    return t->freq;
	}
	t++;
    }
    return -1.0;
}

double GetPeriod(char *name)
{
    TidalConstituent *t = tc;
    while (t->name != NULL) {
	if (!strcmp(t->name, name)) {
	    return t->period;
	}
	t++;
    }
    return -1.0;
}

TideStation *NewTideStation(int nfreq, int type)
{
    TideStation *ts;
    if (nfreq <= 0) {
	return NULL;
    }
    ts = (TideStation *) malloc(sizeof(TideStation));
    ts->omega = (double *) malloc(sizeof(double) * nfreq);
    ts->freqname = (char **) malloc(sizeof(char *) * nfreq);
    ts->nfreq = nfreq;
    ts->elamp = NULL;
    ts->elphase = NULL;
    ts->ampx = NULL;
    ts->phax = NULL;
    ts->ampy = NULL;
    ts->phay = NULL;
    switch (type) {
    case 0:			/* elevations only */
	ts->elamp = (double *) malloc(sizeof(double) * nfreq);
	ts->elphase = (double *) malloc(sizeof(double) * nfreq);
	break;
    case 1:			/* velocities only */
	ts->ampx = (double *) malloc(sizeof(double) * nfreq);
	ts->phax = (double *) malloc(sizeof(double) * nfreq);
	ts->ampy = (double *) malloc(sizeof(double) * nfreq);
	ts->phay = (double *) malloc(sizeof(double) * nfreq);
	break;
    case 2:			/* both */
	ts->elamp = (double *) malloc(sizeof(double) * nfreq);
	ts->elphase = (double *) malloc(sizeof(double) * nfreq);
	ts->ampx = (double *) malloc(sizeof(double) * nfreq);
	ts->phax = (double *) malloc(sizeof(double) * nfreq);
	ts->ampy = (double *) malloc(sizeof(double) * nfreq);
	ts->phay = (double *) malloc(sizeof(double) * nfreq);
	break;
    }
    return ts;
}

void DeleteTideStation(TideStation * ts)
{
    int i;
    if (ts->nfreq <= 0) {
	return;
    }
    free(ts->omega);
    for (i = 0; i < ts->nfreq; i++) {
	free(ts->freqname[i]);
    }
    free(ts->freqname);
    switch (ts->type) {
    case 0:			/* elevations only */
	free(ts->elamp);
	free(ts->elphase);
	break;
    case 1:			/* velocities only */
	free(ts->ampx);
	free(ts->phax);
	free(ts->ampy);
	free(ts->phay);
	break;
    case 2:			/* both */
	free(ts->elamp);
	free(ts->elphase);
	free(ts->ampx);
	free(ts->phax);
	free(ts->ampy);
	free(ts->phay);
	break;
    }
    free(ts);
}

int ReadTideStations(char *fname)
{
    char buf[256];
    int i, npts, type;
    TideStation *ts;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return 1;
    }
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    sscanf(buf, "%d %d", &npts, &type);
    tidestat = (TideStation **) malloc(npts * sizeof(TideStation *));
    ntidestat = npts;
    for (i = 0; i < npts; i++) {
	ts = ReadTideStation(fp);
	tidestat[i] = ts;
/* TODO  doesn't belong here */
	g[cg].tidestat[i].x = tidestat[i]->x;
	g[cg].tidestat[i].y = tidestat[i]->y;
	g[cg].tidestat[i].locx = tidestat[i]->x;
	g[cg].tidestat[i].locy = tidestat[i]->y;
    }
    update_tide(1);
    update_seltide();
    fclose(fp);
    return 1;
}

TideStation *ReadTideStation(FILE * fp)
{
    char buf[255];
    int i, j, nfreq, ftype, type, node;
    int freqtable = 0;		/* frequency from table */
    int unitstype = 1;		/* units for amplitudes */
    int datatype = 0;		/* elevations/velocities/elevations+velocities
				 * */
    int loctype = 0;		/* location type, xy or node */
    int phasetype = 0;		/* phase units, degrees or radians */
    char freqname[32];
    char name[256];
    char gname[256];
    char sfreqtable[32];
    char sunitstype[32];
    char sphasetype[32];
    char sdatatype[32];
    char sloctype[32];
    double freq, x, y;
    double amp, phase, ampx, phax, ampy, phay;
    TideStation *ts;
    fgets(buf, 255, fp);
    buf[strlen(buf) - 1] = 0;
    strcpy(name, buf);
    fgets(buf, 255, fp);
    sscanf(buf, "%d %d %d %d", &ftype, &nfreq, &type, &loctype);
/*
    sscanf(buf, "%d %s %s %s %s", &nfreq,
    	sfreqtable, sunitstype, sphasetype, sdatatype, sloctype);
    if (!strcmp(sfreqtable, "table")) {
	freqtable = 0;
    } else if (!strcmp(sfreqtable, "given")) {
	freqtable = 1;
    } else {
    }
    if (!strcmp(sdatatype, "elevations")) {
	datatype= 0;
    } else if (!strcmp(sdatatable, "velocities")) {
	datatype= 1;
    } else if (!strcmp(sdatatable, "elevations+velocities")) {
	datatype= 2;
    } else {
    }
    if (!strcmp(sunitstype, "cm")) {
	unitstype= 100;
    } else if (!strcmp(sfreqtable, "m")) {
	unitstype= 1;
    } else {
    }
    if (!strcmp(sphasetype, "degrees")) {
	phasetype= 0;
    } else if (!strcmp(sfreqtable, "radians")) {
	phasetype= 1;
    } else {
    }
    if (!strcmp(sloctype, "xy")) {
	loctype= 0;
    } else if (!strcmp(sloctype, "node")) {
	loctype= 1;
    } else {
    }
*/
    ts = NewTideStation(nfreq, type);
    strcpy(ts->name, name);
    ts->loctype = loctype;
    if (loctype == 0) {
	fgets(buf, 255, fp);
	sscanf(buf, "%lf %lf", &x, &y);
	ts->x = x;
	ts->y = y;
    } else {
	fgets(buf, 255, fp);
	sscanf(buf, "%d %s", &node, gname);
	ts->node = node;
	strcpy(ts->gname, gname);
    }
    switch (type) {
    case 0:			/* elevations only */
	for (j = 0; j < nfreq; j++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%s %lf %lf %lf", freqname, &freq, &amp, &phase);
	    ts->freqname[j] = (char *) malloc((strlen(freqname) + 1) * sizeof(char));
	    strcpy(ts->freqname[j], freqname);
	    ts->omega[j] = freq;
	    ts->elamp[j] = amp;
	    ts->elphase[j] = phase;
	}
	break;
    case 1:			/* velocities only */
	for (j = 0; j < nfreq; j++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%s %lf %lf %lf %lf %lf", freqname, &freq, &ampx, &phax, &ampy, &phay);
	    ts->freqname[j] = (char *) malloc((strlen(freqname) + 1) * sizeof(char));
	    strcpy(ts->freqname[j], freqname);
	    ts->omega[j] = freq;
	    ts->ampx[j] = ampx;
	    ts->phax[j] = phax;
	    ts->ampy[j] = ampy;
	    ts->phay[j] = phay;
	}
	break;
    case 2:			/* both */
	for (j = 0; j < nfreq; j++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%s %lf %lf %lf", freqname, &freq, &amp, &phase);
	    ts->freqname[j] = (char *) malloc((strlen(freqname) + 1) * sizeof(char));
	    strcpy(ts->freqname[j], freqname);
	    ts->omega[j] = freq;
	    ts->elamp[j] = amp;
	    ts->elphase[j] = phase;
	}
	for (j = 0; j < nfreq; j++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%s %lf %lf %lf %lf %lf", freqname, &freq, &ampx, &phax, &ampy, &phay);
	    ts->ampx[j] = ampx;
	    ts->phax[j] = phax;
	    ts->ampy[j] = ampy;
	    ts->phay[j] = phay;
	}
	break;
    }
    return ts;
}

int WriteTideStation(TideStation * ts, FILE * fp)
{
    char buf[255];
    int i, j;
    int type = 0;		/* elevations/velocities/elevations+velocities
				 * */
    int loctype = 0;		/* location type, xy or node */
    double x, y;
    fprintf(fp, "%s\n", ts->name);
    fprintf(fp, "%d %d %d %d\n", 0, ts->nfreq, ts->type, ts->loctype);
    if (loctype == 0) {
	fprintf(fp, "%lf %lf\n", ts->x, ts->y);
    } else {
	fprintf(fp, "%d %s\n", ts->node, ts->name);
    }
    switch (type) {
    case 0:			/* elevations only */
	for (j = 0; j < ts->nfreq; j++) {
	    fprintf(fp, "%s %lf %lf %lf\n", ts->freqname[j], ts->omega[j], ts->elamp[j], ts->elphase[j]);
	}
	break;
    case 1:			/* velocities only */
	for (j = 0; j < ts->nfreq; j++) {
	    fprintf(fp, "%s %lf %lf %lf %lf %lf\n", ts->freqname[j], ts->omega[j], ts->ampx[j], ts->phax[j], ts->ampy[j], ts->phay[j]);
	}
	break;
    case 2:			/* both */
	for (j = 0; j < ts->nfreq; j++) {
	    fprintf(fp, "%s %lf %lf %lf\n", ts->freqname[j], ts->omega[j], ts->elamp[j], ts->elphase[j]);
	}
	for (j = 0; j < ts->nfreq; j++) {
	    fprintf(fp, "%s %lf %lf %lf %lf %lf\n", ts->freqname[j], ts->omega[j], ts->ampx[j], ts->phax[j], ts->ampy[j], ts->phay[j]);
	}
	break;
    }
    return 0;
}

/*
 * evaluate elevations from tide station at time t
 * using constituent. If constituent == -1, then use all available
 * frequncies
 */
double EvalTideStationElev(TideStation * tides, double t, int constituent)
{
    int i;
    double e = 0.0;
    for (i = 0; i < tides->nfreq; i++) {
	e += tides->elamp[i] * cos(tides->omega[i] * t - tides->elphase[i]);
    }
    return e;
}

/*
 * evaluate velocities from tide station at time t
 * using constituent. If constituent == -1, then use all available
 * frequencies, return magnitude
 */
double EvalTideStationVel(TideStation * tides, double t, int constituent, double *u, double *v)
{
    double uu, vv;
    char buf[255];
    int i, npts, minc, maxc;
    double x, y;
    uu = vv = 0.0;
    if (constituent == -1) {	/* all components */
	minc = 0;
	maxc = tides->nfreq - 1;
    } else {
	minc = constituent;
	maxc = constituent;
    }
    for (i = minc; i < maxc; i++) {
	uu += tides->ampx[i] * cos(tides->omega[i] * t - tides->phax[i]);
	vv += tides->ampy[i] * cos(tides->omega[i] * t - tides->phay[i]);
    }
    return hypot(uu, vv);
}

void TideStationElevationMinMax(TideStation * tides, double *amin, double *amax)
{
    int i;
    double e = 0.0;
    for (i = 0; i < tides->nfreq; i++) {
	e += tides->elamp[i];
    }
    *amin = -e;
    *amax = e;
}

void find_nearest_tidestation(double wx, double wy, int *ind)
{
    int i;
    double dist, tmp, tmpx, tmpy;
    if (tidestat != NULL && ntidestat > 0) {
	tmpx = tidestat[0]->x - wx;
	tmpy = tidestat[0]->y - wy;
	dist = hypot(tmpx, tmpy);
	*ind = 0;
	for (i = 1; i < ntidestat; i++) {
	    tmpx = tidestat[i]->x - wx;
	    tmpy = tidestat[i]->y - wy;
	    tmp = hypot(tmpx, tmpy);
	    if (tmp < dist) {
		dist = tmp;
		*ind = i;
	    }
	}
    } else {
	*ind = -1;
    }
}
