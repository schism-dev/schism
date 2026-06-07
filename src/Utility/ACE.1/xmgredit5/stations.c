/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2004 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/*
 * utilities for stations
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

/* Station *sta;
 * int nsta;
 */

/* 
 * typedef struct _Station { 
     int active;
     int type;
     char label[256];         
     double x;
     double y;                 
     int display;               
     Props p;
     Isolparms ip;               
 } Station; 
*/

int ReadStations(char *fname, int *n, Station **s);

Station *NewStations(int n)
{
    Station *s;
    s = (Station *) malloc(n * sizeof(Station));
    return s;
}

void DeleteStations(int *n, Station *s)
{
    *n = 0;
    free(s);
}

int ReadStations(char *fname, int *n, Station **s)
{
    char buf[256], label[256];
    int i, nn, type;
    Station *stmp;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return 1;
    }
    fgets(buf, 255, fp);
    sscanf(buf, "%d", &nn);
    stmp = (Station *) malloc(nn * sizeof(Station));
    for (i = 0; i < nn; i++) {
        fgets(buf, 255, fp);
        strcpy(stmp[i].label, buf);
	stmp[i].label[strlen(buf) - 1] = 0;
        fgets(buf, 255, fp);
        sscanf(buf, "%lf %lf", &stmp[i].x, &stmp[i].y);
	//printf("%d, %s: %lf %lf\n", i, stmp[i].label, stmp[i].x, stmp[i].y); 
	stmp[i].type = 0;
	stmp[i].active = ON;
	stmp[i].p.color = 1;
	stmp[i].p.symbol = 2;
    }
    fclose(fp);
    *n = nn;
    *s = stmp;
    return 0;
}

int ReadBuildPoints(char *fname, int *n, Station **s)
{
    char buf[256], label[256];
    int i, nn, type;
    Station *stmp;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return 1;
    }
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    sscanf(buf, "%d", &nn);
    stmp = (Station *) malloc(nn * sizeof(Station));
    for (i = 0; i < nn; i++) {
        sprintf(stmp[i].label, "%d", i + 1);
        fgets(buf, 255, fp);
        sscanf(buf, "%*d %lf %lf", &stmp[i].x, &stmp[i].y);
	//printf("%d, %s: %lf %lf\n", i, stmp[i].label, stmp[i].x, stmp[i].y); 
	stmp[i].type = 0;
	stmp[i].active = ON;
	stmp[i].p.color = 1;
	stmp[i].p.symbol = 2;
    }
    fclose(fp);
    *n = nn;
    *s = stmp;
    return 0;
}

void DrawStations(int n, Station *s)
{
    int i;
    char buf[256];
    for (i=0;i<n;i++) {
	setcolor(s[i].p.color);
	drawpolysym(&s[i].x, &s[i].y, 1, s[i].p.symbol, 0, 1, 0.8);
        buf[0] = ' '; buf[1] = 0;
        strcat(buf, sta[i].label);
	writestr(s[i].x, s[i].y, 0, 0, buf);
    }
}

void DrawStation(int gno, int n)
{
    char buf[256];
    int ix, iy;
    extern int devheight, docoords;
    buf[0] = ' ';
    buf[1] = 0;
    strcat(buf, sta[n].label);
    setcolor(sta[n].p.color);
    drawpolysym(&sta[n].x, &sta[n].y, 1, sta[n].p.symbol, 0, 1, 0.8);
    writestr(sta[n].x, sta[n].y, 0, 0, buf);
    world2deviceabs(sta[n].x, sta[n].y, &ix, &iy);
/*
    if (hardcopyflag && docoords) {
        printf("%s: %lf %lf %d %d\n", buf, sta[n].x, sta[n].y, xconvgd(sta[n].x), devheight - yconvgd(sta[n].y));
    }
*/
}

void FindNearestStation(int n, Station *s, double x, double y, int *ind)
{
    int i;
    double tmp, cmin = 1e307;
    *ind = -1;
    for (i = 0; i < n; i++) {
	if (s[i].active == ON) {
	    tmp = hypot(x - s[i].x, y - s[i].y);
	    if (cmin > tmp) {
		cmin = tmp;
		*ind = i;
	    }
	}
    }
}
