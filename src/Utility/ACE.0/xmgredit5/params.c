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
 * write a parameter file
 *
 */

#ifndef lint
static char RCSid[] = "$Id: params.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";

#endif

#include <stdio.h>
#include "defines.h"
#include "globals.h"

#define PARMS_MAGIC 23456789

void putparms_binary(int gno, char *fname)
{
    int i, j, k, cnt, ming, maxg;
    FILE *pp;
    char buf[256];
    int magic = PARMS_MAGIC;
    if ((pp = fopen(fname, "w")) == NULL) {
	sprintf(buf, "Can't open parameter file %s", fname);
	errwin(buf);
	return;
    }
    fwrite(&magic, sizeof(int), 1, pp);

/* annotation */
    cnt = 0;
    for (k = 0; k < MAXSTR; k++) {
	if (pstr[k].active == ON && (pstr[k].s[0] || pstr[k].type)) {
	    cnt++;
	}
    }
    fwrite(&grid[0].ip, sizeof(Isolparms), 1, pp);
    fwrite(&grid[MAXGRIDS].ip, sizeof(Isolparms), 1, pp);
    fwrite(&cnt, sizeof(int), 1, pp);
    for (k = 0; k < MAXSTR; k++) {
	if (pstr[k].active == ON && (pstr[k].s[0] || pstr[k].type)) {
	    fwrite(&k, sizeof(int), 1, pp);
	    fwrite(&pstr[k], sizeof(plotstr), 1, pp);
	}
    }

    cnt = 0;
    for (k = 0; k < MAXBOXES; k++) {
	if (boxes[k].active == ON) {
	    cnt++;
	}
    }
    fwrite(&cnt, sizeof(int), 1, pp);
    for (k = 0; k < MAXBOXES; k++) {
	if (boxes[k].active == ON) {
	    fwrite(&k, sizeof(int), 1, pp);
	    fwrite(&boxes[k], sizeof(boxtype), 1, pp);
	}
    }

    cnt = 0;
    for (k = 0; k < MAXLINES; k++) {
	if (lines[k].active == ON) {
	    cnt++;
	}
    }
    fwrite(&cnt, sizeof(int), 1, pp);
    for (k = 0; k < MAXLINES; k++) {
	if (lines[k].active == ON) {
	    fwrite(&k, sizeof(int), 1, pp);
	    fwrite(&lines[k], sizeof(linetype), 1, pp);
	}
    }
    fclose(pp);
}

int getparms_binary(int gno, char *fname)
{
    int i, j, k, cnt, ng, itmp, magic;
    char buf[256];
    FILE *pp;
    if ((pp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "Can't open parameter file %s", fname);
	errwin(buf);
	return 0;
    }
    fread(&magic, sizeof(int), 1, pp);
    if (magic != PARMS_MAGIC) {
	errwin("Bad magic in parameter file");
	fclose(pp);
	return 0;
    }
    fread(&grid[0].ip, sizeof(Isolparms), 1, pp);
    fread(&grid[MAXGRIDS].ip, sizeof(Isolparms), 1, pp);

/* annotation */
    fread(&cnt, sizeof(int), 1, pp);
    for (i = 0; i < cnt; i++) {
	fread(&ng, sizeof(int), 1, pp);
	fread(&pstr[ng], sizeof(plotstr), 1, pp);
    }
    fread(&cnt, sizeof(int), 1, pp);
    for (i = 0; i < cnt; i++) {
	fread(&ng, sizeof(int), 1, pp);
	fread(&boxes[ng], sizeof(boxtype), 1, pp);
    }
    fread(&cnt, sizeof(int), 1, pp);
    for (i = 0; i < cnt; i++) {
	fread(&ng, sizeof(int), 1, pp);
	fread(&lines[ng], sizeof(linetype), 1, pp);
    }

    fclose(pp);
    return 1;
}
