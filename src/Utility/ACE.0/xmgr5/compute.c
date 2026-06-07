/* $Id: compute.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * perform math between sets
 *
 */

#include <stdio.h>
#include "globals.h"
#include "noxprotos.h"

void loadset(int gno, int selset, int toval, double startno, double stepno)
{
    int i, lenset;
    double *ltmp;
    double *xtmp, *ytmp;

    if ((lenset = getsetlength(gno, selset)) <= 0) {
	char stmp[60];

	sprintf(stmp, "Length of set %d <= 0", selset);
	errwin(stmp);
	return;
    }
    xtmp = getx(gno, selset);
    ytmp = gety(gno, selset);
    switch (toval) {
    case 1:
	ltmp = xtmp;
	break;
    case 2:
	ltmp = ytmp;
	break;
    case 3:
	ltmp = ax;
	break;
    case 4:
	ltmp = bx;
	break;
    case 5:
	ltmp = cx;
	break;
    case 6:
	ltmp = dx;
	break;
    default:
	return;
    }
    for (i = 0; i < lenset; i++) {
	*ltmp++ = startno + i * stepno;
    }
    updatesetminmax(gno, selset);
    update_set_status(gno, selset);
}

/*
 * evaluate the expression in sscanstr and place the result in selset
 */
int formula(int gno, int selset, char *sscanstr)
{
    char stmp[64], tmpstr[512];
    int i = 0, errpos, lenset;
    double *xtmp, *ytmp;

    if ((lenset = getsetlength(gno, selset)) <= 0) {
	sprintf(stmp, "Length of set %d = 0", selset);
	errwin(stmp);
	return;
    }
    xtmp = getx(gno, selset);
    ytmp = gety(gno, selset);
    strcpy(tmpstr, sscanstr);
    fixupstr(tmpstr);
    scanner(tmpstr, xtmp, ytmp, lenset, ax, bx, cx, dx, MAXARR, i, selset, &errpos);
    updatesetminmax(gno, selset);
    update_set_status(gno, selset);
    return (errpos);
}
