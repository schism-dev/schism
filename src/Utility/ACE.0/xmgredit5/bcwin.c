/*
 * Boundary conditions
 */


#ifndef lint
static char RCSid[] =
    "$Id: bcwin.c,v 1.3 2003/11/04 17:14:42 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

int nopenb;
int adcirc_openb[100][2];
int nlandb;
int adcirc_landb[100][2];
int *adcirc_bound_ind;
extern int nopenb;
extern int adcirc_openb[][2];
extern int nlandb;
extern int adcirc_landb[][2];
extern int *adcirc_bound_ind;

int curadcbc = 0;

int check_boundary_overlap(int gridno);

void accept_adcircbound(Widget w, XtPointer cd);
void adcirc_defineland_proc(void);
void adcirc_defineopen_proc(void);
void adcirc_deleteopen_proc(void);
void create_adcircbound_frame(void);
void define_adcirc_landb(int gridno, int ind1, int ind2);
void define_adcirc_openb(int gridno, int ind1, int ind2);
void draw_adcirc_landb(int gridno);
void draw_adcirc_openb(int gridno);
void write_adcirc_openb(int gridno, FILE * fp);
void adcirc_writeb_proc(Widget w, XtPointer clientd, XtPointer calld);
void create_writeb_frame(void);

/*
 * popup to define boundary segments
 */
void create_adcircbound_frame(void)
{
    extern Widget app_shell;
    static Widget top;
    Widget rc, rc2, wbut;
    int x, y;

    if (!top) {
	XmGetPos(app_shell, 0, &x, &y);

	top =
	    XmCreateDialogShell(app_shell, "ADCIRC Open/Closed Boundaries",
				NULL, 0);
	handle_close(top);
	rc = XmCreateRowColumn(top, "rc", NULL, 0);

	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	wbut = XtVaCreateManagedWidget("Define new open boundary segment",
				       xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback,
		      (XtCallbackProc) adcirc_defineopen_proc,
		      (XtPointer) NULL);
	wbut =
	    XtVaCreateManagedWidget("Define new land boundary segment",
				    xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback,
		      (XtCallbackProc) adcirc_defineland_proc,
		      (XtPointer) NULL);
	wbut =
	    XtVaCreateManagedWidget("Clear all", xmPushButtonWidgetClass,
				    rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback,
		      (XtCallbackProc) adcirc_deleteopen_proc,
		      (XtPointer) NULL);
	wbut =
	    XtVaCreateManagedWidget("Write open/land boundary segments",
				    xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback,
		      (XtCallbackProc) create_writeb_frame,
		      (XtPointer) NULL);
	XtManageChild(rc2);

	rc2 = XmCreateRowColumn(rc, "rc", (XtPointer) NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, (XtPointer) NULL);
	wbut =
	    XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2,
				    (XtPointer) NULL);
	XtAddCallback(wbut, XmNactivateCallback,
		      (XtCallbackProc) create_writeb_frame,
		      (XtPointer) top);
	wbut =
	    XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc2,
				    (XtPointer) NULL);
	XtAddCallback(wbut, XmNactivateCallback,
		      (XtCallbackProc) destroy_dialog, top);
	XtManageChild(rc2);
	XtManageChild(rc);
    }
    XtRaise(top);
    compute_bound();
}

/*
 * select open boundary
 */
void adcirc_defineopen_proc(void)
{
    set_action(0);
    set_action(OPENB1ST);
}

/*
 * select land boundary
 */
void adcirc_defineland_proc(void)
{
    set_action(0);
    set_action(LANDB1ST);
}

/*
 * delete open boundary
 */
void adcirc_deleteopen_proc(void)
{
    nopenb = 0;
    nlandb = 0;
}

/*
 * write open boundary
 */
void adcirc_writeb_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    char *s, buf[1024];
    int i;
    FILE *fp;
    XmFileSelectionBoxCallbackStruct *cbs =
	(XmFileSelectionBoxCallbackStruct *) calld;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin
	    ("Error converting XmString to char string in adcirc_writeb_proc()");
	return;
    }
    strcpy(buf, s);
    XtFree(s);
    fp = fopen(buf, "w");
    if (fp != NULL) {
	write_adcirc_openb(0, fp);
	fclose(fp);
    } else {
	char b[1024];
	sprintf(b, "Unable to open file %s for writing", buf);
	errwin(b);
	return;
    }
}

/* read build points file with locations of start and end locations of open/land boundaries */
int file_bc(char *fname)
{
#define MAXBC 100
    char buf[256];
    FILE *fp;
    int i, n, bno, ires, itmp, ind1, ind2;
    double x[MAXBC], y[MAXBC], d;
    nlandb = nopenb = 0;
    if ((fp = fopen(fname, "r")) == NULL) {
	errwin("file_bc: unable to open file");
	return 0;
    }
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    sscanf(buf, "%d", &n);
    if (n > MAXBC) {
	errwin("file_bc: Too many BC segments");
	fclose(fp);
	return 0;
    }
    for (i = 0; i < n; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("file_bc: Not enough points in file");
	    fclose(fp);
	    return 0;
	}
	sscanf(buf, "%d %lf %lf %lf", &itmp, &x[i], &y[i], &d);
    }
    fclose(fp);
    for (i = 0; i < n; i++) {
	bno = grid[curgrid].boundaries[0];
	find_nearest_boundary_point2(bno, x[i], y[i], &ind1);
	find_nearest_boundary_point2(bno, x[(i + 1) % n], y[(i + 1) % n],
				     &ind2);
/* need to start with ocean boundary */
	if (i % 2 == 0) {
	    //diamond(boundary[bno].boundx[ind2], boundary[bno].boundy[ind2]);
	    //define_adcirc_openb(0, ind1, ind2);
	    adcirc_openb[nopenb][0] = ind1;
	    adcirc_openb[nopenb][1] = ind2;
	    nopenb++;
	    //draw_adcirc_openb(0);
	} else {
	    //diamond(boundary[bno].boundx[ind2], boundary[bno].boundy[ind2]);
	    //define_adcirc_landb(0, ind1, ind2);
	    //draw_adcirc_landb(0);
	    adcirc_landb[nlandb][0] = ind1;
	    adcirc_landb[nlandb][1] = ind2;
	    nlandb++;
	}
    }
}

/*
 * define an open boundary from selections made on the drawing area
 */
void define_adcirc_openb(int gridno, int ind1, int ind2)
{
    adcirc_openb[nopenb][0] = ind1;
    adcirc_openb[nopenb][1] = ind2;
    nopenb++;
    check_boundary_overlap(gridno);
}

/*
 * define a land boundary from selections made on the drawing area
 */
void define_adcirc_landb(int gridno, int ind1, int ind2)
{
    adcirc_landb[nlandb][0] = ind1;
    adcirc_landb[nlandb][1] = ind2;
    nlandb++;
    check_boundary_overlap(gridno);
}

/*
 * draw the open boundaries
 */
void draw_adcirc_openb(int gridno)
{
    int i, j, bno, n;
    bno = grid[gridno].boundaries[0];
    setcolor(4);
    for (j = 0; j < nopenb; j++) {
	if (adcirc_openb[j][0] > adcirc_openb[j][1]) {
	    my_move2(boundary[bno].boundx[adcirc_openb[j][0]],
		     boundary[bno].boundy[adcirc_openb[j][0]]);
	    for (i = adcirc_openb[j][0] + 1; i < boundary[bno].nbpts; i++) {
		my_draw2(boundary[bno].boundx[i], boundary[bno].boundy[i]);
	    }
	    for (i = 0; i <= adcirc_openb[j][1]; i++) {
		my_draw2(boundary[bno].boundx[i], boundary[bno].boundy[i]);
	    }
	} else {
	    my_move2(boundary[bno].boundx[adcirc_openb[j][0]],
		     boundary[bno].boundy[adcirc_openb[j][0]]);
	    for (i = adcirc_openb[j][0]; i <= adcirc_openb[j][1]; i++) {
		my_draw2(boundary[bno].boundx[i], boundary[bno].boundy[i]);
	    }
	}
    }
}

/*
 * draw the land boundaries
 */
void draw_adcirc_landb(int gridno)
{
    int i, j, bno, n;
    bno = grid[gridno].boundaries[0];
    setcolor(3);
    for (j = 0; j < nlandb; j++) {
	if (adcirc_landb[j][0] > adcirc_landb[j][1]) {
	    my_move2(boundary[bno].boundx[adcirc_landb[j][0]],
		     boundary[bno].boundy[adcirc_landb[j][0]]);
	    for (i = adcirc_landb[j][0] + 1; i < boundary[bno].nbpts; i++) {
		my_draw2(boundary[bno].boundx[i], boundary[bno].boundy[i]);
	    }
	    for (i = 0; i <= adcirc_landb[j][1]; i++) {
		my_draw2(boundary[bno].boundx[i], boundary[bno].boundy[i]);
	    }
	} else {
	    my_move2(boundary[bno].boundx[adcirc_landb[j][0]],
		     boundary[bno].boundy[adcirc_landb[j][0]]);
	    for (i = adcirc_landb[j][0]; i <= adcirc_landb[j][1]; i++) {
		my_draw2(boundary[bno].boundx[i], boundary[bno].boundy[i]);
	    }
	}
    }
}

int check_boundary_overlap(int gridno)
{
    int i, j, bno, n, cnt = 0;
    int *b, nl, nln, btmp, ino, inop, ind1, ind2, retval = 0;
    bno = grid[gridno].boundaries[0];
    /* setcolor(4); */
    /* array to mark visited nodes in the external boundary */
    b = (int *) calloc(boundary[bno].nbpts, sizeof(int));
    /* initialize (probably unnecessary) */
    for (i = 0; i < boundary[bno].nbpts; i++) {
	b[i] = 0;
    }
    for (j = 0; j < nopenb; j++) {
	b[adcirc_openb[j][0]] = 3;
	b[adcirc_openb[j][1]] = 3;
	if (adcirc_openb[j][0] > adcirc_openb[j][1]) {
	    for (i = adcirc_openb[j][0]; i < boundary[bno].nbpts; i++) {
		if (b[i] == 1) {
		    printf("Overlap at %d\n", boundary[bno].nodes[i] + 1);
		    retval++;
		} else if (b[i] == 0) {
		    b[i] = 1;
		}
	    }
	    for (i = 0; i <= adcirc_openb[j][1]; i++) {
		if (b[i] == 1) {
		    printf("Overlap at %d\n", boundary[bno].nodes[i] + 1);
		    retval++;
		} else {
		    b[i] = 1;
		}
	    }
	} else {
	    for (i = adcirc_openb[j][0]; i <= adcirc_openb[j][1]; i++) {
		if (b[i] == 1) {
		    printf("Overlap at %d\n", boundary[bno].nodes[i] + 1);
		    retval++;
		} else {
		    b[i] = 1;
		}
	    }
	}
    }
    for (j = 0; j < nlandb; j++) {
	b[adcirc_landb[j][0]] = 4;
	b[adcirc_landb[j][1]] = 4;
	if (adcirc_landb[j][0] > adcirc_landb[j][1]) {
	    for (i = adcirc_landb[j][0]; i < boundary[bno].nbpts; i++) {
		if (b[i] == 1) {
		    printf("Overlap with open boundary at %d\n",
			   boundary[bno].nodes[i] + 1);
		    retval++;
		} else {
		    b[i] = 2;
		}
	    }
	    for (i = 0; i <= adcirc_landb[j][1]; i++) {
		if (b[i] == 1) {
		    printf("Overlap with open boundary at %d\n",
			   boundary[bno].nodes[i] + 1);
		    retval++;
		} else {
		    b[i] = 2;
		}
	    }
	} else {
	    for (i = adcirc_landb[j][0]; i <= adcirc_landb[j][1]; i++) {
		if (b[i] == 1) {
		    printf("Overlap with open boundary at %d\n",
			   boundary[bno].nodes[i] + 1);
		    retval++;
		} else {
		    b[i] = 2;
		}
	    }
	}
    }
    free(b);
    return retval;
}

/*
 * draw the current state of adcirc selections TODO move out of here
 */
void draw_adcirc(void)
{
    if (nlandb) {
	draw_adcirc_landb(0);
	setcolor(1);
    }
    if (nopenb) {
	draw_adcirc_openb(0);
	setcolor(1);
    }
}

/*
 * Read a region popup
 */
void create_writeb_frame(void)
{
    static Widget top;
    int i;
    if (top) {
	XtRaise(top);
	return;
    }
    top =
	XmCreateFileSelectionDialog(app_shell, "Write open/land boundary",
				    NULL, 0);
    XtVaSetValues(XtParent(top), XmNtitle, "Write open/land boundary",
		  NULL);
    XtAddCallback(top, XmNcancelCallback, (XtCallbackProc) destroy_dialog,
		  (XtPointer) top);
    XtAddCallback(top, XmNokCallback, (XtCallbackProc) adcirc_writeb_proc,
		  (XtPointer) top);
    XtManageChild(top);
}

/*
 * write the open boundaries
 */
void write_adcirc_openb(int gridno, FILE * fp)
{
    int i, j, bno, n, cnt = 0;
    int nl, nln, btmp, ino, inop, ind1, ind2;
    bno = grid[gridno].boundaries[0];
    /* setcolor(4); */
    check_boundary_overlap(gridno);
    fprintf(fp, "%d = Number of open boundaries\n", nopenb);
    /* count number of open boundary nodes */
    for (j = 0; j < nopenb; j++) {
	if (adcirc_openb[j][0] > adcirc_openb[j][1]) {
	    cnt +=
		(boundary[bno].nbpts - adcirc_openb[j][0]) +
		adcirc_openb[j][1] + 1;
	} else {
	    cnt += adcirc_openb[j][1] - adcirc_openb[j][0] + 1;
	}
    }
    fprintf(fp, "%d = Total number of open boundary nodes\n", cnt);
    for (j = 0; j < nopenb; j++) {
	if (adcirc_openb[j][0] > adcirc_openb[j][1]) {
	    cnt =
		(boundary[bno].nbpts - adcirc_openb[j][0]) +
		adcirc_openb[j][1] + 1;
	} else {
	    cnt = adcirc_openb[j][1] - adcirc_openb[j][0] + 1;
	}
	fprintf(fp, "%d = Number of nodes for open boundary %d\n", cnt,
		j + 1);
	if (adcirc_openb[j][0] > adcirc_openb[j][1]) {
	    for (i = adcirc_openb[j][0]; i < boundary[bno].nbpts; i++) {
		fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	    }
	    for (i = 0; i <= adcirc_openb[j][1]; i++) {
		fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	    }
	} else {
	    for (i = adcirc_openb[j][0]; i <= adcirc_openb[j][1]; i++) {
		fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	    }
	}
    }
    cnt = 0;
    fprintf(fp, "%d = number of land boundaries\n",
	    nlandb + grid[gridno].nbounds - 1);
    /* count number of land boundary nodes */
    for (j = 0; j < nlandb; j++) {
	if (adcirc_landb[j][0] > adcirc_landb[j][1]) {
	    cnt +=
		(boundary[bno].nbpts - adcirc_landb[j][0]) +
		adcirc_landb[j][1] + 1;
	} else {
	    cnt += adcirc_landb[j][1] - adcirc_landb[j][0] + 1;
	}
    }
    for (j = 1; j < grid[gridno].nbounds; j++) {
	btmp = grid[gridno].boundaries[j];
	cnt += boundary[btmp].nbpts;
    }
    fprintf(fp, "%d = Total number of land boundary nodes\n", cnt);
    for (j = 0; j < nlandb; j++) {
	if (adcirc_landb[j][0] > adcirc_landb[j][1]) {
	    cnt =
		(boundary[bno].nbpts - adcirc_landb[j][0]) +
		adcirc_landb[j][1] + 1;
	} else {
	    cnt = adcirc_landb[j][1] - adcirc_landb[j][0] + 1;
	}
	fprintf(fp, "%d 0 = Number of nodes for land boundary %d\n", cnt,
		j + 1);
	if (adcirc_landb[j][0] > adcirc_landb[j][1]) {
	    for (i = adcirc_landb[j][0]; i < boundary[bno].nbpts; i++) {
		fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	    }
	    for (i = 0; i <= adcirc_landb[j][1]; i++) {
		fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	    }
	} else {
	    for (i = adcirc_landb[j][0]; i <= adcirc_landb[j][1]; i++) {
		fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	    }
	}
    }
    for (j = 1; j < grid[gridno].nbounds; j++) {
	bno = grid[gridno].boundaries[j];
	fprintf(fp, "%d 1 = Number of nodes for island boundary %d\n",
		boundary[bno].nbpts, j);
	for (i = 0; i < boundary[bno].nbpts; i++) {
	    fprintf(fp, "%d\n", boundary[bno].nodes[i] + 1);
	}
    }
}
