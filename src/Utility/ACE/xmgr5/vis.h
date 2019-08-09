/* $Id: vis.h,v 1.1.1.1 2003/07/21 16:18:42 pturner Exp $
 *
 */
/************************************************************************
*
*	OREGON GRADUATE INSTITUTE
*	Center for Coastal and Land Margin Research
*
*************************************************************************
*
*	ACE1.C - Provide definitions and structures for the ace1 model
*                user interfaces.
*
*	Purpose: Structures are defined.
*
*	Method:  White book C conventions are used.
*
*	History: Version 1.0	October 15, 1991	P. J. Turner
*							A. M. Baptista
*                                                       J. R. Hurst
*
*	Copyright:
*
*		Copyright 1991, Center for Coastal and Land Margin Research
*				Oregon Graduate Institute
*				Beaverton, OR  97006
*
*************************************************************************
*/

#define MAXVIS   20		/* max number of visualization sets */
#define MEM   9999		/* max number of visualization sets */

typedef struct {
    int active;
    int type;
    char fname[256];
    FILE *fp;
    int ftype;
    int src;
    int set1, set2;
    int gno;
    int skip;
    int nsteps;
    int npts;
    int nsets;
    int sets[10];
    int slen[10];
    double *x, *y, *z;
    double *tm, **xm, **ym, **zm;
    float *xf, *yf, *zf;
    double cur_time;
    int read_data;
    int display_data;
    plotstr vs;
    int model_flag;
} vis_struct;
