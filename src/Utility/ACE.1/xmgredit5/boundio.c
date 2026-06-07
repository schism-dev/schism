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

 * Read/write boundaries
 *
 */

#ifndef lint
static char RCSid[] = "$Id: boundio.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "globals.h"

/* boundio.c */
int writeboundary(int gridno, char *fname, int type);
int readboundary(int gridno, char *fname, int val);
int readbinboundary(int gridno, char *fname, int val);
int writebinboundary(int gridno, char *fname, int val);
Boundaries *ReadBoundary(char *fname);

/*
 * Write the external and internal boundaries of a grid
 */
int writeboundary(int gridno, char *fname, int type)
{
    int i, j, nb, bno, ntmp;
    FILE *fp;

    if ((fp = fopen(fname, "w")) == NULL) {
	errwin("writeboundary: can't open file for writing");
	return 0;
    }
    if (grid[gridno].nbounds == 0) {
	errwin("writeboundary: No boundaries to write");
	return 0;
    }
    bno = grid[gridno].boundaries[0];
/*
   if (boundary[bno].boundtype != 0) {
   errwin("writeboundary: Assumption about 1st boundary incorrect");
   return 0;
   }
 */
    fprintf(fp, "%s\n", fname);
    fprintf(fp, "%d\n", grid[gridno].nbounds);
    for (i = 0; i < grid[gridno].nbounds; i++) {
	bno = grid[gridno].boundaries[i];
	fprintf(fp, "%d 1\n", boundary[bno].nbpts);
	for (j = 0; j < boundary[bno].nbpts; j++) {
	    if (type == 0) {
		fprintf(fp, "%d\n", boundary[bno].nodes[j] + 1);
	    } else if (type == 1) {
		fprintf(fp, "%lf %lf\n", boundary[bno].boundx[j], boundary[bno].boundy[j]);
	    } else if (type == 2) {
		ntmp = boundary[bno].nodes[j];
		fprintf(fp, "%lf %lf %lf\n", grid[gridno].xord[ntmp], grid[gridno].yord[ntmp], grid[gridno].depth[ntmp]);
	    }
	}
    }
    fclose(fp);
    return 1;
}

/*
 * Read the external and internal boundaries of a grid
 */
int readboundary(int gridno, char *fname, int val)
{
    int i, j, nb, npts, itmp;
    FILE *fp;
    char s[255];
    extern int ebound_defined, ibound_defined;

    if ((fp = fopen(fname, "r")) == NULL) {
	errwin("readboundary: can't open file for reading");
	return 0;
    }
    Free_boundaries(gridno);
    fgets(s, 128, fp);
    fgets(s, 128, fp);
    sscanf(s, "%d", &grid[gridno].nbounds);
    for (i = 0; i < grid[gridno].nbounds; i++) {
	nb = nextboundary();
	if (nb == -1) {
	    errwin("readboundary: Not enough boundaries");
	    fclose(fp);
	    return 0;
	}
	grid[gridno].boundaries[i] = nb;
	fgets(s, 128, fp);
	sscanf(s, "%d", &npts);
	Allocate_boundary(nb, npts, (i != 0));
	for (j = 0; j < boundary[nb].nbpts; j++) {
	    fgets(s, 128, fp);
	    if (val) {
		sscanf(s, "%d", &itmp);
		boundary[nb].boundx[j] = grid[gridno].xord[itmp - 1];
		boundary[nb].boundy[j] = grid[gridno].yord[itmp - 1];
	    } else {
		sscanf(s, "%lf %lf", &boundary[nb].boundx[j], &boundary[nb].boundy[j]);
	    }
	}
    }
    if (grid[gridno].nbounds != 0) {
	if (gridno < MAXGRIDS) {
	    if (grid[gridno].nbounds > 1) {
		ebound_defined = 1;
		ibound_defined = 1;
	    } else {
		ebound_defined = 1;
	    }
	} else {
	}
    }
    fclose(fp);
    check_boundaries(gridno);
    update_boundary_menu(gridno);
    return 1;
}

/*
 * Read a boundary in binary format
 */
int readbinboundary(int gridno, char *fname, int val)
{
    int i, j, n, nb, ib, npts, itmp, magic;
    FILE *fp;
    char s[255];
    float *xbtmp, *ybtmp;
    double *xb, *yb;
    extern int ebound_defined, ibound_defined;
    if ((fp = fopen(fname, "r")) == NULL) {
	errwin("readboundary: can't open file for reading");
	return 0;
    }
/*
   fread(&magic, 1, sizeof(int), fp);
   if (magic != 32) {
   errwin("readbinboundary: Bad magic in boundary file");
   return 0;
   }
 */
    Free_boundaries(gridno);
    fread(&nb, 1, sizeof(int), fp);
    grid[gridno].nbounds = nb;
    for (i = 0; i < nb; i++) {
	fread(&n, 1, sizeof(int), fp);
	ib = nextboundary();
	if (ib == -1) {
	    errwin("readboundary: Not enough boundaries");
	    fclose(fp);
	    return 0;
	}
	grid[gridno].boundaries[i] = ib;
	Allocate_boundary(ib, n, (i != 0));
	xbtmp = (float *) calloc(n, sizeof(float));
	ybtmp = (float *) calloc(n, sizeof(float));
	fread(xbtmp, sizeof(float), n, fp);
	fread(ybtmp, sizeof(float), n, fp);
	xb = boundary[ib].boundx;
	yb = boundary[ib].boundy;
	boundary[ib].boundtype = 3;
	for (j = 0; j < n; j++) {
	    xb[j] = xbtmp[j];
	    yb[j] = ybtmp[j];
	}
	free(xbtmp);
	free(ybtmp);
    }
    fclose(fp);
    if (grid[gridno].nbounds != 0) {
	if (gridno < MAXGRIDS) {
	    if (grid[gridno].nbounds > 1) {
		ebound_defined = 1;
		ibound_defined = 1;
	    } else {
		ebound_defined = 1;
	    }
	} else {
	}
    }
    check_boundaries(gridno);
    update_boundary_menu(gridno);
    return 1;
}

/*
 * Write a boundary in binary format
 */
int writebinboundary(int gridno, char *fname, int val)
{
    int i, j, n, nb, ib, npts, itmp, magic = 32;
    FILE *fp;
    char s[255];
    float *xbtmp, *ybtmp;
    double *xb, *yb;
    if (grid[gridno].nbounds == 0) {
	errwin("writeboundary: No boundaries to write");
	return 0;
    }
    if ((fp = fopen(fname, "w")) == NULL) {
	errwin("writebinboundary: can't open file for writing");
	return 0;
    }
/*
   fwrite(&magic, 1, sizeof(int), fp);
 */
    nb = grid[gridno].nbounds;
    fwrite(&nb, 1, sizeof(int), fp);
    for (i = 0; i < nb; i++) {
	ib = grid[gridno].boundaries[i];
	fwrite(&boundary[ib].nbpts, sizeof(int), 1, fp);
	n = boundary[ib].nbpts;
	xbtmp = (float *) calloc(boundary[ib].nbpts, sizeof(float));
	ybtmp = (float *) calloc(boundary[ib].nbpts, sizeof(float));
	xb = boundary[ib].boundx;
	yb = boundary[ib].boundy;
	for (j = 0; j < n; j++) {
	    xbtmp[j] = xb[j];
	    ybtmp[j] = yb[j];
	}
	fwrite(xbtmp, sizeof(float), n, fp);
	fwrite(ybtmp, sizeof(float), n, fp);
	free(xbtmp);
	free(ybtmp);
    }
    fclose(fp);
    return 1;
}

#ifdef HAVE_NETCDF
/*
 * netcdf file definition
 * netcdf boundary {
 * dimensions:
 *         npoints = 30495 ; // total number of points
 *         npolys = 54 ; // total number of closed polygons
 * variables:
 *         double polyx(npoints) ;
 *         double polyy(npoints) ;
 *         int polynodes(npoints); // boundary as nodes
 *         int polys(npolys) ; // number of nodes in each poly
 * }
 */

#include "netcdf.h"

/* Create the netcdf file */
void CreateBoundarynetcdf(char *fname, int gridno)
{
    int i, j;
    int nb, ib, n;
    int ncid;			/* netCDF id */

    /* dimension ids */
    int npoints_dim;
    int npolys_dim;

    /* dimension lengths */
    size_t npoints_len;
    size_t npolys_len;

    /* variable ids */
    int boundx_id;
    int boundy_id;
    int boundnodes_id;
    int polys_id;

    /* variable shapes */
    int boundx_dims[1];
    int boundy_dims[1];
    int boundnodes_dims[1];
    int polys_dims[1];

    int stat;

    int sum = 0;
    /* 
     * find the total number of points in all boundaries
     */
    nb = grid[gridno].nbounds;
    for (i = 0; i < nb; i++) {
	ib = grid[gridno].boundaries[i];
	sum += boundary[ib].nbpts;
    }
    npoints_len = sum;
    npolys_len = nb;

    /* enter define mode */
    stat = nc_create(fname, NC_CLOBBER, &ncid);
    check_err(stat, __LINE__, __FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "npoints", npoints_len, &npoints_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_dim(ncid, "npolys", npolys_len, &npolys_dim);
    check_err(stat, __LINE__, __FILE__);

    /* define variables */
    boundx_dims[0] = npoints_dim;
    stat = nc_def_var(ncid, "boundx", NC_DOUBLE, 1, boundx_dims, &boundx_id);
    check_err(stat, __LINE__, __FILE__);

    boundy_dims[0] = npoints_dim;
    stat = nc_def_var(ncid, "boundy", NC_DOUBLE, 1, boundy_dims, &boundy_id);
    check_err(stat, __LINE__, __FILE__);

    boundnodes_dims[0] = npoints_dim;
    stat = nc_def_var(ncid, "boundnodes", NC_INT, 1, boundnodes_dims, &boundnodes_id);
    check_err(stat, __LINE__, __FILE__);

    polys_dims[0] = npolys_dim;
    stat = nc_def_var(ncid, "polys", NC_INT, 1, polys_dims, &polys_id);
    check_err(stat, __LINE__, __FILE__);

    /* leave define mode */
    stat = nc_enddef(ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

/* Write the grid to the netcdf file fname */
void WriteBoundarynetcdf(char *fname, int gridno)
{
    int ib, i, j, cnt, *polys, *nodes;
    double *tmpx, *tmpy;

    int ncid;			/* netCDF id */

    /* dimension ids */
    int npoints_dim;
    int npolys_dim;

    /* dimension lengths */
    size_t npoints_len;
    size_t npolys_len;
    size_t start[2], count[2];

    /* variable ids */
    int boundx_id;
    int boundy_id;
    int boundnodes_id;
    int polys_id;

    int stat;

    stat = nc_open(fname, NC_WRITE, &ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "npolys", &npolys_dim);
    stat = nc_inq_dim(ncid, npolys_dim, "npolys", &npolys_len);
    stat = nc_inq_dimid(ncid, "npoints", &npoints_dim);
    stat = nc_inq_dim(ncid, npoints_dim, "npoints", &npoints_len);
    stat = nc_inq_varid(ncid, "boundx", &boundx_id);
    stat = nc_inq_varid(ncid, "boundy", &boundy_id);
    stat = nc_inq_varid(ncid, "boundnodes", &boundnodes_id);
    stat = nc_inq_varid(ncid, "polys", &polys_id);

    polys = (int *) malloc(npolys_len * sizeof(int));
    if (polys == NULL) {
	stat = nc_close(ncid);
	return;
    }
    nodes = (int *) malloc(npoints_len * sizeof(int));
    if (nodes == NULL) {
	stat = nc_close(ncid);
	return;
    }
    tmpx = (double *) malloc(npoints_len * sizeof(double));
    if (tmpx == NULL) {
	stat = nc_close(ncid);
	return;
    }
    tmpy = (double *) malloc(npoints_len * sizeof(double));
    if (tmpy == NULL) {
	free(tmpx);
	stat = nc_close(ncid);
	return;
    }
    cnt = 0;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	ib = grid[gridno].boundaries[i];
	polys[i] = boundary[ib].nbpts;
	for (j = 0; j < boundary[ib].nbpts; j++) {
	    nodes[cnt] = boundary[ib].nodes[j];
	    tmpx[cnt] = boundary[ib].boundx[j];
	    tmpy[cnt++] = boundary[ib].boundy[j];
	}
    }

    start[0] = 0;
    count[0] = cnt;
    stat = nc_put_vara_int(ncid, boundnodes_id, start, count, nodes);
    stat = nc_put_vara_double(ncid, boundx_id, start, count, tmpx);
    stat = nc_put_vara_double(ncid, boundy_id, start, count, tmpy);
    start[0] = 0;
    count[0] = npolys_len;
    stat = nc_put_vara_int(ncid, polys_id, start, count, polys);

    free(nodes);
    free(tmpx);
    free(tmpy);
    free(polys);

    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

/* Read the grid from the netcdf file fname */
void ReadBoundarynetcdf(char *fname, int gridno)
{
    int ib, i, j, cnt;
    int *polys, *nodes;
    double *tmpx, *tmpy;

    int ncid;			/* netCDF id */

    /* dimension ids */
    int npoints_dim;
    int npolys_dim;

    /* dimension lengths */
    size_t npoints_len;
    size_t npolys_len;
    size_t start[1], count[1];

    /* variable ids */
    int boundx_id;
    int boundy_id;
    int boundnodes_id;
    int polys_id;

    int stat;
    Free_boundaries(gridno);

    stat = nc_open(fname, NC_NOWRITE, &ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "npolys", &npolys_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dim(ncid, npolys_dim, "npolys", &npolys_len);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "npoints", &npoints_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dim(ncid, npoints_dim, "npoints", &npoints_len);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_inq_varid(ncid, "boundx", &boundx_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "boundy", &boundy_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "boundnodes", &boundnodes_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "polys", &polys_id);
    check_err(stat, __LINE__, __FILE__);

    polys = (int *) malloc(npolys_len * sizeof(int));
    if (polys == NULL) {
	stat = nc_close(ncid);
	return;
    }
    nodes = (int *) malloc(npoints_len * sizeof(int));
    if (nodes == NULL) {
	stat = nc_close(ncid);
	return;
    }
    tmpx = (double *) malloc(npoints_len * sizeof(double));
    if (tmpx == NULL) {
	free(polys);
	stat = nc_close(ncid);
	return;
    }
    tmpy = (double *) malloc(npoints_len * sizeof(double));
    if (tmpy == NULL) {
	free(polys);
	free(tmpx);
	stat = nc_close(ncid);
	return;
    }
    start[0] = 0;
    count[0] = npoints_len;
    stat = nc_get_vara_double(ncid, boundx_id, start, count, tmpx);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, boundy_id, start, count, tmpy);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_int(ncid, boundnodes_id, start, count, nodes);
    check_err(stat, __LINE__, __FILE__);
    start[0] = 0;
    count[0] = npolys_len;
    stat = nc_get_vara_int(ncid, polys_id, start, count, polys);
    check_err(stat, __LINE__, __FILE__);

    Free_boundaries(gridno);
    grid[gridno].nbounds = npolys_len;
    cnt = 0;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	ib = nextboundary();
	if (ib == -1) {
	    errwin("ReadBoundarynetcdf(): Not enough boundaries");
	    nc_close(ncid);
	    return;
	}
	grid[gridno].boundaries[i] = ib;
	Allocate_boundary(ib, polys[i], (i != 0));
	boundary[ib].boundtype = 3;
	for (j = 0; j < boundary[ib].nbpts; j++) {
	    boundary[ib].boundx[j] = tmpx[cnt + j];
	    boundary[ib].boundy[j] = tmpy[cnt + j];
	    boundary[ib].nodes[j] = nodes[cnt + j];
	}
	cnt += polys[i];
    }

    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
    check_boundaries(gridno);
    update_boundary_menu(gridno);
}

#endif
