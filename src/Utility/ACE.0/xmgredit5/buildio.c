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
 * Read/write build points in ASCII or binary.
 *
 */

#ifndef lint
static char RCSid[] = "$Id: buildio.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "globals.h"

/* buildio.c */
int readbuild(int bno, char *fname);
int writebuild(int bno, char *fname);
int readbuildbinary(int bno, char *fname);
int writebuildbinary(int bno, char *fname);

extern int build_defined;

/*
 * Read a set of build points to build points struct bno
 */
int readbuild(int bno, char *fname)
{
    int i, ires, itmp;
    char buf[256];
    FILE *fp;
    extern int noverride;

    if ((fp = fopen(fname, "r")) == NULL) {
	errwin("readbuild: unable to open file");
	return 0;
    }
    Free_build(bno);
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    if (noverride) {
	build[bno].nbuild = noverride;
    } else {
	sscanf(buf, "%d", &build[bno].nbuild);
    }

    Allocate_build(bno, build[bno].nbuild);
    for (i = 0; i < build[bno].nbuild; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("readbuild: Not enough points in file");
	    Free_build(bno);
	    fclose(fp);
	    return 0;
	}
	if ((ires = sscanf(buf, "%d %lf %lf %lf", &itmp, &build[bno].bx[i], &build[bno].by[i], &build[bno].db[i]) == 3)) {
	    build[bno].db[i] = 1.0;
	}
    }
    fclose(fp);
    if (build[bno].nbuild) {
	build_defined = 1;
    } else {
	build_defined = 0;
    }
    if (inwin) {
	update_fuzz_items();
    }
    return 1;
}

int writebuild(int bno, char *fname)
{
    int i;
    char buf[256];
    FILE *fp;

    if ((fp = fopen(fname, "w")) == NULL) {
	errwin("writebuild: unable to open file");
	return 0;
    }
    fprintf(fp, "%s\n", fname);
    fprintf(fp, "%d\n", build[bno].nbuild);
    for (i = 0; i < build[bno].nbuild; i++) {
	fprintf(fp, "%d %lf %lf %lf\n", i + 1, build[bno].bx[i], build[bno].by[i], build[bno].db[i]);
    }
    fclose(fp);
    return 1;
}

int readbuildbinary(int bno, char *fname)
{
    char buf[256];
    int i, npts;
    int testmagic, magic = 40;
    float *ftmp;
    FILE *fp;

    if ((fp = fopen(fname, "w")) == NULL) {
	errwin("writebuild: unable to open file");
	return 0;
    }
    fread(&testmagic, sizeof(int), 1, fp);
    if (testmagic != magic) {
	fclose(fp);
	return 0;
    }
    fread(&npts, sizeof(int), 1, fp);
    Free_build(bno);
    Allocate_build(bno, npts);
    ftmp = (float *) malloc(npts * sizeof(float));
    fread(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < npts; i++) {
	build[bno].bx[i] = ftmp[i];
    }
    fread(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < build[bno].nbuild; i++) {
	build[bno].by[i] = ftmp[i];
    }
    fread(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < build[bno].nbuild; i++) {
	build[bno].db[i] = ftmp[i];
    }
    free(ftmp);
    fclose(fp);
    if (build[bno].nbuild) {
	build_defined = 1;
    } else {
	build_defined = 0;
    }
    if (inwin) {
	update_fuzz_items();
    }
    return 1;
}

int writebuildbinary(int bno, char *fname)
{
    int i, npts;
    char buf[256];
    FILE *fp;
    float *ftmp;
    int magic = 40;

    if ((fp = fopen(fname, "w")) == NULL) {
	errwin("writebuild: unable to open file");
	return 0;
    }
    npts = build[bno].nbuild;
    ftmp = (float *) malloc(build[bno].nbuild * sizeof(float));
    fwrite(&magic, sizeof(int), 1, fp);
    fwrite(&build[bno].nbuild, sizeof(int), 1, fp);
    for (i = 0; i < build[bno].nbuild; i++) {
	ftmp[i] = build[bno].bx[i];
    }
    fwrite(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < build[bno].nbuild; i++) {
	ftmp[i] = build[bno].by[i];
    }
    fwrite(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < build[bno].nbuild; i++) {
	ftmp[i] = build[bno].db[i];
    }
    fwrite(ftmp, sizeof(float), npts, fp);
    free(ftmp);
    fclose(fp);
    return 1;
}

#ifdef HAVE_NETCDF
/*
 * netcdf file definition
 * netcdf build {
 * dimensions:
 *         nbuild = 30495 ; // total number of points
 * variables:
 *         double buildx(nbuild);
 *         double buildy(nbuild);
 *         double buildz(nbuild);
 * }
 */

#include "netcdf.h"

/* Create the netcdf file */
void CreateBuildnetcdf(char *fname, int bno)
{
    int i;
    int ncid;			/* netCDF id */

    /* dimension ids */
    int nbuild_dim;

    /* dimension lengths */
    size_t nbuild_len = build[bno].nbuild;

    /* variable ids */
    int x_id;
    int y_id;
    int z_id;

    /* variable shapes */
    int x_dims[1];
    int y_dims[1];
    int z_dims[1];
   
    int stat;

    /* enter define mode */
    stat = nc_create(fname, NC_CLOBBER, &ncid);
    check_err(stat, __LINE__, __FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "nbuild", nbuild_len, &nbuild_dim);
    check_err(stat, __LINE__, __FILE__);

    /* define variables */
    x_dims[0] = nbuild_dim;
    stat = nc_def_var(ncid, "buildx", NC_DOUBLE, 1, x_dims, &x_id);
    check_err(stat, __LINE__, __FILE__);

    y_dims[0] = nbuild_dim;
    stat = nc_def_var(ncid, "buildy", NC_DOUBLE, 1, y_dims, &y_id);
    check_err(stat, __LINE__, __FILE__);

    z_dims[0] = nbuild_dim;
    stat = nc_def_var(ncid, "buildz", NC_DOUBLE, 1, z_dims, &z_id);
    check_err(stat, __LINE__, __FILE__);

    /* leave define mode */
    stat = nc_enddef(ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

/* Write the grid to the netcdf file fname */
void WriteBuildnetcdf(char *fname, int bno)
{
    int i, j, cnt;
    int *ellist;

    int ncid;			/* netCDF id */

    /* dimension ids */
    int nbuild_dim;

    /* dimension lengths */
    size_t nbuild_len;
    size_t start[1], count[1];

    /* variable ids */
    int x_id;
    int y_id;
    int z_id;

    int stat;

    stat = nc_open(fname, NC_WRITE, &ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "nbuild", &nbuild_dim);
    stat = nc_inq_dim(ncid, nbuild_dim, "nbuild", &nbuild_len);
    stat = nc_inq_varid(ncid, "buildx", &x_id);
    stat = nc_inq_varid(ncid, "buildy", &y_id);
    stat = nc_inq_varid(ncid, "buildz", &z_id);

    start[0] = 0;
    count[0] = nbuild_len;
    stat = nc_put_vara_double(ncid, x_id, start, count, build[bno].bx);
    stat = nc_put_vara_double(ncid, y_id, start, count, build[bno].by);
    stat = nc_put_vara_double(ncid, z_id, start, count, build[bno].db);

    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

/* Read the grid from the netcdf file fname */
void ReadBuildnetcdf(char *fname, int bno)
{
    int i, j, cnt, nel, newelement;
    int *ellist;

    int ncid;			/* netCDF id */

    /* dimension ids */
    int nbuild_dim;

    /* dimension lengths */
    size_t nbuild_len;
    size_t start[1], count[1];

    /* variable ids */
    int x_id;
    int y_id;
    int z_id;

    int stat;

    stat = nc_open(fname, NC_NOWRITE, &ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "nbuild", &nbuild_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dim(ncid, nbuild_dim, "nbuild", &nbuild_len);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "buildx", &x_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "buildy", &y_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "buildz", &z_id);
    check_err(stat, __LINE__, __FILE__);

    Free_build(bno);
    Allocate_build(bno, nbuild_len);

/*
 * Read build points
 */
    start[0] = 0;
    count[0] = nbuild_len;
    stat = nc_get_vara_double(ncid, x_id, start, count, build[bno].bx);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, y_id, start, count, build[bno].by);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, z_id, start, count, build[bno].db);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

#endif
