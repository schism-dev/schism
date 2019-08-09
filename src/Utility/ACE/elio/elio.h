/*! @file $RCSfile: elio.h,v $
 * 
 * @version $Revision: 1.34 $
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 * Header files for I/O routines for new file format
 * pturner 11-2002
 *
 */
 
#ifndef ELIO_H
#define ELIO_H

/*! Error codes. */
/*! No errors */
#define ELIO_OK 0 
/*! Generic error return code */
#define ELIO_ERR 1 
/*! Error return code for fopen() */
#define ELIO_FOPEN_ERR 2 
/*! Error return code for fseek() */
#define ELIO_FSEEK_ERR 3 
/*! Error return code for fread() */
#define ELIO_FREAD_ERR 4 
/*! Error return code for fwrite() */
#define ELIO_FWRITE_ERR 5 
/*! Error return code for fclose() */
#define ELIO_FCLOSE_ERR 6 
/*! Error return code for ftell() */
#define ELIO_FTELL_ERR 7 
/*! Compression types, COMPRESS_NONE == no compression */
#define COMPRESS_NONE 0 
/*! COMPRESS_C16 == Fixed point 16 bit integer compression. */
#define COMPRESS_C16 1 
/*! ZLEVEL == 0 */
#define ZLEVEL 0 
/*! SIGMA_S0 == sigma grid scheme 0 */
#define SIGMA_S0 1 
/*! SIGMA_S1 == sigma grid scheme 1 */
#define SIGMA_S1 2 
/*! SIGMA_S2 == sigma grid scheme 2 */
#define SIGMA_S2 3 
/*! SIGMA_S3 == sigma grid scheme 3 */
#define SIGMA_S3 4 
/*! SIGMA_S4 == sigma grid scheme 4 */
#define SIGMA_S4 5 
/*! Missing data value for compressed files */
#define MISSING_DATA -32767

extern int ElioErr; /*! Elcirc error ID */
extern int ElioSysErr; /*! system error number */
extern int ElioLineErr; /*! Line at which error occurred in elio.c */
extern char *ElioMsg; /*! Error message */

/*! ElcircHeader:
 * A struct used to hold the header information from the data file. 
 * Some elements of ElcircHeader are derived rather than actually being 
 * present in the file. An example of a derived variable is the member 
 * 'no' which provides 
 * offsets into the file for a given time step to access nodes directly.
 * The elements of the 'bi' member used to hold bottom indices are one less 
 * than the on disk header as a convenience for C coders.
 */
typedef struct {
/*! magic string */
    char magic[49]; 
/*! File format version ID */
    char version[49]; 
/*! Start time of the run DD/MM/YYYY HH:MM:SS TZ */
    char start_time[49]; 
/*! Variable name for data in this file */
    char variable_nm[49]; 
/*! Variable dimension */
    char variable_dim[49]; 
/*! version number == 2 or == 3 or == 4, derived from 'version' */
    int v; 
/*! Compression type 0 == no compression */
    int compress; 
/*! If 16bit fixed point compression, the slope of the line */
    double a; 
/*! If 16bit fixed point compression, the intercept of the line */
    double b; 
/*! Number of time steps in the file, NOTE: this may not be the actual number of time steps in the file */
    int nsteps; 
/*! Model time step */
    float timestep; 
/*! Number of time steps skipped */
    int skip; 
/*! Vector or scalar data, values are 1 = scalar or 2 = vector */
    int ivs; 
/*! i23d = 2 => 2d (such as elevations or wind) or i23d = 3 => 3d data (such as velocities and salinity) */
    int i23d; 
/*! 1.0 or 0.5 depending whether variable is defined on the level or half level */
    float vpos; 
/*! correction to mean sea level */
    float zmsl; 
/*! Version 4 sigma */
    int ivcor;
/*! Version 4 sigma */
    float h0;
/*! Version 5 sigma */
    float hs;
/*! Version 4 sigma */
    float hc;
/*! Version 4 sigma */
    float thetab;
/*! Version 4 sigma */
    float thetaf;
/*! Number of vertical levels. Note: this is set to 1 for 2d variables where the on-disk header will have the number of levels as in a 3d variable file */
    int nvrt; 
/*! Number of vertical levels for Z - version 5.0 */
    int kz; 
/*! Number of vertical levels for sigma - version 5.0 */
    int ks; 
/*! number of elements in a time step */
    int nitems; 
/*! Size of the header on disk */
    int hsize; 
/*! Time step size */
    int ssize; 
/*! Sigma grid if true */
    int sigma; 
/*! depth of each level or sigma coords */
    float *zcor; 
/*! Number of nodes in the grid */
    int np; 
/*! Number of elements in the grid */
    int ne; 
/*! X, Y for nodes, eventually doubles in some future version. */
    float *x; 
    float *y; 
/*! The grid depth */
    float *d; 
/*! Holds the indexes of the bottom level for each node starting from 0. The on-disk header stores these indexes starting from 1. */ 
    int *bi; 
/*! Save the original bottom index as 0 can indicate dry - this is a special purpose thing done for MPI output */ 
    int *bisave; 
/*! Node offsets into time step for each node starting from 0 */
    int *no; 
/*! Array holding the element type for each element. Each element is assigned either a 3 (tris) or 4 (quads) */
    int *etype; 
/*! Pointers to table of elements for the grid. */
    int *icon[4]; 
} ElcircHeader;

typedef struct {
/*! Model time */
    float t; 
/*! Model time step number */
    int it; 
/*! Surface indices, these are not transformed but are as read from the file so the maximum value could be nvrt (as opposed to the bottom indices where are decremented by one after being read in). */
    int *surfind; 
/*! In version 4, the elevation for each time step. */
    float *e; 
/*! Data for all nodes and levels for the time step. */
    float *d; 
} ElcircTimeStep;

typedef struct {
    float t;
    int it;
    int istart;
    int istop;
    int npts;
    float *d;
} ElcircTimeStepAtNode;

/*! ElioGrid:
 *
 * Structure for a stand alone grid. Useful where no header exists.
 */
typedef struct {
/*! Number of elements */
    int ne;
/*! Number of nodes */
    int np;
/*! X */
    double *x;
/*! Y */
    double *y;
/*! Depth */
    double *d;
/*! Element type (3 or 4) */
    int *etype;
/*! Table of elements */
    int *icon[4];
/* number of open boundaries */
    int nopen;
    int ntotalopen;
    int *oseg;
/* number of land boundaries */
    int nland;
    int ntotalland;
    int *lseg;
/*! boundary */
    int **b;
} ElioGrid;

ElioGrid *ElioAllocateGrid(int np, int ne);
int ElioReadGrid(char *fname, ElioGrid *g);
int ElioFindElementInGrid(ElioGrid *g, double x, double y);
int ElioGetCoefficientsGrid(ElioGrid *g, int elem, double xp, double yp, double *w);
int ElioGridFindNearestNode(ElioGrid * g, double x, double y);
int ElioFindNearestNode(ElcircHeader * h, double x, double y);
int ElioGetGridElementCenter(ElioGrid *g, int elem, double *cx, double *cy);
double ElioGetGridElementArea(ElioGrid *g, int elem);

/*!
 *  Functions to get the header, print the header to stdout, free memory associated with a header.
 */

/*!
 * Read the header from an ELCIRC output file. Returns ELIO_OK on success and ELIO_ error codes otherwise.
 */  
int ElioGetHeader(char *fname, ElcircHeader *h);
/*! Allocate memory allocated for a header */
int ElioAllocateHeader(ElcircHeader *h);
/*! Free memory allocated for a header */
void ElioFreeHeader(ElcircHeader *h);
/*! write a summary of the header to stdout */
void ElioPrintHeader(ElcircHeader *h);

/*!
 * Functions to get information about the data and to read data as chunks or from nodes.
 */
int ElioGetNStepsInFile(char *fname, ElcircHeader *h);
int ElioGetTimeStep(FILE *fp, int step, ElcircHeader *h, ElcircTimeStep *t);
void ElioFreeTimeStep(ElcircTimeStep *t);
int ElioAllocateTimeStep(ElcircHeader *h, ElcircTimeStep *t);
int ElioGetNode(FILE *fp, int step, int node, ElcircHeader *h, float *t, int *it, int *bind, int *sind, float *d);
int ElioGetNodeOld(FILE *fp, int step, int node, ElcircHeader *h, float *t, int *it, float *d);
int ElioExtractNode(ElcircTimeStep *t, ElcircHeader *h, int node, int bind, int sind, float *d);
int ElioGetXYData(FILE * fp, int step, double x, double y, ElcircHeader *h, float *t, int *it, int *bind, int *sind, float *d);
int ElioGetXYData2(FILE * fp, int step, int elem,  ElcircHeader *h, double *hh, float *t, int *it, int *bind, int *sind, float *d);
int ElioInterpTimeStep(ElcircHeader *h, int elem, double x, double y, double *hh, ElcircTimeStep *t, int *bind, int *sind, float *d);

int ElioGetPoint(FILE *fp, int step, int node, int level, ElcircHeader *h, float *t, int *it, float *d);
/*! create a scalar data set compatible with the old format */
int ElioMakeScalarsOld(ElcircHeader *h, float *d, float *dd);
/*! create a vector data set compatible with the old format */
int ElioMakeVectorsOld(ElcircHeader *h, float *d, float *dd);

/*! extract a grid from the header and make a new header */
int ElioExtractGrid(ElcircHeader *h1, int n, double *xb, double *yb, int *isin, int invert, int most, ElcircHeader *h2);
/*! extract data from a timestamp for a sampled data set */
int ElioExtractData(ElcircHeader *h1, ElcircHeader *h2, int *isin, ElcircTimeStep t1, ElcircTimeStep * t2);
int ElioIntersectToLeft(double x, double y, double x1, double y1, double x2, double y2);
int ElioInPolygon(double x, double y, int n, double *xlist, double *ylist);
int ElioFindElementXY(ElcircHeader *h, double x, double y);
int ElioInsideElement(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3);
int ElioInsideElement4(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

/*! Interpolators */
int ElioEvalFlowXY(int n, double *h, double *u, double *v, double *uret, double *vret);
int ElioEvalScalarXY(int , double *, double *u, double *uret);
int ElioEval(int , double *, double *u, double *uret);
int ElioGetCoefficients(ElcircHeader *h, int elem, double xp, double yp, double *w);
int ibilinear(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double x, double y, double *xi, double *eta, double *shapef);
/*! Get the depths for an interpolated point. */
int ElioGetZPos(ElcircHeader *h, double d, double e, int *ilev1, int *ilev2, double *zpos);

/*!
 * Extract or compute special data sets from a time step
 */
/* Surface values */
int ElioGetSurfaceHeader(ElcircHeader *h1, ElcircHeader *h2);
int ElioGetSurfaceStep(ElcircHeader *h1, ElcircHeader *h2, ElcircTimeStep t1, ElcircTimeStep *t2);
/* Bottom values */
int ElioGetBottomHeader(ElcircHeader *h1, ElcircHeader *h2);
int ElioGetBottomStep(ElcircHeader *h1, ElcircHeader *h2, ElcircTimeStep t1, ElcircTimeStep *t2);
/* Transect values */
int ElioGetTransectHeader(ElcircHeader *h1, int n, double *x, double *y, ElcircHeader *h2);
int ElioGetTransectStep(ElcircHeader *h1, ElcircHeader *h2, int step, int n, double *x, double *y, ElcircTimeStep t1, ElcircTimeStep *t2);
/* A single level */
int ElioGetLevelHeader(ElcircHeader *h1, ElcircHeader *h2);
int ElioGetLevelStep(ElcircHeader *h1, ElcircHeader *h2, ElcircTimeStep t1, int level, ElcircTimeStep *t2);
/* A level at arbitrary depth */
int ElioGetZLevelHeader(ElcircHeader *h1, ElcircHeader *h2);
int ElioGetZLevelStep(ElcircHeader *h1, ElcircHeader *h2, ElcircTimeStep t1, float z, ElcircTimeStep *t2);

/*!
 * Functions to write to a file data sets in the new format
 */
/* Write a header to a file */
int ElioPutHeader(FILE *fp, ElcircHeader *h);
int ElioPutTimeStep(FILE *fp, ElcircHeader *h, ElcircTimeStep t);
int ElioPutHeaderOld(FILE *fp, ElcircHeader *h);
int ElioPutTimeStepOld(FILE *fp, ElcircHeader *h, ElcircTimeStep t);
/* Write to a single step .61 file */
int ElioPut61File(char *fn, ElcircHeader * h, double *v);
int ElioGetFileType(FILE *fp);

/*!
 * Misc. functions.
 */
double ElioGetElementArea(ElcircHeader * h, int elem);
int ElioMinMax(int n, double *ar, double *dmin, double *dmax);
int ElioIntMin(int n, int *iar);
int ElioIntMax(int n, int *iar);
/* For interpolation */
int ElioFindIndex(int n, double *data, double d);
void ElioInterpolateArray(int n, double *z, double *s, int np, double *zp, double *sp);
double ElioInterpolate(int n, double *z, double *s, double zp);
double ElioInterpolateAtIndex(int n, double *z, double *s, int ind, double zp);

/*!
 * Time functions
 */
double ElioGetDay(int mon, int day, int year, int h, int mi, double se);
void ElioGetYearDay(double jd, int *y, int *yd, int *h, int *mi, double *sec);
void ElioSetCorieTime(int b);

#endif /* ELIO_H */
