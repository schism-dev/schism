/* $Id: globals.h,v 1.3 2006/02/18 16:44:07 pturner Exp $
 *
 * globals.h - global variables for gredit
 *
 */

extern double hypot(double x, double y);

extern char ps_prstr[], noprint[];

extern int use_colors;
extern int revflag;
extern int doclear;
extern int noerase;
extern int overlay;
extern int bgcolor;
extern int inwin;
extern int backingstore;
extern int win_h, win_w;
extern int save_images;
extern int save_images_type;
extern int save_images_count;
extern char save_images_fname[];

#define EDIT_GRID 0
#define BACKGROUND_GRID 1
#define BACKGROUND_BOUNDARY 2
#define EDIT_BOUNDARY 3
#define EDIT_GRID_NODE_NUMBERS 4
#define EDIT_GRID_ELEMENT_NUMBERS 5
#define EDIT_GRID_FILLED 6
#define EDIT_GRID_ISOLINES 7
#define BACKGROUND_GRID_ISOLINES 8
#define BUILD_POINTS 9
#define EDIT_GRID_DEPTHS 10
#define EDIT_GRID_NODES 11
#define CONV_GRID MAXGRIDS+1

#ifdef MAIN

int debuglevel = 0;

Isolparms curip;

int swapBytes = 0;

/* projection parameters */
double proj_clat;
double proj_clon;
int proj_project = -1;
int proj_inv;

Grid grid[MAXGRIDS + 3];
int curgrid = 0;
int backgrid = MAXGRIDS;
int convgrid = MAXGRIDS + 1;
int buildgrid = MAXGRIDS + 2;

Boundary boundary[MAXBOUNDS];

Slicepts slices[MAXSLICES];

Buildpts build[MAXBUILD];
Element *tmptable;
int curbuild = 0;
double mindist = 0.0;

Boundaries *cb;			/* coastal boundaries */
int ncb;

double genfd_dx, genfd_dy;
int genfd_nx, genfd_ny;
double genfd_lx, genfd_ly;
double genfd_sx, genfd_sy;
double genfd_rot;
double genfd_x1, genfd_y1;
double genfd_x2, genfd_y2;
double genfd_x3, genfd_y3;
double genfd_x4, genfd_y4;
int genfdflag;

int curnlist = 0;
NodeList nodelist[MAXNLISTS];

int noverride = 0;
double belel_tol = -1e-04;
double qcutoff = 45.0;
int fixed_scale = 1;

double xg1 = 0.0, yg1 = 0.0, xg2 = 1.0, yg2 = 1.0;
double sxg1 = 0.0, syg1 = 0.0, sxg2 = 1.0, syg2 = 1.0;
double dx = 1.0, dy = 1.0;
double wfact = 0.3;
double redfact = 1.0;
double rdel = 0.2;
double gscale = 1.0;
/*
 * flags for display objects
 */
int display_flags[30], display_order[30], ndisplay = 11;
int display_edit_isolines = 0;
int display_back_isolines = 0;
int display_back_grid = 0;
int display_back_boundary = 0;
int display_edit_grid_boundary = 1;
int display_edit_grid = 1;
int display_node_numbers = 0;
int display_elem_numbers = 0;
int display_node_depths = 0;
int display_build_points = 1;
int display_build_line = 0;

/*
 * flags for fuzzing out items in the UI
 */
int ebound_defined = 0;
int ibound_defined = 0;
int build_defined = 0;

int ibuf[MAXSCRATCH];
int nibuf;
int ipolyx[MAXSCRATCH];
int ipolyy[MAXSCRATCH];
int nipoly;
double polyx[MAXSCRATCH];
double polyy[MAXSCRATCH];
int npoly;

double regionx[MAXSCRATCH];
double regiony[MAXSCRATCH];
int nregion;
int region_flag = 0;

int cutgrid_type = 0;

int nels = 0;

plotstr pstr[MAXSTR];		/* strings */
boxtype boxes[MAXBOXES];	/* boxes */
linetype lines[MAXLINES];	/* lines */

/* lines and boxes flags */
int box_color = 1;
int box_lines = 1;
int box_linew = 1;
int box_fill = NONE;
int box_fillpat = 1;
int box_fillcolor = 1;
int box_loctype = WORLD;

int line_color = 1;
int line_arrow = 0;
int line_lines = 1;
int line_linew = 1;
int line_loctype = WORLD;
double line_asize = 1.0;
int line_atype = 0;

/* default string parameters */
int string_color = 1;
int string_linew = 1;
int string_font = 4;
int string_rot = 0;
int string_just = 0;
int string_fill = NONE;
int string_fillpat = 1;
int string_fillcolor = 1;
int string_loctype = WORLD;
double string_size = 1.0;
int string_sym = 0;
int string_type = 0;
int string_grid = 0;
int string_prec = 0;
int string_format = 0;
double string_symsize = 1.0;

Station *sta;
int nsta;

#else				/* in other than main.c */

extern int debuglevel;

extern int swapBytes;

/* projection parameters */
extern double proj_clat;
extern double proj_clon;
extern int proj_project;
extern int proj_inv;

extern Isolparms curip;

extern Grid grid[];
extern int curgrid;
extern int backgrid;
extern int convgrid;
extern int buildgrid;

extern Boundary boundary[];

extern Buildpts build[];
extern Element *tmptable;
extern int curbuild;
extern double mindist;

extern Boundaries *cb;		/* coastal boundaries */
extern int ncb;

extern double genfd_dx, genfd_dy;
extern int genfd_nx, genfd_ny;
extern double genfd_lx, genfd_ly;
extern double genfd_sx, genfd_sy;
extern double genfd_x1, genfd_y1;
extern double genfd_x2, genfd_y2;
extern double genfd_x3, genfd_y3;
extern double genfd_x4, genfd_y4;
extern double genfd_rot;
extern int genfdflag;

extern int curnlist;
extern NodeList nodelist[];

extern int noverride;
extern double belel_tol;
extern double qcutoff;
extern int fixed_scale;

extern double *buildx, *buildy;
extern double nbuild;

extern double xg1, yg1, xg2, yg2;
extern double sxg1, syg1, sxg2, syg2;
extern double dx, dy;
extern double wfact;
extern double redfact;
extern double rdel;
extern double gscale;

extern int display_flags[], display_order[], ndisplay;
extern int display_edit_isolines;
extern int display_back_isolines;
extern int display_back_grid;
extern int display_back_boundary;
extern int display_edit_grid_boundary;
extern int display_edit_grid;
extern int display_node_numbers;
extern int display_elem_numbers;
extern int display_node_depths;
extern int display_build_points;
extern int display_build_line;

extern int ebound_defined;
extern int ibound_defined;
extern int build_defined;

extern int ibuf[];
extern int nibuf;
extern int ipolyx[];
extern int ipolyy[];
extern int nipoly;
extern double polyx[];
extern double polyy[];
extern int npoly;

extern double regionx[];
extern double regiony[];
extern int nregion;
extern int region_flag;
extern int cutgrid_type;

extern int nels;

extern plotstr pstr[];		/* strings */
extern boxtype boxes[];		/* boxes */
extern linetype lines[];	/* lines */

extern int box_color;
extern int box_lines;
extern int box_linew;
extern int box_fill;
extern int box_fillpat;
extern int box_fillcolor;
extern int box_loctype;

extern int line_color;
extern int line_arrow;
extern int line_lines;
extern int line_linew;
extern int line_loctype;
extern double line_asize;
extern int line_atype;

extern int string_color;
extern int string_linew;
extern int string_font;
extern int string_rot;
extern int string_just;
extern int string_fill;
extern int string_fillpat;
extern int string_fillcolor;
extern int string_loctype;
extern double string_size;
extern int string_sym;
extern int string_type;
extern int string_grid;
extern int string_prec;
extern int string_format;
extern double string_symsize;

extern double charsize, xlibcharsize;	/* declared in draw.c and xlib.c
					 * resp. */
extern Station *sta;
extern int nsta;

#endif				/* ifdef MAIN */
