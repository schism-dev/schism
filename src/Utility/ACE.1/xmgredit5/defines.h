/* $Id: defines.h,v 1.4 2006/02/18 16:44:05 pturner Exp $
 *
 * defines and typedefs
 */

#define MAXGRIDS 5
#define MAXBOUNDS 20000      /*max. # of bnd nodes on each segment*/
#define MAXGRIDBOUNDS 10000 /*max. # of bnd segments*/
#define MAXBUILD 5
#define MAXSLICES 15
#define MAXSCRATCH 5000
#define MAXNLISTS 150

#define MAXGRAPH 10		/* max number of graphs */
#define MAXPLOT 30		/* max number of XY data sets */
#define MAXFIT 15		/* max number degree of polynomial fit */
#define MAX_ZOOM_STACK 20	/* max stack depth for world stack */
#define MAXREGION 5		/* max number of regions */
#define MAX_TICK_LABELS 35	/* max number of user defined ticks/labels */
#define MAX_STRING_LENGTH 256	/* max length for strings */
#define MAX_LINESTYLE 5		/* max number of linestyles */
#define MAX_LINEWIDTH 10	/* max number of line widths */
#define MAX_SET_COLS 6		/* max number of columns per XY set */
#define MAXZOOMBOXES 10		/* max number of zoom boxes */


#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef M_PI
#define M_PI       3.1415926535897931160E0
#endif

#ifndef M_PI_2
#define M_PI_2     1.5707963267948965580E0
#endif

#ifndef M_2_PI
#define M_2_PI     6.3661977236758134308E-1
#endif

#define HORIZONTAL 1
#define VERTICAL 2

#define MAX_STRING_LENGTH 256	/* max length for strings */
#define MAXBOXES 150		/* max number of boxes */
#define MAXLINES 150		/* max number of lines */
#define MAXSTR 200		/* max number of strings */
#define MAXSTRLEN 140
#define MAXISOLINES 20		/* max number of isolines */

#define SELECT_REGION 1
#define RUBBER_BAND 2
#define CERC_RENUM 3
#define CERC_RENUM1 4
#define CERC_RENUM2 5
#define CERC_WRITE 6
#define DEFINE_REGION 7
#define SELECT_SNAP 8
#define SNAP_LINE1ST 9
#define SNAP_LINE2ND 10
#define DEF_EXT_BOUND 11
#define DEF_INT_BOUND 12
#define MOVE_EXT_BOUND1ST 13
#define MOVE_EXT_BOUND2ND 14
#define MOVE_INT_BOUND1ST 15
#define MOVE_INT_BOUND2ND 16
#define DEL_EXT_BOUNDPT 17
#define DEL_INT_BOUNDPT 18
#define DEL_INT_BOUND 19
#define DEL_BOUND_PT 20
#define ADD_BOUND1 21
#define ADD_BOUND2 22
#define ADD_BOUND3 23
#define ADD_NODE_BOUND1 24
#define ADD_NODE_BOUND2 25
#define ADD_NODE_BOUND3 26
#define PLACE_BUILD 27
#define PLACE_BUILD_ARC 28
#define DELETE_BUILD 29
#define MOVE_BUILD1ST 30
#define MOVE_BUILD2ND 31
#define ADD_BUILD1 32
#define ADD_BUILD2 33
#define ADD_BUILD3 34
#define ADD_NODE 35
#define GET_NEAREST_NODE 36
#define GET_NEAREST_ELEMENT 37
#define GET_ELEMENT 38
#define GET_GRID_DEPTH 39
#define GET_BACK_DEPTH 40
#define GET_ALL 41
#define MOVE_NODE1ST 42
#define MOVE_NODE2ND 43
#define SWAPLINE_1ST 44
#define SWAPLINE_2ND 45
#define CUT_GRID1ST 46
#define CUT_GRID2ND 47
#define ADD_ELEMENT1 48
#define ADD_ELEMENT2 49
#define ADD_ELEMENT3 50
#define DELETE_ELEMENT 51
#define DELETE_ELEMENTS 52
#define SPLIT_ELEMENT3 53
#define SPLIT_ELEMENT4 54
#define SPLIT_ELEMENTS3 55
#define SPLIT_ELEMENTS4 56
#define ZOOM_1ST 57
#define ZOOM_2ND 58
#define COMP_AREA 59
#define COMP_PERIMETER 60
#define SEL_POINT 61
#define SLICE_BATH1ST 62
#define SLICE_BATH2ND 63
#define PICK_BATH 64
#define SLICE_GRIDBATH1ST 65
#define SLICE_GRIDBATH2ND 66
#define SLICE_GRIDBATH3RD 67
#define PICK_ISTOK 68
#define PICK_PROP_NODE 69
#define PICK_PROP_ELEM 70
#define PLACE_MAPSCALE 71
#define PLACE_GRIDBATH_LEGEND 72
#define PLACE_BACKBATH_LEGEND 73
#define PLACE_CONC_LEGEND 74
#define STR_LOC 75
#define MOVE_OBJECT_1ST 76
#define MOVE_OBJECT_2ND 77
#define DEL_OBJECT 78
#define MAKE_BOX_1ST 79
#define MAKE_BOX_2ND 80
#define EDIT_OBJECT 81
#define STR_LOC1ST 82
#define STR_LOC2ND 83
#define MAKE_LINE_1ST 84
#define MAKE_LINE_2ND 85
#define STR_EDIT 86
#define PLACE_VSCALE 87
#define PLACE_GRAD_SCALE 88
#define PLACE_SLOPE_LEGEND 89
#define FIND_GRAD 90
#define FIND_SLOPE 91
#define OPENB1ST 111
#define OPENB2ND 112
#define LANDB1ST 113
#define LANDB2ND 114
#define EXTRACT_GRID 118
#define ADD_QUAD1 137
#define ADD_QUAD2 138
#define ADD_QUAD3 139
#define ADD_QUAD4 140
#define CONV_4TO3 142
#define CONV_3TO4_1ST 143
#define CONV_3TO4_2ND 144
#define EDIT_NODE 148
#define CONV_NODE2COL 149
#define ADD_NODELIST 154
#define ADD_NODELIST1ST 155
#define ADD_NODELIST2ND 156
#define DELETE_NODE 157
#define ADD_BCNODELIST 158
#define ADD_BCNODELIST1ST 159
#define ADD_BCNODELIST2ND 160
#define PICK_CHAIN 163
#define PICK_LIST 164
#define PICK_ELEM_LIST 165
#define PICK_NODE_LIST 166
#define GENFD_1ST 175
#define GENFD_2ND 176
#define GENFD_3RD 177
#define PLACE_ISOLINE_LABEL 178
#define GET_NEAREST_BUILDPT 179
#define GET_DEPTH_ALL 180
#define STRIKE_1ST 181
#define STRIKE_2ND 182
#define PLACE_EPICENTER 183
#define PLACE_ISOLINES_LEGEND 188
#define SELECT_BCNODE 199
#define VIEW_1ST 204
#define VIEW_2ND 205
#define QUERY_PROP_NODE 208
#define QUERY_PROP_ELEM 209

#define NODE_LIST 1
#define ELEM_LIST 2
#define CHAIN_LIST 3

#define MIDDLE 1
#define CORNER 2
#define COLOCATED 3

#define STRING 0
#define BOX  1
#define LINE 2

#define NONE 0
#define COLOR 1
#define PATTERN 2
#define WORLD 3
#define VIEW 4
#define ON 5
#define OFF 6
#define DECIMAL 7
#define EXPONENTIAL 8
#define GENERAL 9

/* graphics output to the following */
#define GR_X            0
#define GR_PS_L         1	/* PostScript landscape */
#define GR_PS_P         2	/* PostScript portrait */
#define GR_XWD          9	/* X window dump */

/* set HDEV to the default hardcopy device */
#ifndef HDEV
#       define HDEV GR_PS_L
#endif

/* TDEV is always X */
#define TDEV GR_X

#define FREE 0
#define LANDSCAPE 1
#define PORTRAIT 2

/* Element types */
#define LINEAR1D 2
#define QUAD1D 1
#define COLOC 0
#define TRANSITION 5
#define LINEAR2D 3
#define QUAD2D 6
#define QUADRANGLE4 4
#define QUADRANGLE8 8

typedef struct _Props {
    int color;
    int linew;
    int lines;
    int font;
    int format;
    int prec;
    int points;
    double charsize;
    int symbol;
    int symsize;
    int fill;
    int fillusing;
    int fillcol;
    int fillpat;
    int arrow;
    int atype;
    double asize;
} Props;

typedef struct isolparms {
    int active;
    int type;
    Props p;
    int nisol;
    int nsplits;
    int marker;
    int markstep;
    int numbs;
    int loctype;
    double x;
    double y;
    int xlen;			/* in character widths */
    int ylen;			/* in character widths */
    int xgap;			/* in character widths */
    int ygap;			/* in character widths */
    int lactive;		/* legend active */
    int layout;			/* legend HORIZONTAL or VERTICAL */
    int llabels;		/* legend labels on or off */
    int isoltype;
    double cmin;
    double cmax;
    double cint;
    double cis[MAXISOLINES];
    int color[MAXISOLINES];
    int linew[MAXISOLINES];
    int lines[MAXISOLINES];
} Isolparms;

/* nodes defining a linear element */
typedef struct {
    int type;			/* general type */
    int nn;			/* number of nodes */
    int wetdry;			/* for finite differece grid */
    int ngeom;			/* number of nodes describing geometry */
    int nl[6];			/* node numbers - geometry first */
} Element;

/* nodal properties */
typedef struct {
    int type;			/* general type */
    int utype;			/* user type */
    int col;			/* collocated node element */
    int n;			/* number of points for x-section */
    double *x, *y;		/* locations of points for x-section */
} NodeProps;

/* connectivities */
typedef struct {
    int cnt;
    int *el;
} conlist;

/* boundary definition */
typedef struct {
    int nbpts;			/* number of points in this boundary */
    int boundtype;		/* 0 if external 1 if internal */
    int bactive;		/* boundary is active */
    double *boundx;		/* boundary x */
    double *boundy;		/* boundary y */
    int *nodes;			/* boundary as node numbers */
    int *btype;			/* boundary node type */
} Boundary;

typedef struct _Boundaries {
    int type;
    int n;
    int display;
    int color;
    int linew;
    int lines;
    int sym;
    char name[256];
    Boundary *b;
} Boundaries;

/* build points definition */
typedef struct {
    int nbuild;			/* number of build points */
    int buildactive;		/* build points active */
    int *fromdb;		/* point from database */
    double *bx;			/* build x */
    double *by;			/* build y */
    double *db;			/* depth at build points */
    int color;
    int sym;
    double size;
    int linew;
    int isol;
} Buildpts;

/* list of nodes */
typedef struct {
    int type;
    int active;
    int gridno;
    int n;
    int *nodes;
    int display;
    int color;
    int linew;
    char name[64];
} NodeList;

/* slice points definition */
typedef struct {
    int nslice;			/* number of slice points */
    int sliceactive;		/* slice points active */
    double *sx;			/* slice x */
    double *sy;			/* slice y */
    double *d;			/* depth at slice points */
    int color;
    int linew;
} Slicepts;

/* grid definition */
typedef struct {
    int type;
    char name[256];
    char fname[256];
    char alphid[256];		/* study name */
    int nmel;			/* number of elements */
    int nmnp;			/* number of nodes */
    int gridtype;		/* 0 undefined, 1 linear, 2 quadratic */
    int readonly;		/* 0 editable, 1 not */
    int transform;		/* 0 not known, 1 meters, 2, lon/lat */
    int depthflag;		/* depths at nodes or elements */
    double *xord;		/* node x */
    double *yord;		/* node y */
    double *depth;		/* node depth */
    double *edepth;		/* element depth */
    Element *icon;		/* table of elements */
    conlist *nodecon;		/* elements associated with a given node */
    int *ntype;
    int *ellist;
    int *nlist;
    int nnodep;
    NodeProps *np;
    int nbounds;		/* number of boundaries */
    int boundaries[MAXGRIDBOUNDS];	/* boundaries - 0 =  external
					 * boundary */
    double xmin, xmax, ymin, ymax, dmin, dmax;
    double scalef;
    int gactive;		/* grid is active */
    int bactive;		/* boundary is active */
    int gcolor;
    int glines;
    int glinew;
    int bcolor;
    int blines;
    int blinew;
    int ibcolor;
    int iblines;
    int iblinew;
    Isolparms ip;
    char projection[128];
    int proj;
    double lon;
    double lat;
    int invproj;
    int fdgrid; /* true if grid is finite difference */
    int nx;
    int ny;
} Grid;

/*
 * typedefs for objects
 */
typedef struct {
    int active;
    int loctype;
    int gno;
    double x1, y1, x2, y2;
    int lines;
    int linew;
    int color;
    int fill;
    int fillcolor;
    int fillpattern;
} boxtype;

typedef struct {
    int active;
    int loctype;
    int gno;
    double x1, y1, x2, y2;
    int lines;
    int linew;
    int color;
    int arrow;
    int atype;
    double asize;
} linetype;

typedef struct {
    int active;
    int type;
    int loctype;
    int gno;
    int el;
    double x, y, d;
    int lines;
    int linew;
    int color;
    int rot;
    int font;
    int just;
    double charsize;
    int sym;
    double symsize;
    char s[MAXSTRLEN + 1];
} plotstr;

typedef struct {
    int display;		/* tide station is displayed and how */
    int display_marker;		/* marker is displayed */
    int type;			/* fixed, strip chart, elev marker */
    int attach;
    Props p;			/* box props */
    Props rp;			/* point props */
    int loctype;
    double x;
    double y;
    double locx;
    double locy;
    double delt;		/* if strip chart history box */
    double wx1;
    double wy1;			/* world coordinates */
    double wx2;
    double wy2;			/* world coordinates */
    double vx;
    double vy;			/* viewport delta */
    double xtickm;
    double ytickm;
    int nsteps;
    double start;
    double step;
    double stop;
    int precx;
    int precy;
} DisplayTideStation;

typedef struct {
    int active;			/* on or off */
    int type;			/* grid or at stations */
    char name[256];		/* file name */
    char sname[256];		/* station name */
    double dcor;		/* depth correction */
    int loctype;		/* 0 = x, y, 1 = attach to grid */
    char gname[256];		/* grid name to attach to for loctype */
    int grid;			/* number of grid */
    int node;			/* use nodes or x, y */
    double x;			/* location */
    double y;			/* for stations */
    int nfreq;			/* total number of frequencies */
    char **freqname;		/* frequency names */
    double *omega;		/* frequencies in rad/sec */
    double *elamp;		/* amplitudes of elevation */
    double *elphase;		/* phases */
    double *ampx;		/* amplitudes of velocity */
    double *ampy;		/* amplitudes of velocity */
    double *phax;		/* phases of velocity */
    double *phay;		/* phases of velocity */
    double emin;		/* extrema */
    double emax;
    double umin;
    double umax;
    double vmin;
    double vmax;
    DisplayTideStation t;
} TideStation;

typedef struct _TidalConstituent {
    int index;
    double freq;
    char *name;
    double period;
} TidalConstituent;

typedef struct {		/* a station */
    int active;
    int type;
    char fname[256];		/* file name */
    char label[256];		/* label for marker */
    double x;
    double y;			/* location and hotspot */
    int display;		/* displayed and how */
    Props p;
    Isolparms ip;		/* for isolines */
} Station;

typedef struct {
    int display;		/* grid is displayed and how */
    int display_bath;		/* grid bathymetry is displayed and how */
    int display_boundary;	/* grid boundary is displayed and how */
    int display_nodes;		/* display node numbers */
    int display_elements;	/* display element numbers */
    int display_depths;		/* display depths */
    int display_gridf;		/* display grid filled */
    int display_courant;	/* display courant numbers */
    int display_dimw;		/* display dimensionless wavelength */
    int display_flags[50];	/* which items to display - for the side bar */
    int step;			/* step to display */
    double redfact;		/* reduction factor for elements */
    double widthfact;		/* for RITA grids */
    Props p;			/* grid props */
    Props bp;			/* boundary props */
    Isolparms ip;		/* for isolines of bathymetry */
    Isolparms cip;		/* for courant numbers */
    Isolparms dip;		/* for dimensionless wl */
} DisplayGrid;

typedef struct {
    int display;		/* boundary is displayed and how */
    int close;			/* connect first and last point */
    Props p;			/* properties */
} DisplayBoundary;

typedef struct {
    int active;
    int loctype;
    int gno;
    double x;
    double y;
    int lines;
    int linew;
    int color;
    int rot;
    int font;
    int just;
    double charsize;
    char s[MAXSTRLEN + 1];
} plotsym;

typedef struct {
    int active;
    int type;
    Props p;
    int loctype;
    double x;
    double y;
    double scale;
    int units;
    double unitfac;
    double len;
    double len2;
    plotstr s;			/* legend string for velocity legend */
} velocity_scale;

typedef struct {
    int active;
    int type;
    Props p;
    int loctype;
    double x;
    double y;
    double scale;
    int units;
    double unitfac;
    double len;
    double len2;
    plotstr s;			/* legend string for velocity legend */
} flux_scale;

typedef struct {
    int active;
    int type;
    Props p;
    int loctype;
    double x;
    double y;
    double scale;
    int units;
    double unitfac;
    double len;
} map_scale;

typedef struct {
    int active;
    int type;
    Props p;
    int loctype;
    double x;
    double y;
    double angle;
} north_indicator;

typedef struct {
    int active;
    int type;
    int display;		/* display region */
    int display_marker;		/* display marker */
    Props rp;			/* props for the region */
    Props p;			/* props for the zoom box */
    int loctype;
    double x;
    double y;
    double locx;
    double locy;
    int attach;			/* which corner to attach the connecting line */
    double expand;		/* how many times to zoom */
    double wx1;
    double wy1;			/* world coordinates */
    double wx2;
    double wy2;			/* world coordinates */
    double vx;
    double vy;			/* viewport delta */
    double xtickm;
    double ytickm;
    int precx;
    int precy;
} Zoom_box;

typedef struct {
    int active;
    int type;			/* line or polyline */
    int npts;			/* number of points */
    double x1;
    double y1;
    double x2;
    double y2;			/* if by line */
    double *sx;
    double *sy;			/* if by polyline */
    int display;		/* ON or other */
    int display_marker;		/* ON or other */
    int display_slice;		/* ON or other */
    Props p;
    int loctype;
    double x;
    double y;
    double locx;
    double locy;
    int attach;
    double wx1;
    double wy1;			/* world coordinates */
    double wx2;
    double wy2;			/* world coordinates */
    double vx;
    double vy;			/* viewport delta */
    double xtickm;
    double ytickm;
    int precx;
    int precy;
    int bath[MAXGRIDS];
    Props bp[MAXGRIDS];
} DisplaySlice;

typedef struct {
    char *s;
    int type;
} symtab_entry;

typedef struct {
    double xg1;
    double xg2;
    double yg1;
    double yg2;			/* window into world coords */
} world;

typedef struct {
    double xv1;
    double xv2;
    double yv1;
    double yv2;			/* device viewport */
} view;

typedef struct {
    world w;			/* current world */
    world t[3];			/* current tick spacing */
} world_stack;

typedef struct {
    plotstr title;		/* graph title */
    plotstr stitle;		/* graph subtitle */
} labels;

typedef struct {
    int axis;			/* which axis */
    int active;			/* active or not */
    int alt;			/* alternate map if TRUE */
    double tmin;
    double tmax;		/* mapping for alternate tickmarks */
    double tmajor;
    double tminor;		/* major, minor tick divisions */
    double offsx;
    double offsy;		/* offset of axes in viewport coords */
    plotstr label;		/* graph axis label */
    int label_layout;		/* axis label orientation (h or v) */
    int label_place;		/* axis label placement (specfied or auto) */
    int tl_flag;		/* toggle tickmark labels on or off */
    int mtl_flag;		/* toggle minor tickmark labels on or off */
    int tl_type;		/* either auto or specified (below) */
    int tl_layout;		/* horizontal, vertical, or specified */
    int tl_angle;		/* angle to draw labels if layout is
				 * specified */
    int tl_sign;		/* tick labels normal, absolute value, or
				 * negate */
    int tl_just;		/* justification of ticklabel and type of
				 * anchor point */
    int tl_prec;		/* places to right of decimal point */
    int tl_format;		/* decimal or exponential ticmark labels .. */
    int tl_skip;		/* tick labels to skip */
    int tl_staggered;		/* tick labels staggered */
    int tl_starttype;		/* start at graphmin or use tl_start/stop */
    int tl_stoptype;		/* start at graphmax or use tl_start/stop */
    double tl_start;		/* value of x to begin tick labels and major
				 * ticks */
    double tl_stop;		/* value of x to begin tick labels and major
				 * ticks */
    int tl_op;			/* tick labels on opposite side or both */
    double tl_vgap;		/* tick label to tickmark distance vertically */
    double tl_hgap;		/* tick label to tickmark distance
				 * horizontally */
    int tl_font;		/* font to use for labels */
    double tl_charsize;		/* character size for labels */
    int tl_color;		/* color */
    int tl_linew;		/* line width for labels */
    char tl_appstr[256];	/* append string to tick label */
    char tl_prestr[256];	/* prepend string to tick label */
    int t_type;			/* type of tickmarks, usual, xticstart, or
				 * specified */
    int t_flag;			/* toggle tickmark display */
    int t_mflag;		/* toggle minor tickmark display */
    int t_integer;		/* major tic marks on integer divisions */
    int t_num;			/* approximate default number of X-axis ticks */
    int t_inout;		/* ticks inward, outward or both */
    int t_log;			/* logarithmic ticmarks */
    int t_op;			/* ticks on opposite side */
    int t_color;		/* colors and linestyles */
    int t_lines;
    int t_linew;
    int t_mcolor;
    int t_mlines;
    int t_mlinew;		/* minor grid colors and linestyles */
    double t_size;		/* length of tickmarks */
    double t_msize;		/* length of minor tickmarks */
    int t_drawbar;		/* draw a bar connecting tick marks */
    int t_drawbarcolor;		/* color of bar */
    int t_drawbarlines;		/* linestyle of bar */
    int t_drawbarlinew;		/* line width of bar */
    int t_gridflag;		/* grid lines at major tick marks */
    int t_mgridflag;		/* grid lines at minor tick marks */
    int t_spec;			/* number of ticks at specified locations */
    double t_specloc[MAX_TICK_LABELS];
    plotstr t_speclab[MAX_TICK_LABELS];
    int spec_font;
    double spec_charsize;
    int spec_color;
    int spec_linew;
} tickmarks;

typedef struct {
    int active;			/* region on or off */
    int type;			/* region type */
    int color;			/* region color */
    int lines;			/* region linestyle */
    int linew;			/* region line width */
    int linkto[MAXGRAPH];	/* associated with graphs in linkto */
    int n;			/* number of points if type is POLY */
    double *x;
    double *y;			/* coordinates if type is POLY */
    double x1;
    double y1;
    double x2;
    double y2;			/* starting and ending points if type is not
				 * POLY */
} region;

typedef struct {
    int active;			/* frame on or off */
    int type;			/* frame type */
    int color;			/* frame color */
    int lines;			/* frame linestyle */
    int linew;			/* frame line width */
    int fillbg;			/* fill background */
    int bgcolor;		/* background color inside frame */
} framep;

typedef struct {
    int active;			/* active flag */
    int type;			/* dataset type */
    int deact;			/* deactivated set */
    double missing;		/* value for missing data */
    double *ex[MAX_SET_COLS];	/* x, y, dx, z, r, hi depending on dataset
				 * type */
    char **s;			/* pointer to strings */
    int font;			/* font for strings */
    int format;			/* format for drawing values */
    int prec;			/* precision for drawing values */
    double xmin;
    double xmax;		/* min max for x */
    double ymin;
    double ymax;		/* min max for y */
    int len;			/* set length */
    int sym;			/* set plot symbol */
    char symchar;		/* character for symbol */
    int symskip;		/* How many symbols to skip */
    int symfill;		/* Symbol fill type */
    int symdot;			/* Symbol dot in center */
    double symsize;		/* size of symbols */
    int lines;			/* set line style */
    int linew;			/* line width */
    int color;			/* color for symbol and linestyle */
    int fill;			/* fill type */
    int fillusing;		/* fill using color or pattern */
    int fillcolor;		/* fill color */
    int fillpattern;		/* fill pattern */
    int errbar;			/* if type is _DX, _DY, _DXDY and errbar =
				 * TRUE */
    int errbarxy;		/* type of error bar */
    int errbar_linew;		/* error bar line width */
    int errbar_lines;		/* error bar line style */
    int errbar_riser;		/* connecting line between error limits */
    int errbar_riser_linew;	/* connecting line between error limits line
				 * width */
    int errbar_riser_lines;	/* connecting line between error limits line
				 * style */
    double errbarper;		/* length of error bar */
    double hilowper;		/* length of hi-low */
    int density_plot;		/* if type is XYZ then density_plot  = 1 */
    double zmin;
    double zmax;		/* min max for density plots */
    char comments[256];		/* how did this set originate */
} plotarr;

typedef struct {
    int active;			/* legend on or off */
    int loctype;		/* locate in world or viewport coords */
    int layout;			/* verticle or horizontal */
    int vgap;			/* verticle gap between entries */
    int hgap;			/* horizontal gap between entries */
    int len;			/* length of line to draw */
    int box;			/* box around legend on or off */
    double legx;		/* location on graph */
    double legy;
    int font;
    double charsize;
    int color;
    int linew;
    int lines;
    int boxfill;		/* legend frame fill toggle */
    int boxfillusing;		/* legend frame fill type */
    int boxfillcolor;		/* legend frame fill color */
    int boxfillpat;		/* legend frame fill pattern */
    int boxlcolor;		/* legend frame line color */
    int boxlinew;		/* legend frame line width */
    int boxlines;		/* legend frame line style */
    plotstr str[MAXPLOT];	/* legend for each of MAXPLOT sets */
} legend;

typedef struct {
    int active;			/* alive or dead */
    int hidden;			/* display or not */
    int label;			/* label graph */
    int type;			/* type of graph */
    int noauto_world;		/* only time this is used is at startup */
    int noauto_tics;		/* only time this is used is at startup */
    int auto_type;		/* */
    int parmsread;		/* was a paramter file read for this graph */
    int revx;
    int revy;			/* reverse mapping for x and y if true */
    int maxplot;		/* max number of sets for this graph */
    int displayin[MAXGRAPH];	/* display this graph's world in graph
				 * displayin */
    Props dp;			/* props to use for displayin */
    world w;			/* world */
    view v;			/* world/view */
    labels labs;		/* title, subtitle, axes labels */
    tickmarks t[6];		/* flags etc. for tickmarks for all axes */
    framep f;			/* type of box around plot */
    legend l;			/* legend for XY data */
    int pointset;		/* if (dsx, dsy) have been set */
    int pt_type;		/* type of locator display */
    double dsx;
    double dsy;			/* locator fixed point */
    int fx;
    int fy;			/* locator format type */
    int px;
    int py;			/* locator precision */
    world_stack ws[MAX_ZOOM_STACK];	/* zoom stack */
    int ws_top;			/* stack pointer */
    int curw;			/* for cycling through the stack */
    int display_flags[50];	/* which items to display - for the side bar */
    int display_fish;		/* display fishes */
    int use_timeclock;		/* timeclock override */
    velocity_scale vl;
    velocity_scale fl;		/* flux */
    map_scale mapscale;
    north_indicator north;
    DisplayGrid grid[MAXGRIDS];
    DisplayBoundary bounds[MAXBOUNDS];
    DisplaySlice sbox[MAXSLICES];
    Zoom_box zbox[MAXZOOMBOXES];
    plotarr p[5];

    DisplayTideStation *tidestat;

    int curtidestat;
    int curtimehist;

    int curgrid;
    int curset;
    int curslice;
    int curflux;
    int curzoom;
    int curbound;

    Isolparms curip;
    int curplaceitem;

} graph;
