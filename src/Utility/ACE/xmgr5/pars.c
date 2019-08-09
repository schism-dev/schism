/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 1 "pars.y"


/* $Id: pars.y,v 1.8 2004/06/16 18:34:10 pturner Exp $
 * 
 * evaluate expressions, commands, parameter files
 * 
 */

#define PARS			/* to overide some defines in defines.h */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#ifndef WIN32
#include <sys/param.h>
#endif
#include <stdarg.h>
 
#include "globals.h"
#include "noxprotos.h"

void set_prop(int gno, ...);
void set_monitor(int monitor, ...);

#ifndef M_PI
#     define M_PI  3.14159265358979323846
#endif

#ifndef TRUE
#     define TRUE 1
#endif

#ifndef FALSE
#     define FALSE 0
#endif

double result, resx, resy;	/* return value if expression */

double drand48(void);
long lrand48(void);
double erf(double arg); /* doesn't seem to be in ANSI C */
double erfc(double arg); /* doesn't seem to be in ANSI C */
double rnorm(double mean, double sdev);
double fx(double x);
double *getvptr(int gno, int setno, int v);
double vmin(double *x, int n);
double vmax(double *x, int n);
void yyerror(char *s);

static int interr;

static double *freelist[100]; 	/* temporary vectors */
static int fcnt;		/* number allocated */

int naxis = 0;	/* current axis */
int curline, curbox, curstring, curleg, curgrid;

int gotbatch, gotparams, gotread; /* these guys attempt to avoid reentrancy problems */
int readtype, readsrc;
char batchfile[256], paramfile[256], readfile[256];

static char f_string[512];	/* buffer for string to parse */
static int pos = 0;
static double *aa, *bb, *cc, *dd, *xx, *yy;
static int setindex, lxy, ls;
static int setsetno;
static int whichgraph;
static int whichset;

extern int change_gno;
extern int change_type;


/* Line 371 of yacc.c  */
#line 144 "y.tab.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ABS = 258,
     ACOS = 259,
     ALL = 260,
     ASIN = 261,
     ATAN = 262,
     ATAN2 = 263,
     CEIL = 264,
     COS = 265,
     DEG = 266,
     DX = 267,
     DY = 268,
     ERF = 269,
     ERFC = 270,
     EXP = 271,
     FLOOR = 272,
     HYPOT = 273,
     INDEX = 274,
     INT = 275,
     INVN = 276,
     INVT = 277,
     IRAND = 278,
     LGAMMA = 279,
     LN = 280,
     LOG = 281,
     LOGISTIC = 282,
     MAXP = 283,
     MINP = 284,
     MOD = 285,
     NORM = 286,
     NORMP = 287,
     PI = 288,
     RAD = 289,
     RAND = 290,
     RNORM = 291,
     SETNO = 292,
     SIN = 293,
     SQR = 294,
     SQRT = 295,
     TAN = 296,
     INUM = 297,
     CTD = 298,
     ADP = 299,
     TPC = 300,
     VX1 = 301,
     VX2 = 302,
     VY1 = 303,
     VY2 = 304,
     WX1 = 305,
     WX2 = 306,
     WY1 = 307,
     WY2 = 308,
     DELAYP = 309,
     DOUBLEBUFFER = 310,
     DOWN = 311,
     ABSOLUTE = 312,
     ABORT = 313,
     ACTIVATE = 314,
     ACTIVE = 315,
     ALT = 316,
     ALTERNATE = 317,
     ALTXAXIS = 318,
     ALTYAXIS = 319,
     ANGLE = 320,
     ANNOTATE = 321,
     APPEND = 322,
     AREA = 323,
     ARROW = 324,
     AUTO = 325,
     AUTOSCALE = 326,
     AUTOTICKS = 327,
     AVG = 328,
     AXIS = 329,
     AXES = 330,
     BACKBUFFER = 331,
     BACKGROUND = 332,
     BAR = 333,
     BATCH = 334,
     BLOCK = 335,
     BIN = 336,
     BOTH = 337,
     BOTTOM = 338,
     BOX = 339,
     BOXPLOT = 340,
     BP = 341,
     CD = 342,
     CELLS = 343,
     CENTER = 344,
     CHAR = 345,
     CHRSTR = 346,
     CLEAR = 347,
     CLICK = 348,
     CMAP = 349,
     COLOR = 350,
     COMMENT = 351,
     COPY = 352,
     CORIE = 353,
     CYCLE = 354,
     DB = 355,
     DECIMAL = 356,
     DEF = 357,
     DEFAULT = 358,
     DELETE = 359,
     DEVICE = 360,
     DFT = 361,
     DIFFERENCE = 362,
     DISK = 363,
     DRAW2 = 364,
     DROP = 365,
     DXDX = 366,
     DXP = 367,
     DYDY = 368,
     DYP = 369,
     ECHO = 370,
     EDIT = 371,
     ELSE = 372,
     END = 373,
     ERRORBAR = 374,
     EXIT = 375,
     EXPONENTIAL = 376,
     FALSEP = 377,
     FFT = 378,
     FILEP = 379,
     FILL = 380,
     FIND = 381,
     FIXEDPOINT = 382,
     FLUSH = 383,
     FOCUS = 384,
     FOLLOWS = 385,
     FONTP = 386,
     FOREGROUND = 387,
     FORMAT = 388,
     FRONTBUFFER = 389,
     FRAMEP = 390,
     GETP = 391,
     GIFL = 392,
     GIFP = 393,
     GRAPH = 394,
     GRAPHNO = 395,
     GRAPHS = 396,
     GRAPHTYPE = 397,
     GRID = 398,
     HARDCOPY = 399,
     HBAR = 400,
     HBOXPLOT = 401,
     HGAP = 402,
     HIDDEN = 403,
     HORIZONTAL = 404,
     HISTO = 405,
     IF = 406,
     IGNORE = 407,
     IHL = 408,
     IMAGE = 409,
     IN = 410,
     INIT = 411,
     INITGRAPHICS = 412,
     INOUT = 413,
     INTEGRATE = 414,
     INTERP = 415,
     INVDFT = 416,
     INVFFT = 417,
     ISOLINE = 418,
     ISOLINES = 419,
     JUST = 420,
     KILL = 421,
     LABEL = 422,
     LAYOUT = 423,
     LEAVE = 424,
     LEAVEGRAPHICS = 425,
     LEFT = 426,
     LEGEND = 427,
     LENGTH = 428,
     LEVEL = 429,
     LEVELS = 430,
     LINE = 431,
     LINESTYLE = 432,
     LINETO = 433,
     LINEWIDTH = 434,
     LINK = 435,
     LOAD = 436,
     LOCATOR = 437,
     LOCTYPE = 438,
     LOGX = 439,
     LOGY = 440,
     LOGXY = 441,
     MAJOR = 442,
     MIFL = 443,
     MIFP = 444,
     MINOR = 445,
     MISSINGP = 446,
     MONITOR = 447,
     MOVE = 448,
     MOVE2 = 449,
     MOVETO = 450,
     NEGATE = 451,
     NO = 452,
     NONE = 453,
     NORMAL = 454,
     NXY = 455,
     OFF = 456,
     OFFSETX = 457,
     OFFSETY = 458,
     ON = 459,
     OP = 460,
     ORIENT = 461,
     OUT = 462,
     PAGE = 463,
     PARA = 464,
     PARALLEL = 465,
     PARAMETERS = 466,
     PARAMS = 467,
     PATTERN = 468,
     PERIMETER = 469,
     PERP = 470,
     PERPENDICULAR = 471,
     PIE = 472,
     PIPE = 473,
     PLACE = 474,
     POINT = 475,
     POLAR = 476,
     POWER = 477,
     PREC = 478,
     PREPEND = 479,
     PRINT = 480,
     PS = 481,
     PSCOLORP = 482,
     PSMONOP = 483,
     PSCOLORL = 484,
     PSMONOL = 485,
     PUSH = 486,
     POP = 487,
     PUTP = 488,
     READ = 489,
     REDRAW = 490,
     REGRESS = 491,
     REGNUM = 492,
     REGIONS = 493,
     RENDER = 494,
     REVERSE = 495,
     RIGHT = 496,
     RISER = 497,
     ROT = 498,
     RUNAVG = 499,
     RUNMED = 500,
     RUNSTD = 501,
     RUNMIN = 502,
     RUNMAX = 503,
     SAMPLE = 504,
     SAVEALL = 505,
     SCALE = 506,
     SCIENTIFIC = 507,
     SET = 508,
     SETNUM = 509,
     SETS = 510,
     SIGN = 511,
     SIZE = 512,
     SKIP = 513,
     SLEEP = 514,
     SLICE = 515,
     SOURCE = 516,
     SPEC = 517,
     SPECIFIED = 518,
     SPECTRUM = 519,
     STACK = 520,
     STACKEDBAR = 521,
     STACKEDHBAR = 522,
     STACKEDLINE = 523,
     STAGGER = 524,
     START = 525,
     STARTTYPE = 526,
     STATUS = 527,
     STEP = 528,
     STOP = 529,
     STRING = 530,
     SUBTITLE = 531,
     SWAPBUFFER = 532,
     SYMBOL = 533,
     TICKP = 534,
     TICKLABEL = 535,
     TICKMARKS = 536,
     TITLE = 537,
     TO = 538,
     TOP = 539,
     TRUEP = 540,
     TYPE = 541,
     UP = 542,
     VELOCITY = 543,
     VERTICAL = 544,
     VGAP = 545,
     VIEW = 546,
     WITH = 547,
     WORLD = 548,
     WRITE = 549,
     X = 550,
     X0 = 551,
     X1 = 552,
     XAXES = 553,
     XAXIS = 554,
     XCOR = 555,
     XMAX = 556,
     XMIN = 557,
     FEGRID = 558,
     RECTGRID = 559,
     XY = 560,
     XYARC = 561,
     XYBOX = 562,
     XYBOXPLOT = 563,
     XYFIXED = 564,
     XYHILO = 565,
     XYRT = 566,
     XYSEG = 567,
     XYSTRING = 568,
     XYDX = 569,
     XYDY = 570,
     XYDXDX = 571,
     XYDYDY = 572,
     XYDXDY = 573,
     XYX2Y2 = 574,
     XYXX = 575,
     XYYY = 576,
     XYZ = 577,
     XYZW = 578,
     XYUV = 579,
     Y = 580,
     Y0 = 581,
     Y1 = 582,
     Y2 = 583,
     Y3 = 584,
     Y4 = 585,
     Y5 = 586,
     YAXES = 587,
     YAXIS = 588,
     YES = 589,
     YMAX = 590,
     YMIN = 591,
     ZEROXAXIS = 592,
     ZEROYAXIS = 593,
     ABOVE = 594,
     BELOW = 595,
     POLYI = 596,
     POLYO = 597,
     GENERAL = 598,
     DDMMYY = 599,
     MMDDYY = 600,
     MMYY = 601,
     MMDD = 602,
     MONTHDAY = 603,
     DAYMONTH = 604,
     MONTHS = 605,
     MONTHL = 606,
     DDMONTHSYYHHMMSS = 607,
     DDMONTHSYY = 608,
     DAYOFWEEKS = 609,
     DAYOFWEEKL = 610,
     DAYOFYEAR = 611,
     HMS = 612,
     HH = 613,
     MMDDHMS = 614,
     MMDDYYHMS = 615,
     DEGREESLON = 616,
     DEGREESMMLON = 617,
     DEGREESMMSSLON = 618,
     MMSSLON = 619,
     DEGREESLAT = 620,
     DEGREESMMLAT = 621,
     DEGREESMMSSLAT = 622,
     MMSSLAT = 623,
     DOT = 624,
     STAR = 625,
     PLUS = 626,
     CROSS = 627,
     CIRCLE = 628,
     SQUARE = 629,
     DIAMOND = 630,
     TRIANGLE1 = 631,
     TRIANGLE2 = 632,
     TRIANGLE3 = 633,
     TRIANGLE4 = 634,
     SPLINE = 635,
     LANDSCAPE = 636,
     PORTRAIT = 637,
     FREE = 638,
     FIXED = 639,
     STATUSBAR = 640,
     LOCATORBAR = 641,
     TOOLBAR = 642,
     ELCIRC = 643,
     SCALAR = 644,
     VECTOR = 645,
     HEAT = 646,
     HISTORY = 647,
     PROFILE = 648,
     NODE = 649,
     VAR = 650,
     NUMBER = 651,
     OR = 652,
     AND = 653,
     NE = 654,
     EQ = 655,
     GE = 656,
     LE = 657,
     LT = 658,
     GT = 659,
     NOT = 660,
     UMINUS = 661
   };
#endif
/* Tokens.  */
#define ABS 258
#define ACOS 259
#define ALL 260
#define ASIN 261
#define ATAN 262
#define ATAN2 263
#define CEIL 264
#define COS 265
#define DEG 266
#define DX 267
#define DY 268
#define ERF 269
#define ERFC 270
#define EXP 271
#define FLOOR 272
#define HYPOT 273
#define INDEX 274
#define INT 275
#define INVN 276
#define INVT 277
#define IRAND 278
#define LGAMMA 279
#define LN 280
#define LOG 281
#define LOGISTIC 282
#define MAXP 283
#define MINP 284
#define MOD 285
#define NORM 286
#define NORMP 287
#define PI 288
#define RAD 289
#define RAND 290
#define RNORM 291
#define SETNO 292
#define SIN 293
#define SQR 294
#define SQRT 295
#define TAN 296
#define INUM 297
#define CTD 298
#define ADP 299
#define TPC 300
#define VX1 301
#define VX2 302
#define VY1 303
#define VY2 304
#define WX1 305
#define WX2 306
#define WY1 307
#define WY2 308
#define DELAYP 309
#define DOUBLEBUFFER 310
#define DOWN 311
#define ABSOLUTE 312
#define ABORT 313
#define ACTIVATE 314
#define ACTIVE 315
#define ALT 316
#define ALTERNATE 317
#define ALTXAXIS 318
#define ALTYAXIS 319
#define ANGLE 320
#define ANNOTATE 321
#define APPEND 322
#define AREA 323
#define ARROW 324
#define AUTO 325
#define AUTOSCALE 326
#define AUTOTICKS 327
#define AVG 328
#define AXIS 329
#define AXES 330
#define BACKBUFFER 331
#define BACKGROUND 332
#define BAR 333
#define BATCH 334
#define BLOCK 335
#define BIN 336
#define BOTH 337
#define BOTTOM 338
#define BOX 339
#define BOXPLOT 340
#define BP 341
#define CD 342
#define CELLS 343
#define CENTER 344
#define CHAR 345
#define CHRSTR 346
#define CLEAR 347
#define CLICK 348
#define CMAP 349
#define COLOR 350
#define COMMENT 351
#define COPY 352
#define CORIE 353
#define CYCLE 354
#define DB 355
#define DECIMAL 356
#define DEF 357
#define DEFAULT 358
#define DELETE 359
#define DEVICE 360
#define DFT 361
#define DIFFERENCE 362
#define DISK 363
#define DRAW2 364
#define DROP 365
#define DXDX 366
#define DXP 367
#define DYDY 368
#define DYP 369
#define ECHO 370
#define EDIT 371
#define ELSE 372
#define END 373
#define ERRORBAR 374
#define EXIT 375
#define EXPONENTIAL 376
#define FALSEP 377
#define FFT 378
#define FILEP 379
#define FILL 380
#define FIND 381
#define FIXEDPOINT 382
#define FLUSH 383
#define FOCUS 384
#define FOLLOWS 385
#define FONTP 386
#define FOREGROUND 387
#define FORMAT 388
#define FRONTBUFFER 389
#define FRAMEP 390
#define GETP 391
#define GIFL 392
#define GIFP 393
#define GRAPH 394
#define GRAPHNO 395
#define GRAPHS 396
#define GRAPHTYPE 397
#define GRID 398
#define HARDCOPY 399
#define HBAR 400
#define HBOXPLOT 401
#define HGAP 402
#define HIDDEN 403
#define HORIZONTAL 404
#define HISTO 405
#define IF 406
#define IGNORE 407
#define IHL 408
#define IMAGE 409
#define IN 410
#define INIT 411
#define INITGRAPHICS 412
#define INOUT 413
#define INTEGRATE 414
#define INTERP 415
#define INVDFT 416
#define INVFFT 417
#define ISOLINE 418
#define ISOLINES 419
#define JUST 420
#define KILL 421
#define LABEL 422
#define LAYOUT 423
#define LEAVE 424
#define LEAVEGRAPHICS 425
#define LEFT 426
#define LEGEND 427
#define LENGTH 428
#define LEVEL 429
#define LEVELS 430
#define LINE 431
#define LINESTYLE 432
#define LINETO 433
#define LINEWIDTH 434
#define LINK 435
#define LOAD 436
#define LOCATOR 437
#define LOCTYPE 438
#define LOGX 439
#define LOGY 440
#define LOGXY 441
#define MAJOR 442
#define MIFL 443
#define MIFP 444
#define MINOR 445
#define MISSINGP 446
#define MONITOR 447
#define MOVE 448
#define MOVE2 449
#define MOVETO 450
#define NEGATE 451
#define NO 452
#define NONE 453
#define NORMAL 454
#define NXY 455
#define OFF 456
#define OFFSETX 457
#define OFFSETY 458
#define ON 459
#define OP 460
#define ORIENT 461
#define OUT 462
#define PAGE 463
#define PARA 464
#define PARALLEL 465
#define PARAMETERS 466
#define PARAMS 467
#define PATTERN 468
#define PERIMETER 469
#define PERP 470
#define PERPENDICULAR 471
#define PIE 472
#define PIPE 473
#define PLACE 474
#define POINT 475
#define POLAR 476
#define POWER 477
#define PREC 478
#define PREPEND 479
#define PRINT 480
#define PS 481
#define PSCOLORP 482
#define PSMONOP 483
#define PSCOLORL 484
#define PSMONOL 485
#define PUSH 486
#define POP 487
#define PUTP 488
#define READ 489
#define REDRAW 490
#define REGRESS 491
#define REGNUM 492
#define REGIONS 493
#define RENDER 494
#define REVERSE 495
#define RIGHT 496
#define RISER 497
#define ROT 498
#define RUNAVG 499
#define RUNMED 500
#define RUNSTD 501
#define RUNMIN 502
#define RUNMAX 503
#define SAMPLE 504
#define SAVEALL 505
#define SCALE 506
#define SCIENTIFIC 507
#define SET 508
#define SETNUM 509
#define SETS 510
#define SIGN 511
#define SIZE 512
#define SKIP 513
#define SLEEP 514
#define SLICE 515
#define SOURCE 516
#define SPEC 517
#define SPECIFIED 518
#define SPECTRUM 519
#define STACK 520
#define STACKEDBAR 521
#define STACKEDHBAR 522
#define STACKEDLINE 523
#define STAGGER 524
#define START 525
#define STARTTYPE 526
#define STATUS 527
#define STEP 528
#define STOP 529
#define STRING 530
#define SUBTITLE 531
#define SWAPBUFFER 532
#define SYMBOL 533
#define TICKP 534
#define TICKLABEL 535
#define TICKMARKS 536
#define TITLE 537
#define TO 538
#define TOP 539
#define TRUEP 540
#define TYPE 541
#define UP 542
#define VELOCITY 543
#define VERTICAL 544
#define VGAP 545
#define VIEW 546
#define WITH 547
#define WORLD 548
#define WRITE 549
#define X 550
#define X0 551
#define X1 552
#define XAXES 553
#define XAXIS 554
#define XCOR 555
#define XMAX 556
#define XMIN 557
#define FEGRID 558
#define RECTGRID 559
#define XY 560
#define XYARC 561
#define XYBOX 562
#define XYBOXPLOT 563
#define XYFIXED 564
#define XYHILO 565
#define XYRT 566
#define XYSEG 567
#define XYSTRING 568
#define XYDX 569
#define XYDY 570
#define XYDXDX 571
#define XYDYDY 572
#define XYDXDY 573
#define XYX2Y2 574
#define XYXX 575
#define XYYY 576
#define XYZ 577
#define XYZW 578
#define XYUV 579
#define Y 580
#define Y0 581
#define Y1 582
#define Y2 583
#define Y3 584
#define Y4 585
#define Y5 586
#define YAXES 587
#define YAXIS 588
#define YES 589
#define YMAX 590
#define YMIN 591
#define ZEROXAXIS 592
#define ZEROYAXIS 593
#define ABOVE 594
#define BELOW 595
#define POLYI 596
#define POLYO 597
#define GENERAL 598
#define DDMMYY 599
#define MMDDYY 600
#define MMYY 601
#define MMDD 602
#define MONTHDAY 603
#define DAYMONTH 604
#define MONTHS 605
#define MONTHL 606
#define DDMONTHSYYHHMMSS 607
#define DDMONTHSYY 608
#define DAYOFWEEKS 609
#define DAYOFWEEKL 610
#define DAYOFYEAR 611
#define HMS 612
#define HH 613
#define MMDDHMS 614
#define MMDDYYHMS 615
#define DEGREESLON 616
#define DEGREESMMLON 617
#define DEGREESMMSSLON 618
#define MMSSLON 619
#define DEGREESLAT 620
#define DEGREESMMLAT 621
#define DEGREESMMSSLAT 622
#define MMSSLAT 623
#define DOT 624
#define STAR 625
#define PLUS 626
#define CROSS 627
#define CIRCLE 628
#define SQUARE 629
#define DIAMOND 630
#define TRIANGLE1 631
#define TRIANGLE2 632
#define TRIANGLE3 633
#define TRIANGLE4 634
#define SPLINE 635
#define LANDSCAPE 636
#define PORTRAIT 637
#define FREE 638
#define FIXED 639
#define STATUSBAR 640
#define LOCATORBAR 641
#define TOOLBAR 642
#define ELCIRC 643
#define SCALAR 644
#define VECTOR 645
#define HEAT 646
#define HISTORY 647
#define PROFILE 648
#define NODE 649
#define VAR 650
#define NUMBER 651
#define OR 652
#define AND 653
#define NE 654
#define EQ 655
#define GE 656
#define LE 657
#define LT 658
#define GT 659
#define NOT 660
#define UMINUS 661



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 387 of yacc.c  */
#line 77 "pars.y"

    double val;
    long ival;
    double *ptr;
    long func;
    long pset;
    char *str;


/* Line 387 of yacc.c  */
#line 1009 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 1037 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  513
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   7338

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  423
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  49
/* YYNRULES -- Number of rules.  */
#define YYNRULES  733
/* YYNRULES -- Number of states.  */
#define YYNSTATES  1647

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   661

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     414,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,   410,     2,     2,
     417,   418,   408,   406,   415,   407,   416,   409,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   422,     2,
       2,   397,     2,   421,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   419,     2,   420,   411,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   398,   399,   400,   401,   402,   403,   404,   405,
     412,   413
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    16,    19,    22,
      25,    28,    31,    34,    37,    40,    43,    47,    52,    56,
      59,    63,    65,    69,    71,    73,    75,    77,    79,    81,
      83,    85,    88,    92,    96,   100,   104,   114,   120,   125,
     127,   131,   134,   137,   141,   150,   153,   155,   159,   162,
     165,   168,   173,   178,   185,   188,   191,   194,   202,   205,
     208,   212,   216,   219,   222,   225,   227,   230,   233,   236,
     239,   242,   245,   248,   253,   260,   264,   270,   274,   279,
     286,   292,   300,   306,   314,   319,   328,   333,   342,   345,
     348,   351,   354,   356,   364,   371,   378,   385,   396,   407,
     414,   419,   421,   424,   427,   430,   433,   436,   439,   442,
     445,   448,   451,   454,   457,   461,   465,   470,   474,   478,
     483,   487,   492,   512,   521,   534,   543,   547,   553,   556,
     561,   566,   573,   576,   578,   580,   582,   585,   603,   606,
     609,   612,   616,   619,   622,   631,   635,   639,   643,   647,
     651,   656,   661,   664,   667,   671,   674,   677,   680,   689,
     693,   697,   701,   705,   709,   714,   719,   722,   725,   728,
     732,   735,   738,   743,   747,   751,   755,   759,   763,   767,
     772,   776,   780,   784,   788,   793,   797,   802,   807,   816,
     820,   824,   828,   832,   841,   845,   849,   853,   857,   860,
     864,   868,   872,   876,   879,   883,   887,   891,   895,   900,
     903,   907,   911,   915,   919,   923,   927,   932,   938,   944,
     949,   954,   959,   964,   968,   972,   977,   981,   985,   989,
     993,   998,  1001,  1005,  1010,  1015,  1020,  1024,  1030,  1034,
    1039,  1045,  1050,  1055,  1061,  1067,  1073,  1080,  1087,  1090,
    1094,  1098,  1102,  1106,  1110,  1115,  1118,  1122,  1127,  1132,
    1136,  1140,  1144,  1149,  1154,  1158,  1164,  1171,  1178,  1183,
    1188,  1193,  1209,  1222,  1236,  1255,  1261,  1265,  1269,  1276,
    1281,  1286,  1290,  1294,  1298,  1302,  1307,  1309,  1311,  1313,
    1315,  1317,  1319,  1321,  1323,  1325,  1327,  1329,  1331,  1333,
    1335,  1337,  1339,  1341,  1343,  1345,  1347,  1349,  1351,  1353,
    1355,  1357,  1359,  1361,  1363,  1365,  1367,  1369,  1371,  1373,
    1375,  1377,  1379,  1381,  1383,  1385,  1387,  1389,  1391,  1393,
    1395,  1397,  1399,  1401,  1403,  1406,  1408,  1410,  1413,  1416,
    1419,  1422,  1425,  1429,  1433,  1437,  1441,  1445,  1449,  1453,
    1457,  1460,  1463,  1467,  1471,  1474,  1478,  1482,  1486,  1490,
    1494,  1499,  1504,  1509,  1512,  1514,  1517,  1519,  1522,  1526,
    1529,  1531,  1533,  1535,  1537,  1539,  1541,  1544,  1547,  1550,
    1552,  1555,  1558,  1561,  1564,  1568,  1571,  1574,  1577,  1580,
    1583,  1585,  1587,  1590,  1592,  1595,  1598,  1601,  1604,  1607,
    1610,  1613,  1616,  1619,  1622,  1624,  1627,  1630,  1634,  1638,
    1641,  1644,  1648,  1652,  1656,  1660,  1664,  1668,  1672,  1676,
    1679,  1682,  1685,  1688,  1692,  1694,  1697,  1699,  1702,  1705,
    1708,  1711,  1714,  1717,  1720,  1723,  1726,  1729,  1732,  1735,
    1738,  1741,  1744,  1747,  1750,  1753,  1757,  1761,  1765,  1769,
    1772,  1775,  1779,  1782,  1785,  1788,  1792,  1794,  1797,  1800,
    1803,  1806,  1811,  1814,  1818,  1821,  1824,  1827,  1829,  1832,
    1835,  1838,  1842,  1844,  1846,  1849,  1852,  1855,  1857,  1859,
    1861,  1863,  1865,  1867,  1869,  1871,  1873,  1875,  1877,  1879,
    1881,  1883,  1885,  1887,  1889,  1891,  1893,  1895,  1897,  1899,
    1901,  1903,  1905,  1907,  1909,  1911,  1913,  1915,  1917,  1919,
    1921,  1923,  1925,  1927,  1929,  1931,  1933,  1935,  1937,  1939,
    1941,  1943,  1945,  1947,  1949,  1951,  1953,  1955,  1957,  1959,
    1961,  1963,  1965,  1967,  1969,  1971,  1973,  1975,  1977,  1979,
    1981,  1983,  1985,  1987,  1989,  1991,  1993,  1995,  1997,  1999,
    2001,  2003,  2005,  2007,  2009,  2011,  2013,  2015,  2017,  2019,
    2021,  2023,  2025,  2027,  2029,  2036,  2043,  2052,  2058,  2069,
    2075,  2079,  2083,  2089,  2097,  2101,  2105,  2111,  2119,  2121,
    2123,  2127,  2133,  2135,  2139,  2143,  2147,  2151,  2155,  2159,
    2163,  2167,  2171,  2175,  2179,  2183,  2187,  2191,  2195,  2199,
    2203,  2207,  2211,  2215,  2220,  2225,  2230,  2235,  2240,  2247,
    2252,  2257,  2259,  2264,  2269,  2274,  2279,  2286,  2293,  2300,
    2307,  2309,  2311,  2316,  2321,  2326,  2331,  2336,  2345,  2352,
    2359,  2366,  2371,  2373,  2375,  2377,  2384,  2391,  2398,  2405,
    2410,  2415,  2420,  2425,  2431,  2435,  2439,  2443,  2447,  2451,
    2455,  2459,  2463,  2466,  2470,  2473,  2475,  2480,  2484,  2489,
    2496,  2505,  2511,  2519,  2523,  2529,  2533,  2537,  2541,  2545,
    2549,  2553,  2558,  2563,  2568,  2573,  2580,  2585,  2590,  2592,
    2597,  2602,  2607,  2612,  2619,  2623,  2627,  2631,  2635,  2639,
    2643,  2647,  2651,  2653,  2655,  2657,  2659,  2661,  2663,  2665,
    2667,  2669,  2671,  2673,  2675,  2680,  2685,  2690,  2695,  2700,
    2709,  2716,  2723,  2730,  2735,  2737,  2739,  2741,  2748,  2753,
    2758,  2763,  2768,  2774,  2778,  2782,  2786,  2790,  2794,  2798,
    2802,  2806,  2809,  2813
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     424,     0,    -1,    -1,   467,   414,    -1,   469,   414,    -1,
     468,   414,    -1,   471,   414,    -1,   470,   414,    -1,   428,
     414,    -1,   427,   414,    -1,   438,   414,    -1,   430,   414,
      -1,   429,   414,    -1,   436,   414,    -1,   425,   414,    -1,
       1,   414,    -1,   225,   426,    91,    -1,   225,   283,   426,
      91,    -1,   225,   283,   426,    -1,   105,   396,    -1,   144,
     105,   396,    -1,   144,    -1,   226,   179,   396,    -1,   228,
      -1,   230,    -1,   189,    -1,   188,    -1,   138,    -1,   137,
      -1,   144,    -1,   124,    -1,   237,   451,    -1,   237,   286,
     434,    -1,   237,    95,   396,    -1,   237,   177,   396,    -1,
     237,   179,   396,    -1,   237,   176,   471,   415,   471,   415,
     471,   415,   471,    -1,   237,   305,   471,   415,   471,    -1,
     180,   237,   283,   140,    -1,   235,    -1,    70,   235,   451,
      -1,    87,    91,    -1,   115,    91,    -1,    77,    95,   396,
      -1,    94,   396,   415,   396,   415,   396,   415,   396,    -1,
      98,   460,    -1,   120,    -1,   208,   168,   433,    -1,   385,
     451,    -1,   387,   451,    -1,   386,   451,    -1,   109,   471,
     415,   471,    -1,   194,   471,   415,   471,    -1,   278,   471,
     415,   471,   415,   396,    -1,    95,   396,    -1,   179,   396,
      -1,   177,   396,    -1,   294,   275,   471,   415,   471,   415,
      91,    -1,   208,   464,    -1,   208,   396,    -1,   208,   158,
     396,    -1,   180,   208,   451,    -1,    55,   460,    -1,   134,
     460,    -1,    76,   460,    -1,   277,    -1,   259,   396,    -1,
      54,   396,    -1,    58,   460,    -1,   136,    91,    -1,   233,
      91,    -1,   292,   140,    -1,   292,   254,    -1,   254,   180,
     455,    91,    -1,   140,   416,   254,   180,   455,    91,    -1,
     254,   180,   451,    -1,   140,   416,   254,   180,   451,    -1,
      59,   254,   396,    -1,    59,   254,   431,   396,    -1,    59,
     140,   416,   254,   431,   396,    -1,   254,   220,   471,   415,
     471,    -1,   140,   416,   254,   220,   471,   415,   471,    -1,
     254,   110,   396,   415,   396,    -1,   140,   416,   254,   110,
     396,   415,   396,    -1,    97,   254,   283,   254,    -1,    97,
     140,   416,   254,   283,   140,   416,   254,    -1,   193,   254,
     283,   254,    -1,   193,   140,   416,   254,   283,   140,   416,
     254,    -1,   166,   254,    -1,   166,   255,    -1,   166,   140,
      -1,   166,   141,    -1,   128,    -1,   181,   395,   396,   415,
     471,   415,   471,    -1,   236,   417,   254,   415,   396,   418,
      -1,   453,   417,   254,   415,   396,   418,    -1,   454,   417,
     254,   415,   396,   418,    -1,   380,   417,   254,   415,   471,
     415,   471,   415,   396,   418,    -1,   150,   417,   254,   415,
     471,   415,   471,   415,   396,   418,    -1,   107,   417,   254,
     415,   396,   418,    -1,    20,   417,   254,   418,    -1,    71,
      -1,    71,   298,    -1,    71,   332,    -1,    71,   254,    -1,
     182,   451,    -1,   129,   140,    -1,   129,   451,    -1,   129,
     253,    -1,   129,   130,    -1,   129,    93,    -1,   261,   455,
      -1,   286,   431,    -1,   234,    91,    -1,   234,    79,    91,
      -1,   234,    80,    91,    -1,   234,    80,   455,    91,    -1,
      80,   431,    91,    -1,   234,   431,    91,    -1,   234,   431,
     455,    91,    -1,   234,   154,    91,    -1,   234,    43,     5,
      91,    -1,   234,   388,   140,   416,   254,    91,   435,   392,
     396,   415,   396,   415,   396,   415,   396,   415,   396,   191,
     396,    -1,   234,   388,   254,    91,   435,   392,   415,   396,
      -1,   234,   388,   140,   416,   254,    91,   435,   393,   415,
     396,   191,   396,    -1,   234,   388,   254,    91,   435,   392,
     415,   396,    -1,   294,   154,    91,    -1,   154,   305,   396,
     415,   396,    -1,   294,   254,    -1,   294,   254,   133,    91,
      -1,   294,   254,   124,    91,    -1,   294,   254,   124,    91,
     133,    91,    -1,   250,    91,    -1,   231,    -1,   232,    -1,
      99,    -1,   265,   396,    -1,   265,   293,   471,   415,   471,
     415,   471,   415,   471,   279,   471,   415,   471,   415,   471,
     415,   471,    -1,    92,   265,    -1,    92,    84,    -1,   292,
      84,    -1,   292,    84,   396,    -1,    84,   451,    -1,    84,
     140,    -1,    84,   471,   415,   471,   415,   471,   415,   471,
      -1,    84,   183,   465,    -1,    84,   177,   396,    -1,    84,
     179,   396,    -1,    84,    95,   396,    -1,    84,   125,   456,
      -1,    84,   125,    95,   396,    -1,    84,   125,   213,   396,
      -1,    84,   102,    -1,   292,   176,    -1,   292,   176,   396,
      -1,    92,   176,    -1,   176,   451,    -1,   176,   140,    -1,
     176,   471,   415,   471,   415,   471,   415,   471,    -1,   176,
     183,   465,    -1,   176,   179,   396,    -1,   176,   177,   396,
      -1,   176,    95,   396,    -1,   176,    69,   396,    -1,   176,
      69,   257,   396,    -1,   176,    69,   286,   396,    -1,   176,
     102,    -1,    92,   275,    -1,   292,   275,    -1,   292,   275,
     396,    -1,   275,   451,    -1,   275,   140,    -1,   275,   471,
     415,   471,    -1,   275,   183,   465,    -1,   275,   179,   396,
      -1,   275,    95,   396,    -1,   275,   243,   396,    -1,   275,
     131,   396,    -1,   275,   165,   396,    -1,   275,    90,   257,
     396,    -1,   275,   102,    91,    -1,   103,   177,   396,    -1,
     103,   179,   396,    -1,   103,    95,   396,    -1,   103,    90,
     257,   396,    -1,   103,   131,   396,    -1,   103,   131,   261,
     396,    -1,   103,   278,   257,   396,    -1,   293,   471,   415,
     471,   415,   471,   415,   471,    -1,   293,   302,   471,    -1,
     293,   301,   471,    -1,   293,   336,   471,    -1,   293,   335,
     471,    -1,   291,   471,   415,   471,   415,   471,   415,   471,
      -1,   291,   302,   396,    -1,   291,   301,   396,    -1,   291,
     336,   396,    -1,   291,   335,   396,    -1,   282,    91,    -1,
     282,   131,   396,    -1,   282,   257,   396,    -1,   282,    95,
     396,    -1,   282,   179,   396,    -1,   276,    91,    -1,   276,
     131,   396,    -1,   276,   257,   396,    -1,   276,    95,   396,
      -1,   276,   179,   396,    -1,   141,    28,   255,   396,    -1,
     172,   451,    -1,   172,   183,   465,    -1,   172,   168,   396,
      -1,   172,   290,   396,    -1,   172,   147,   396,    -1,   172,
     173,   396,    -1,   172,    84,   451,    -1,   172,    84,   125,
     451,    -1,   172,    84,   125,   292,   452,    -1,   172,    84,
     125,   452,   396,    -1,   172,    84,    95,   396,    -1,   172,
      84,   177,   396,    -1,   172,    84,   179,   396,    -1,   172,
     471,   415,   471,    -1,   172,   297,   471,    -1,   172,   327,
     471,    -1,   172,    90,   257,   396,    -1,   172,   131,   396,
      -1,   172,   177,   396,    -1,   172,   179,   396,    -1,   172,
      95,   396,    -1,   172,   275,   396,    91,    -1,   164,   396,
      -1,   164,   286,   396,    -1,   163,   396,   415,   471,    -1,
     163,   396,   179,   396,    -1,   163,   396,   177,   396,    -1,
     163,    95,   396,    -1,   164,   270,   471,   273,   471,    -1,
     164,   172,   451,    -1,   164,   172,   183,   465,    -1,   164,
     172,   396,   415,   471,    -1,   164,   172,   168,   149,    -1,
     164,   172,   168,   289,    -1,   164,   172,   167,   133,   462,
      -1,   164,   172,   167,   223,   396,    -1,   164,   172,   167,
     257,   396,    -1,   164,   172,   257,   396,   415,   396,    -1,
     164,   172,   147,   396,   415,   396,    -1,   135,   451,    -1,
     135,   286,   396,    -1,   135,   177,   396,    -1,   135,   179,
     396,    -1,   135,    95,   396,    -1,   135,   125,   451,    -1,
     135,    77,    95,   396,    -1,   140,   451,    -1,   140,   167,
     451,    -1,   140,    71,   286,    70,    -1,   140,    71,   286,
     262,    -1,   140,    71,   460,    -1,   140,   148,   460,    -1,
     140,   286,   432,    -1,   140,    78,   257,   396,    -1,   140,
     266,   257,   396,    -1,   140,   127,   451,    -1,   140,   127,
     133,   462,   462,    -1,   140,   127,   223,   396,   415,   396,
      -1,   140,   127,   305,   471,   415,   471,    -1,   140,   127,
     286,   396,    -1,   140,    28,   255,   396,    -1,   100,    91,
     415,    91,    -1,   100,   140,   416,   254,    91,   415,    91,
     415,    91,   415,    91,   415,   396,   415,   396,    -1,   100,
      91,   415,    91,   415,    91,   415,    91,   415,   396,   415,
     396,    -1,   100,   140,    91,   415,    91,   415,    91,   415,
      91,   415,   396,   415,   396,    -1,   100,   140,   416,   254,
      44,    91,   415,    91,   415,    91,   415,    91,   415,   396,
     415,   396,   415,   396,    -1,   100,   140,   416,   254,    91,
      -1,   288,   251,   396,    -1,   288,   172,   451,    -1,   288,
     172,   305,   396,   415,   396,    -1,   288,   172,   173,   396,
      -1,   288,   172,   167,    91,    -1,   288,   179,   396,    -1,
     288,   177,   396,    -1,   288,    95,   396,    -1,   288,   183,
     465,    -1,   288,    69,   286,   396,    -1,   305,    -1,   306,
      -1,   307,    -1,   308,    -1,   310,    -1,   311,    -1,   312,
      -1,   313,    -1,   314,    -1,   315,    -1,   316,    -1,   317,
      -1,   318,    -1,   319,    -1,   320,    -1,   321,    -1,   322,
      -1,   323,    -1,   324,    -1,   304,    -1,   303,    -1,   200,
      -1,    81,    -1,   305,    -1,   184,    -1,   185,    -1,   186,
      -1,    78,    -1,   145,    -1,   266,    -1,   267,    -1,   221,
      -1,   309,    -1,    85,    -1,   146,    -1,   383,    -1,   381,
      -1,   382,    -1,   384,    -1,   339,    -1,   340,    -1,   171,
      -1,   241,    -1,   341,    -1,   342,    -1,   389,    -1,   390,
      -1,   391,    -1,   449,   437,    -1,   451,    -1,   152,    -1,
     286,   431,    -1,   191,   471,    -1,   223,   396,    -1,   133,
     462,    -1,   278,   471,    -1,   278,   125,   396,    -1,   278,
      89,   460,    -1,   278,   257,   471,    -1,   278,    90,   396,
      -1,   278,   258,   396,    -1,   278,    95,   471,    -1,   278,
     179,   396,    -1,   278,   177,   396,    -1,   450,   396,    -1,
     125,   396,    -1,   125,   292,   452,    -1,   125,   452,   396,
      -1,   258,   396,    -1,   119,   173,   396,    -1,   119,   286,
     457,    -1,   119,   179,   396,    -1,   119,   177,   396,    -1,
     119,   242,   451,    -1,   119,   242,   179,   396,    -1,   119,
     242,   177,   396,    -1,   322,   471,   415,   471,    -1,    96,
      91,    -1,   439,    -1,   439,   442,    -1,   440,    -1,   141,
     439,    -1,   141,   439,   442,    -1,   141,   440,    -1,   299,
      -1,   333,    -1,    63,    -1,    64,    -1,   337,    -1,   338,
      -1,    75,   441,    -1,   298,   441,    -1,   332,   441,    -1,
     451,    -1,    95,   396,    -1,   179,   396,    -1,   177,   396,
      -1,   131,   396,    -1,    90,   257,   396,    -1,   143,   451,
      -1,   279,   443,    -1,   280,   445,    -1,   167,   447,    -1,
      78,   448,    -1,   451,    -1,   444,    -1,   443,   444,    -1,
     451,    -1,   187,   451,    -1,   190,   451,    -1,   187,   471,
      -1,   190,   471,    -1,   202,   471,    -1,   203,   471,    -1,
      61,   451,    -1,    29,   471,    -1,    28,   471,    -1,   103,
     396,    -1,   461,    -1,    26,   451,    -1,   257,   396,    -1,
     187,   257,   396,    -1,   190,   257,   396,    -1,    95,   396,
      -1,   179,   396,    -1,   187,    95,   396,    -1,   190,    95,
     396,    -1,   187,   179,   396,    -1,   190,   179,   396,    -1,
     187,   177,   396,    -1,   190,   177,   396,    -1,   187,   143,
     451,    -1,   190,   143,   451,    -1,   205,   457,    -1,   286,
      70,    -1,   286,   262,    -1,   262,   396,    -1,   396,   415,
     471,    -1,   446,    -1,   445,   446,    -1,   451,    -1,   286,
      70,    -1,   286,   262,    -1,   223,   396,    -1,   133,   462,
      -1,   133,   396,    -1,    67,    91,    -1,   224,    91,    -1,
     168,   149,    -1,   168,   289,    -1,   168,   262,    -1,    65,
     396,    -1,   165,   458,    -1,   258,   396,    -1,   269,   396,
      -1,   205,   457,    -1,   256,   463,    -1,   270,   471,    -1,
     274,   471,    -1,   270,   286,   262,    -1,   270,   286,    70,
      -1,   274,   286,   262,    -1,   274,   286,    70,    -1,   290,
     396,    -1,   147,   396,    -1,    90,   257,   396,    -1,   131,
     396,    -1,    95,   396,    -1,   179,   396,    -1,   396,   415,
      91,    -1,    91,    -1,   168,   215,    -1,   168,   209,    -1,
     219,    70,    -1,   219,   262,    -1,   219,   396,   415,   396,
      -1,   165,   458,    -1,    90,   257,   396,    -1,   131,   396,
      -1,    95,   396,    -1,   179,   396,    -1,   451,    -1,    95,
     396,    -1,   177,   396,    -1,   179,   396,    -1,   140,   416,
     254,    -1,   254,    -1,   255,    -1,   140,   255,    -1,   141,
     255,    -1,   141,   254,    -1,   177,    -1,   179,    -1,   131,
      -1,    95,    -1,   257,    -1,   204,    -1,   201,    -1,    95,
      -1,   213,    -1,   244,    -1,   246,    -1,   245,    -1,   248,
      -1,   247,    -1,   106,    -1,   123,    -1,   161,    -1,   162,
      -1,   108,    -1,   218,    -1,   213,    -1,    95,    -1,   198,
      -1,   284,    -1,    83,    -1,   171,    -1,   241,    -1,    82,
      -1,   241,    -1,   171,    -1,    89,    -1,    29,    -1,    28,
      -1,   285,    -1,   122,    -1,   155,    -1,   207,    -1,    82,
      -1,   101,    -1,   121,    -1,   222,    -1,   343,    -1,   344,
      -1,   345,    -1,   346,    -1,   347,    -1,   348,    -1,   349,
      -1,   352,    -1,   353,    -1,   350,    -1,   351,    -1,   354,
      -1,   355,    -1,   356,    -1,   357,    -1,   358,    -1,   359,
      -1,   360,    -1,   361,    -1,   362,    -1,   363,    -1,   364,
      -1,   365,    -1,   366,    -1,   367,    -1,   368,    -1,   199,
      -1,    57,    -1,   196,    -1,   287,    -1,    56,    -1,   241,
      -1,   171,    -1,   155,    -1,   207,    -1,   293,    -1,   291,
      -1,   295,    -1,   325,    -1,   296,    -1,   326,    -1,   327,
      -1,   328,    -1,   329,    -1,   330,    -1,   331,    -1,   395,
     419,   471,   420,   397,   471,    -1,   466,   419,   471,   420,
     397,   471,    -1,   254,   416,   466,   419,   471,   420,   397,
     471,    -1,   254,   416,   278,   397,   471,    -1,   140,   416,
     254,   416,   466,   419,   471,   420,   397,   471,    -1,   237,
     416,   466,   397,   471,    -1,   395,   397,   470,    -1,   466,
     397,   470,    -1,   254,   416,   466,   397,   470,    -1,   140,
     416,   254,   416,   466,   397,   470,    -1,   395,   397,   471,
      -1,   466,   397,   471,    -1,   254,   416,   466,   397,   471,
      -1,   140,   416,   254,   416,   466,   397,   471,    -1,   395,
      -1,   466,    -1,   254,   416,   466,    -1,   140,   416,   254,
     416,   466,    -1,   471,    -1,   471,   406,   471,    -1,   470,
     406,   470,    -1,   471,   406,   470,    -1,   470,   406,   471,
      -1,   471,   407,   471,    -1,   470,   407,   470,    -1,   471,
     407,   470,    -1,   470,   407,   471,    -1,   471,   408,   471,
      -1,   470,   408,   470,    -1,   471,   408,   470,    -1,   470,
     408,   471,    -1,   471,   409,   471,    -1,   470,   409,   470,
      -1,   471,   409,   470,    -1,   470,   409,   471,    -1,   471,
     411,   471,    -1,   471,   411,   470,    -1,   470,   411,   471,
      -1,   470,   411,   470,    -1,     3,   417,   471,   418,    -1,
       3,   417,   470,   418,    -1,     4,   417,   470,   418,    -1,
       6,   417,   470,   418,    -1,     7,   417,   470,   418,    -1,
       8,   417,   470,   415,   470,   418,    -1,     9,   417,   470,
     418,    -1,    10,   417,   470,   418,    -1,    11,    -1,    14,
     417,   470,   418,    -1,    15,   417,   470,   418,    -1,    16,
     417,   470,   418,    -1,    17,   417,   470,   418,    -1,    18,
     417,   470,   415,   470,   418,    -1,    18,   417,   471,   415,
     470,   418,    -1,    18,   417,   470,   415,   471,   418,    -1,
      18,   417,   471,   415,   471,   418,    -1,    19,    -1,    37,
      -1,    20,   417,   470,   418,    -1,    23,   417,   396,   418,
      -1,    24,   417,   470,   418,    -1,    25,   417,   470,   418,
      -1,    26,   417,   470,   418,    -1,    27,   417,   470,   415,
     471,   415,   471,   418,    -1,    28,   417,   470,   415,   470,
     418,    -1,    29,   417,   470,   415,   470,   418,    -1,    30,
     417,   470,   415,   470,   418,    -1,    31,   417,   470,   418,
      -1,    33,    -1,    34,    -1,    35,    -1,    36,   417,   470,
     415,   470,   418,    -1,    36,   417,   471,   415,   470,   418,
      -1,    36,   417,   470,   415,   471,   418,    -1,    36,   417,
     471,   415,   471,   418,    -1,    38,   417,   470,   418,    -1,
      39,   417,   470,   418,    -1,    40,   417,   470,   418,    -1,
      41,   417,   470,   418,    -1,   470,   421,   470,   422,   470,
      -1,   470,   405,   470,    -1,   470,   404,   470,    -1,   470,
     403,   470,    -1,   470,   402,   470,    -1,   470,   401,   470,
      -1,   470,   400,   470,    -1,   470,   399,   470,    -1,   470,
     398,   470,    -1,   412,   470,    -1,   417,   470,   418,    -1,
     407,   470,    -1,   396,    -1,   395,   419,   471,   420,    -1,
     237,   416,   466,    -1,   466,   419,   471,   420,    -1,   254,
     416,   466,   419,   471,   420,    -1,   140,   416,   254,   416,
     466,   419,   471,   420,    -1,   254,   416,   466,   416,   459,
      -1,   140,   416,   254,   416,   466,   416,   459,    -1,   254,
     416,   173,    -1,   254,   416,   466,   416,    73,    -1,   471,
     406,   471,    -1,   471,   407,   471,    -1,   471,   408,   471,
      -1,   471,   409,   471,    -1,   471,   410,   471,    -1,   471,
     411,   471,    -1,     3,   417,   471,   418,    -1,     4,   417,
     471,   418,    -1,     6,   417,   471,   418,    -1,     7,   417,
     471,   418,    -1,     8,   417,   471,   415,   471,   418,    -1,
       9,   417,   471,   418,    -1,    10,   417,   471,   418,    -1,
      11,    -1,    14,   417,   471,   418,    -1,    15,   417,   471,
     418,    -1,    16,   417,   471,   418,    -1,    17,   417,   471,
     418,    -1,    18,   417,   471,   415,   471,   418,    -1,   140,
     416,    46,    -1,   140,   416,    47,    -1,   140,   416,    48,
      -1,   140,   416,    49,    -1,   140,   416,    50,    -1,   140,
     416,    51,    -1,   140,   416,    52,    -1,   140,   416,    53,
      -1,    46,    -1,    47,    -1,    48,    -1,    49,    -1,    50,
      -1,    51,    -1,    52,    -1,    53,    -1,    12,    -1,    13,
      -1,    19,    -1,    37,    -1,    20,   417,   471,   418,    -1,
      23,   417,   396,   418,    -1,    24,   417,   471,   418,    -1,
      25,   417,   471,   418,    -1,    26,   417,   471,   418,    -1,
      27,   417,   471,   415,   471,   415,   471,   418,    -1,    28,
     417,   471,   415,   471,   418,    -1,    29,   417,   471,   415,
     471,   418,    -1,    30,   417,   471,   415,   471,   418,    -1,
      31,   417,   471,   418,    -1,    33,    -1,    34,    -1,    35,
      -1,    36,   417,   471,   415,   471,   418,    -1,    38,   417,
     471,   418,    -1,    39,   417,   471,   418,    -1,    40,   417,
     471,   418,    -1,    41,   417,   471,   418,    -1,   471,   421,
     471,   422,   471,    -1,   471,   405,   471,    -1,   471,   404,
     471,    -1,   471,   403,   471,    -1,   471,   402,   471,    -1,
     471,   401,   471,    -1,   471,   400,   471,    -1,   471,   399,
     471,    -1,   471,   398,   471,    -1,   412,   471,    -1,   417,
     471,   418,    -1,   407,   471,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   531,   531,   532,   533,   534,   535,   538,   541,   542,
     543,   544,   545,   546,   547,   548,   554,   563,   572,   580,
     583,   586,   589,   597,   598,   599,   600,   601,   602,   603,
     604,   608,   611,   614,   617,   620,   623,   630,   644,   650,
     653,   656,   668,   677,   680,   683,   686,   689,   693,   697,
     701,   706,   712,   717,   725,   730,   735,   740,   743,   766,
     769,   772,   775,   778,   781,   784,   787,   790,   792,   794,
     800,   813,   818,   823,   827,   831,   834,   838,   841,   848,
     855,   858,   862,   870,   878,   881,   884,   887,   890,   894,
     901,   905,   909,   913,   920,   925,   929,   933,   937,   941,
     945,   949,   963,   975,   987,  1001,  1006,  1016,  1019,  1022,
    1025,  1028,  1031,  1035,  1043,  1049,  1054,  1059,  1064,  1072,
    1080,  1085,  1089,  1092,  1095,  1098,  1101,  1106,  1111,  1115,
    1120,  1125,  1131,  1137,  1140,  1143,  1146,  1150,  1154,  1157,
    1160,  1163,  1166,  1169,  1172,  1181,  1184,  1187,  1190,  1193,
    1196,  1199,  1202,  1214,  1217,  1220,  1223,  1226,  1229,  1236,
    1239,  1242,  1245,  1248,  1251,  1254,  1257,  1269,  1272,  1275,
    1278,  1281,  1284,  1289,  1292,  1295,  1298,  1301,  1304,  1307,
    1310,  1322,  1325,  1328,  1331,  1334,  1337,  1340,  1343,  1350,
    1353,  1356,  1359,  1362,  1369,  1372,  1375,  1378,  1381,  1385,
    1388,  1391,  1394,  1398,  1402,  1406,  1409,  1413,  1417,  1420,
    1423,  1426,  1429,  1432,  1435,  1438,  1441,  1444,  1447,  1455,
    1458,  1461,  1464,  1468,  1471,  1474,  1477,  1480,  1483,  1486,
    1489,  1493,  1496,  1499,  1502,  1504,  1506,  1508,  1512,  1515,
    1518,  1522,  1525,  1528,  1531,  1534,  1537,  1541,  1545,  1548,
    1551,  1554,  1557,  1560,  1563,  1566,  1569,  1572,  1575,  1578,
    1581,  1584,  1587,  1590,  1593,  1596,  1601,  1606,  1611,  1614,
    1621,  1626,  1633,  1640,  1647,  1654,  1663,  1664,  1665,  1669,
    1670,  1671,  1672,  1673,  1674,  1675,  1679,  1682,  1685,  1688,
    1691,  1694,  1697,  1700,  1703,  1706,  1709,  1712,  1715,  1718,
    1721,  1724,  1727,  1730,  1733,  1736,  1739,  1742,  1745,  1751,
    1754,  1757,  1760,  1763,  1766,  1769,  1772,  1775,  1778,  1781,
    1784,  1790,  1793,  1796,  1799,  1805,  1806,  1807,  1808,  1809,
    1810,  1814,  1815,  1816,  1820,  1835,  1838,  1841,  1844,  1847,
    1850,  1853,  1856,  1859,  1862,  1865,  1868,  1871,  1874,  1877,
    1880,  1883,  1886,  1889,  1892,  1895,  1898,  1901,  1904,  1907,
    1910,  1913,  1916,  1919,  1926,  1927,  1928,  1929,  1930,  1931,
    1935,  1936,  1937,  1938,  1939,  1940,  1944,  1945,  1946,  1950,
    1953,  1956,  1959,  1962,  1965,  1968,  1974,  1975,  1976,  1977,
    1978,  1984,  1985,  1989,  1994,  1997,  2000,  2003,  2006,  2009,
    2012,  2015,  2018,  2021,  2024,  2027,  2030,  2033,  2036,  2039,
    2042,  2045,  2048,  2051,  2054,  2057,  2060,  2063,  2066,  2069,
    2072,  2075,  2078,  2081,  2087,  2088,  2092,  2095,  2098,  2101,
    2104,  2107,  2110,  2114,  2118,  2121,  2124,  2127,  2130,  2133,
    2136,  2139,  2142,  2145,  2148,  2151,  2154,  2157,  2160,  2163,
    2166,  2169,  2172,  2175,  2178,  2181,  2188,  2192,  2195,  2198,
    2201,  2204,  2208,  2211,  2214,  2217,  2220,  2226,  2229,  2232,
    2235,  2241,  2246,  2251,  2256,  2261,  2266,  2274,  2277,  2280,
    2283,  2286,  2292,  2295,  2301,  2304,  2310,  2313,  2316,  2319,
    2322,  2328,  2331,  2334,  2337,  2343,  2346,  2352,  2355,  2358,
    2364,  2367,  2370,  2373,  2376,  2382,  2385,  2388,  2394,  2397,
    2403,  2406,  2412,  2415,  2418,  2424,  2427,  2430,  2433,  2436,
    2439,  2442,  2445,  2448,  2451,  2454,  2457,  2460,  2463,  2466,
    2469,  2472,  2475,  2478,  2481,  2484,  2487,  2490,  2493,  2496,
    2499,  2502,  2505,  2508,  2514,  2517,  2520,  2526,  2529,  2532,
    2535,  2538,  2541,  2547,  2550,  2556,  2557,  2558,  2559,  2560,
    2561,  2562,  2563,  2564,  2568,  2579,  2594,  2609,  2621,  2639,
    2650,  2658,  2681,  2703,  2726,  2734,  2757,  2780,  2806,  2815,
    2829,  2843,  2857,  2866,  2875,  2884,  2893,  2902,  2911,  2920,
    2929,  2938,  2947,  2956,  2965,  2974,  2987,  3002,  3017,  3030,
    3039,  3048,  3057,  3066,  3075,  3084,  3093,  3102,  3111,  3120,
    3129,  3138,  3147,  3156,  3165,  3174,  3183,  3192,  3201,  3210,
    3219,  3228,  3237,  3246,  3255,  3264,  3273,  3282,  3291,  3300,
    3309,  3318,  3327,  3336,  3345,  3354,  3363,  3372,  3381,  3390,
    3399,  3408,  3417,  3426,  3438,  3447,  3456,  3465,  3474,  3483,
    3492,  3501,  3510,  3519,  3528,  3538,  3539,  3542,  3545,  3555,
    3565,  3575,  3590,  3605,  3608,  3621,  3624,  3627,  3630,  3639,
    3642,  3645,  3648,  3651,  3654,  3657,  3660,  3663,  3666,  3669,
    3672,  3675,  3678,  3681,  3684,  3687,  3690,  3693,  3696,  3699,
    3702,  3705,  3708,  3711,  3714,  3717,  3720,  3723,  3726,  3729,
    3732,  3735,  3738,  3741,  3744,  3747,  3750,  3753,  3756,  3759,
    3763,  3766,  3769,  3772,  3775,  3778,  3781,  3784,  3787,  3790,
    3793,  3796,  3799,  3806,  3809,  3812,  3815,  3818,  3821,  3824,
    3827,  3830,  3833,  3836
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ABS", "ACOS", "ALL", "ASIN", "ATAN",
  "ATAN2", "CEIL", "COS", "DEG", "DX", "DY", "ERF", "ERFC", "EXP", "FLOOR",
  "HYPOT", "INDEX", "INT", "INVN", "INVT", "IRAND", "LGAMMA", "LN", "LOG",
  "LOGISTIC", "MAXP", "MINP", "MOD", "NORM", "NORMP", "PI", "RAD", "RAND",
  "RNORM", "SETNO", "SIN", "SQR", "SQRT", "TAN", "INUM", "CTD", "ADP",
  "TPC", "VX1", "VX2", "VY1", "VY2", "WX1", "WX2", "WY1", "WY2", "DELAYP",
  "DOUBLEBUFFER", "DOWN", "ABSOLUTE", "ABORT", "ACTIVATE", "ACTIVE", "ALT",
  "ALTERNATE", "ALTXAXIS", "ALTYAXIS", "ANGLE", "ANNOTATE", "APPEND",
  "AREA", "ARROW", "AUTO", "AUTOSCALE", "AUTOTICKS", "AVG", "AXIS", "AXES",
  "BACKBUFFER", "BACKGROUND", "BAR", "BATCH", "BLOCK", "BIN", "BOTH",
  "BOTTOM", "BOX", "BOXPLOT", "BP", "CD", "CELLS", "CENTER", "CHAR",
  "CHRSTR", "CLEAR", "CLICK", "CMAP", "COLOR", "COMMENT", "COPY", "CORIE",
  "CYCLE", "DB", "DECIMAL", "DEF", "DEFAULT", "DELETE", "DEVICE", "DFT",
  "DIFFERENCE", "DISK", "DRAW2", "DROP", "DXDX", "DXP", "DYDY", "DYP",
  "ECHO", "EDIT", "ELSE", "END", "ERRORBAR", "EXIT", "EXPONENTIAL",
  "FALSEP", "FFT", "FILEP", "FILL", "FIND", "FIXEDPOINT", "FLUSH", "FOCUS",
  "FOLLOWS", "FONTP", "FOREGROUND", "FORMAT", "FRONTBUFFER", "FRAMEP",
  "GETP", "GIFL", "GIFP", "GRAPH", "GRAPHNO", "GRAPHS", "GRAPHTYPE",
  "GRID", "HARDCOPY", "HBAR", "HBOXPLOT", "HGAP", "HIDDEN", "HORIZONTAL",
  "HISTO", "IF", "IGNORE", "IHL", "IMAGE", "IN", "INIT", "INITGRAPHICS",
  "INOUT", "INTEGRATE", "INTERP", "INVDFT", "INVFFT", "ISOLINE",
  "ISOLINES", "JUST", "KILL", "LABEL", "LAYOUT", "LEAVE", "LEAVEGRAPHICS",
  "LEFT", "LEGEND", "LENGTH", "LEVEL", "LEVELS", "LINE", "LINESTYLE",
  "LINETO", "LINEWIDTH", "LINK", "LOAD", "LOCATOR", "LOCTYPE", "LOGX",
  "LOGY", "LOGXY", "MAJOR", "MIFL", "MIFP", "MINOR", "MISSINGP", "MONITOR",
  "MOVE", "MOVE2", "MOVETO", "NEGATE", "NO", "NONE", "NORMAL", "NXY",
  "OFF", "OFFSETX", "OFFSETY", "ON", "OP", "ORIENT", "OUT", "PAGE", "PARA",
  "PARALLEL", "PARAMETERS", "PARAMS", "PATTERN", "PERIMETER", "PERP",
  "PERPENDICULAR", "PIE", "PIPE", "PLACE", "POINT", "POLAR", "POWER",
  "PREC", "PREPEND", "PRINT", "PS", "PSCOLORP", "PSMONOP", "PSCOLORL",
  "PSMONOL", "PUSH", "POP", "PUTP", "READ", "REDRAW", "REGRESS", "REGNUM",
  "REGIONS", "RENDER", "REVERSE", "RIGHT", "RISER", "ROT", "RUNAVG",
  "RUNMED", "RUNSTD", "RUNMIN", "RUNMAX", "SAMPLE", "SAVEALL", "SCALE",
  "SCIENTIFIC", "SET", "SETNUM", "SETS", "SIGN", "SIZE", "SKIP", "SLEEP",
  "SLICE", "SOURCE", "SPEC", "SPECIFIED", "SPECTRUM", "STACK",
  "STACKEDBAR", "STACKEDHBAR", "STACKEDLINE", "STAGGER", "START",
  "STARTTYPE", "STATUS", "STEP", "STOP", "STRING", "SUBTITLE",
  "SWAPBUFFER", "SYMBOL", "TICKP", "TICKLABEL", "TICKMARKS", "TITLE", "TO",
  "TOP", "TRUEP", "TYPE", "UP", "VELOCITY", "VERTICAL", "VGAP", "VIEW",
  "WITH", "WORLD", "WRITE", "X", "X0", "X1", "XAXES", "XAXIS", "XCOR",
  "XMAX", "XMIN", "FEGRID", "RECTGRID", "XY", "XYARC", "XYBOX",
  "XYBOXPLOT", "XYFIXED", "XYHILO", "XYRT", "XYSEG", "XYSTRING", "XYDX",
  "XYDY", "XYDXDX", "XYDYDY", "XYDXDY", "XYX2Y2", "XYXX", "XYYY", "XYZ",
  "XYZW", "XYUV", "Y", "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "YAXES",
  "YAXIS", "YES", "YMAX", "YMIN", "ZEROXAXIS", "ZEROYAXIS", "ABOVE",
  "BELOW", "POLYI", "POLYO", "GENERAL", "DDMMYY", "MMDDYY", "MMYY", "MMDD",
  "MONTHDAY", "DAYMONTH", "MONTHS", "MONTHL", "DDMONTHSYYHHMMSS",
  "DDMONTHSYY", "DAYOFWEEKS", "DAYOFWEEKL", "DAYOFYEAR", "HMS", "HH",
  "MMDDHMS", "MMDDYYHMS", "DEGREESLON", "DEGREESMMLON", "DEGREESMMSSLON",
  "MMSSLON", "DEGREESLAT", "DEGREESMMLAT", "DEGREESMMSSLAT", "MMSSLAT",
  "DOT", "STAR", "PLUS", "CROSS", "CIRCLE", "SQUARE", "DIAMOND",
  "TRIANGLE1", "TRIANGLE2", "TRIANGLE3", "TRIANGLE4", "SPLINE",
  "LANDSCAPE", "PORTRAIT", "FREE", "FIXED", "STATUSBAR", "LOCATORBAR",
  "TOOLBAR", "ELCIRC", "SCALAR", "VECTOR", "HEAT", "HISTORY", "PROFILE",
  "NODE", "VAR", "NUMBER", "'='", "OR", "AND", "NE", "EQ", "GE", "LE",
  "LT", "GT", "'+'", "'-'", "'*'", "'/'", "'%'", "'^'", "NOT", "UMINUS",
  "'\\n'", "','", "'.'", "'('", "')'", "'['", "']'", "'?'", "':'",
  "$accept", "list", "setprint", "printer", "regionset", "parmset", "db",
  "setvelocity", "xytype", "graphtype", "pagelayout", "regiontype",
  "elcirctype", "set_setprop", "setprop", "setaxis", "axis", "allaxes",
  "axesprops", "axisfeature", "tickdesc", "tickattr", "ticklabeldesc",
  "ticklabelattr", "axislabeldesc", "axisbardesc", "selectsets", "prop",
  "onoff", "colpat", "runtype", "ffttype", "sourcetype", "filltype",
  "opchoice", "justchoice", "extremetype", "torf", "inoutchoice",
  "formatchoice", "signchoice", "direction", "worldview", "vector", "asgn",
  "rasgn", "vasgn", "vexpr", "expr", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   614,
     615,   616,   617,   618,   619,   620,   621,   622,   623,   624,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     635,   636,   637,   638,   639,   640,   641,   642,   643,   644,
     645,   646,   647,   648,   649,   650,   651,    61,   652,   653,
     654,   655,   656,   657,   658,   659,    43,    45,    42,    47,
      37,    94,   660,   661,    10,    44,    46,    40,    41,    91,
      93,    63,    58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   423,   424,   424,   424,   424,   424,   424,   424,   424,
     424,   424,   424,   424,   424,   424,   425,   425,   425,   425,
     425,   425,   425,   426,   426,   426,   426,   426,   426,   426,
     426,   427,   427,   427,   427,   427,   427,   427,   427,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     428,   428,   428,   428,   428,   428,   428,   428,   428,   428,
     429,   429,   429,   429,   429,   429,   430,   430,   430,   430,
     430,   430,   430,   430,   430,   430,   431,   431,   431,   431,
     431,   431,   431,   431,   431,   431,   431,   431,   431,   431,
     431,   431,   431,   431,   431,   431,   431,   431,   431,   432,
     432,   432,   432,   432,   432,   432,   432,   432,   432,   432,
     432,   433,   433,   433,   433,   434,   434,   434,   434,   434,
     434,   435,   435,   435,   436,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   437,   437,   437,   437,   437,   437,
     437,   437,   437,   437,   438,   438,   438,   438,   438,   438,
     439,   439,   439,   439,   439,   439,   440,   440,   440,   441,
     441,   441,   441,   441,   441,   441,   442,   442,   442,   442,
     442,   443,   443,   444,   444,   444,   444,   444,   444,   444,
     444,   444,   444,   444,   444,   444,   444,   444,   444,   444,
     444,   444,   444,   444,   444,   444,   444,   444,   444,   444,
     444,   444,   444,   444,   445,   445,   446,   446,   446,   446,
     446,   446,   446,   446,   446,   446,   446,   446,   446,   446,
     446,   446,   446,   446,   446,   446,   446,   446,   446,   446,
     446,   446,   446,   446,   446,   446,   447,   447,   447,   447,
     447,   447,   447,   447,   447,   447,   447,   448,   448,   448,
     448,   449,   449,   449,   449,   449,   449,   450,   450,   450,
     450,   450,   451,   451,   452,   452,   453,   453,   453,   453,
     453,   454,   454,   454,   454,   455,   455,   456,   456,   456,
     457,   457,   457,   457,   457,   458,   458,   458,   459,   459,
     460,   460,   461,   461,   461,   462,   462,   462,   462,   462,
     462,   462,   462,   462,   462,   462,   462,   462,   462,   462,
     462,   462,   462,   462,   462,   462,   462,   462,   462,   462,
     462,   462,   462,   462,   463,   463,   463,   464,   464,   464,
     464,   464,   464,   465,   465,   466,   466,   466,   466,   466,
     466,   466,   466,   466,   467,   467,   467,   467,   467,   468,
     469,   469,   469,   469,   469,   469,   469,   469,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     3,     4,     3,     2,
       3,     1,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     3,     3,     3,     3,     9,     5,     4,     1,
       3,     2,     2,     3,     8,     2,     1,     3,     2,     2,
       2,     4,     4,     6,     2,     2,     2,     7,     2,     2,
       3,     3,     2,     2,     2,     1,     2,     2,     2,     2,
       2,     2,     2,     4,     6,     3,     5,     3,     4,     6,
       5,     7,     5,     7,     4,     8,     4,     8,     2,     2,
       2,     2,     1,     7,     6,     6,     6,    10,    10,     6,
       4,     1,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     3,     3,     4,     3,     3,     4,
       3,     4,    19,     8,    12,     8,     3,     5,     2,     4,
       4,     6,     2,     1,     1,     1,     2,    17,     2,     2,
       2,     3,     2,     2,     8,     3,     3,     3,     3,     3,
       4,     4,     2,     2,     3,     2,     2,     2,     8,     3,
       3,     3,     3,     3,     4,     4,     2,     2,     2,     3,
       2,     2,     4,     3,     3,     3,     3,     3,     3,     4,
       3,     3,     3,     3,     4,     3,     4,     4,     8,     3,
       3,     3,     3,     8,     3,     3,     3,     3,     2,     3,
       3,     3,     3,     2,     3,     3,     3,     3,     4,     2,
       3,     3,     3,     3,     3,     3,     4,     5,     5,     4,
       4,     4,     4,     3,     3,     4,     3,     3,     3,     3,
       4,     2,     3,     4,     4,     4,     3,     5,     3,     4,
       5,     4,     4,     5,     5,     5,     6,     6,     2,     3,
       3,     3,     3,     3,     4,     2,     3,     4,     4,     3,
       3,     3,     4,     4,     3,     5,     6,     6,     4,     4,
       4,    15,    12,    13,    18,     5,     3,     3,     6,     4,
       4,     3,     3,     3,     3,     4,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     2,     1,     1,     2,     2,     2,
       2,     2,     3,     3,     3,     3,     3,     3,     3,     3,
       2,     2,     3,     3,     2,     3,     3,     3,     3,     3,
       4,     4,     4,     2,     1,     2,     1,     2,     3,     2,
       1,     1,     1,     1,     1,     1,     2,     2,     2,     1,
       2,     2,     2,     2,     3,     2,     2,     2,     2,     2,
       1,     1,     2,     1,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     1,     2,     2,     3,     3,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       2,     2,     2,     3,     1,     2,     1,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     3,     3,     3,     3,     2,
       2,     3,     2,     2,     2,     3,     1,     2,     2,     2,
       2,     4,     2,     3,     2,     2,     2,     1,     2,     2,
       2,     3,     1,     1,     2,     2,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     6,     6,     8,     5,    10,     5,
       3,     3,     5,     7,     3,     3,     5,     7,     1,     1,
       3,     5,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     4,     4,     4,     4,     4,     6,     4,
       4,     1,     4,     4,     4,     4,     6,     6,     6,     6,
       1,     1,     4,     4,     4,     4,     4,     8,     6,     6,
       6,     4,     1,     1,     1,     6,     6,     6,     6,     4,
       4,     4,     4,     5,     3,     3,     3,     3,     3,     3,
       3,     3,     2,     3,     2,     1,     4,     3,     4,     6,
       8,     5,     7,     3,     5,     3,     3,     3,     3,     3,
       3,     4,     4,     4,     4,     6,     4,     4,     1,     4,
       4,     4,     4,     6,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     4,     4,     4,     4,     4,     8,
       6,     6,     6,     4,     1,     1,     1,     6,     4,     4,
       4,     4,     5,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     3,     2
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,   611,
     700,   701,     0,     0,     0,     0,     0,   620,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   632,   633,
     634,     0,   621,     0,     0,     0,     0,   692,   693,   694,
     695,   696,   697,   698,   699,     0,     0,     0,     0,   372,
     373,     0,   101,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   135,     0,     0,     0,   491,     0,
       0,     0,    46,   492,    92,     0,     0,     0,     0,     0,
       0,    21,     0,     0,   493,   494,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   133,   134,     0,     0,    39,     0,     0,   486,   488,
     487,   490,   489,     0,   472,   473,     0,     0,     0,     0,
       0,    65,     0,     0,     0,     0,     0,     0,     0,     0,
     555,   557,     0,   370,   556,   558,   559,   560,   561,   562,
     563,     0,   371,   374,   375,     0,     0,     0,     0,   578,
     655,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   364,   366,     0,     0,     0,   579,     0,     0,
       0,     0,     0,    15,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    67,   511,   510,    62,    68,     0,     0,     0,   104,
     102,   103,     0,     0,     0,     0,     0,     0,   483,   482,
     376,   379,    64,     0,   308,   307,   306,   305,   286,   287,
     288,   289,   290,   291,   292,   293,   294,   295,   296,   297,
     298,   299,   300,   301,   302,   303,   304,     0,     0,     0,
       0,     0,     0,     0,     0,   678,     0,     0,     0,     0,
       0,   702,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   714,   715,   716,     0,   703,     0,     0,     0,
       0,     0,   152,     0,   143,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   142,     0,     0,    41,   139,   155,
     138,   167,     0,    54,     0,     0,    45,     0,     0,     0,
       0,     0,     0,     0,     0,    19,     0,     0,     0,    42,
     110,   109,   106,   108,   107,    63,     0,     0,     0,     0,
       0,     0,   248,    69,     0,     0,     0,     0,     0,     0,
     474,     0,     0,     0,   255,     0,   476,   475,   367,   369,
       0,     0,     0,     0,     0,     0,     0,     0,   231,    90,
      91,    88,    89,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   559,   209,     0,     0,
       0,   166,   157,     0,     0,     0,   156,     0,    56,    55,
       0,     0,     0,   105,     0,     0,     0,   548,   551,     0,
       0,   550,   552,   549,   547,    59,    58,    30,    28,    27,
      29,    26,    25,    23,    24,     0,     0,     0,    70,     0,
       0,     0,   113,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    31,   132,     0,     0,     0,     0,
      66,   495,   496,   111,     0,   136,     0,     0,     0,     0,
     171,     0,     0,     0,     0,   170,     0,   203,     0,     0,
       0,     0,     0,   198,     0,     0,     0,     0,   112,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   140,    71,   153,    72,   168,     0,     0,     0,     0,
       0,     0,   128,     0,   377,   378,     0,    48,    50,    49,
       0,     0,     0,     0,     0,   578,   579,   654,   582,   652,
     582,     0,     0,     1,    14,     9,     8,    12,    11,    13,
      10,     0,     0,     0,     0,   365,   390,   480,     0,     0,
       0,   479,     0,   336,   477,   478,     0,     0,   481,     0,
       0,     0,     0,   334,     0,   335,     0,     0,     0,     0,
       3,     5,     4,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     7,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     6,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    77,
       0,    40,     0,   380,   383,   385,   382,   381,    43,   117,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   148,   498,   499,
     497,   149,     0,   146,   147,   554,   553,   145,     0,     0,
       0,   733,   731,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   183,     0,
     185,   181,   182,     0,     0,     0,     0,   252,   253,   250,
     251,   249,     0,     0,   259,     0,     0,     0,     0,     0,
     264,   260,   256,     0,   313,   319,   314,   320,   310,   311,
     312,   317,   315,   316,   309,   318,   261,   684,   685,   686,
     687,   688,   689,   690,   691,   471,     0,   368,    20,     0,
       0,   236,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   238,     0,   232,     0,     0,     0,     0,   215,     0,
     229,   226,   213,   211,   214,   227,   228,   210,     0,   212,
     223,   224,     0,     0,     0,   163,   162,   161,   160,   159,
       0,    61,     0,     0,     0,     0,     0,    60,   322,   323,
     321,   324,    47,    18,    16,    22,     0,   114,   115,     0,
     120,     0,     0,   118,     0,     0,    33,     0,    34,    35,
     327,   328,   325,   326,   329,   330,    32,     0,   657,     0,
      75,     0,     0,   663,     0,   580,     0,     0,   175,   180,
     177,   178,   174,   173,   176,     0,   206,   204,   207,   205,
       0,   201,   199,   202,   200,     0,   283,     0,     0,     0,
     277,   282,   281,   284,   276,   195,   194,   197,   196,     0,
     141,   154,   169,   190,   189,   192,   191,     0,   126,     0,
       0,     0,     0,   570,   574,     0,     0,     0,     0,   653,
     732,     0,     0,     0,   389,   467,     0,   456,     0,     0,
       0,     0,     0,     0,   388,     0,     0,     0,     0,   514,
       0,     0,   512,     0,     0,     0,     0,     0,     0,   513,
       0,     0,     0,     0,   386,   391,   393,   404,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   387,
     424,   426,   363,     0,     0,     0,     0,     0,   484,   485,
       0,   351,     0,   515,   516,   517,   518,   519,   520,   521,
     522,   523,   524,   527,   528,   525,   526,   529,   530,   531,
     532,   533,   534,   535,   536,   537,   538,   539,   540,   541,
     542,   543,   340,   338,   339,   354,     0,     0,     0,     0,
       0,     0,     0,     0,   341,   337,     0,   350,     0,     0,
     571,   575,     0,   651,   582,   650,   649,   648,   647,   646,
     645,   644,   584,   582,   588,   582,   592,   582,   596,   582,
     602,   582,     0,   730,   729,   728,   727,   726,   725,   724,
     723,   585,   582,   589,   582,   593,   582,   597,   582,   669,
     600,   582,     0,   604,   603,   605,   672,   606,   673,   607,
     674,     0,     0,   609,   676,   610,   677,   612,   679,   613,
     680,   614,   681,   615,   682,     0,     0,   100,   622,   704,
     623,   624,   706,   625,   707,   626,   708,     0,     0,     0,
       0,     0,     0,     0,     0,   631,   713,     0,     0,   639,
     718,   640,   719,   641,   720,   642,   721,     0,    78,   384,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   150,   151,     0,
     657,     0,     0,     0,   665,   666,   667,   668,   670,     0,
       0,     0,    84,   270,     0,     0,   184,   186,   187,     0,
      51,   254,   269,   257,   258,   262,     0,     0,   268,     0,
     263,     0,     0,     0,     0,   208,     0,     0,   235,   234,
     233,     0,     0,     0,     0,   241,   242,   239,     0,     0,
       0,   219,     0,   216,     0,   220,   221,   225,   230,   222,
     164,   165,     0,    38,     0,     0,    86,    52,    17,   121,
     116,     0,     0,   119,     0,     0,     0,     0,     0,    73,
       0,     0,     0,     0,     0,     0,   179,   172,     0,   285,
     280,   279,     0,     0,     0,   130,   129,     0,     0,   656,
       0,   580,   468,   469,   470,     0,   465,   464,   507,   506,
     505,   462,   458,   457,   466,   459,   460,     0,   405,   402,
     401,   400,   409,   403,   410,     0,     0,     0,     0,     0,
     394,   396,     0,     0,     0,     0,     0,   395,   397,   398,
     399,   504,   501,   502,   503,   500,   419,   406,   422,   420,
     421,     0,   392,   437,   432,     0,   453,   452,   431,   430,
     450,   438,   434,   436,   435,   454,   441,   429,   433,   545,
     546,   544,   442,   439,   440,     0,   443,     0,   444,   427,
     428,   449,     0,   425,   355,   358,   357,     0,     0,   359,
     356,   352,   353,   343,   345,   347,   342,   349,   348,   344,
     346,     0,     0,     0,   658,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   671,     0,   705,     0,
       0,     0,   656,   658,     0,     0,     0,     0,     0,     0,
     275,     0,   265,     0,     0,     0,    76,     0,     0,   581,
       0,   127,     0,   243,   244,   245,     0,   240,   237,   217,
     218,     0,     0,     0,     0,   331,   332,   333,     0,     0,
       0,    37,   569,    82,    80,   567,   572,   576,   509,   508,
     664,   661,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   463,     0,   411,   417,   415,   413,   407,
     412,   418,   416,   414,   408,   423,   451,   446,   445,   448,
     447,   455,   361,   360,   362,     0,     0,     0,   643,   722,
     608,   675,   616,   618,   617,   619,     0,     0,   628,   710,
     629,   711,   630,   712,   635,   637,   636,   638,    79,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      99,   266,   267,     0,    74,     0,     0,     0,     0,     0,
     247,   246,     0,     0,     0,     0,     0,    94,     0,   659,
       0,    53,   278,     0,     0,   131,     0,     0,   564,   581,
     461,    95,    96,   565,     0,     0,   683,   717,     0,   659,
       0,     0,     0,     0,     0,     0,     0,    83,    81,   573,
     577,   662,     0,     0,     0,    93,     0,     0,     0,     0,
       0,     0,     0,     0,    57,     0,   627,   709,     0,   144,
      44,    85,     0,     0,     0,     0,   660,     0,   158,    87,
       0,     0,   123,     0,   566,     0,   193,   188,     0,   660,
       0,     0,     0,     0,     0,     0,     0,     0,    36,     0,
       0,     0,     0,     0,     0,   568,    98,     0,     0,     0,
      97,     0,     0,     0,     0,     0,     0,     0,   272,     0,
       0,     0,     0,   124,     0,   273,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   271,     0,     0,     0,     0,
       0,     0,     0,   137,   274,     0,   122
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   154,   155,   416,   156,   157,   158,   159,   247,   746,
     812,   836,  1418,   160,   543,   161,   162,   163,   220,   525,
     934,   935,   959,   960,   914,   904,   164,   544,   221,   972,
     165,   166,   443,   681,  1296,  1261,  1431,   204,   937,  1002,
    1322,   406,   687,   295,   168,   169,   170,   609,  1024
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -882
static const yytype_int16 yypact[] =
{
    1635,  -363,  -275,  -271,  -267,  -152,  -136,  -107,   -76,   -28,
    -882,  -882,   -64,   -46,    -2,    25,    28,     3,    34,    41,
      47,    49,    65,    67,    75,    78,    80,    89,    61,   102,
     104,    92,   116,   117,   120,   121,   122,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,   115,   -78,   -78,  -106,  -882,
    -882,   293,  -192,   503,   -78,   446,   682,  1014,   455,    32,
     151,   157,   -90,   -78,  -882,   -45,   181,   160,  -882,   133,
    4149,   467,  -882,  -882,  -882,   502,   -78,   304,   468,    20,
       5,   461,   150,   282,  -882,  -882,   -85,  -143,   186,  2047,
    2628,   205,   206,  -150,   194,   -92,   -29,  4149,   111,   366,
     425,  -882,  -882,   527,   658,  -882,   203,   144,  -882,  -882,
    -882,  -882,  -882,   531,   -99,  -882,   227,   154,  -239,  2376,
     237,  -882,  4149,   239,   682,   271,  3484,   193,  3638,   171,
    -882,  -882,   503,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,   503,  -882,  -882,  -882,   208,   -92,   -92,   -92,  -331,
    -882,  4303,  4303,  4303,   635,   224,   225,   226,   229,   230,
     231,   233,    94,  -882,  1629,   246,   262,  -295,   234,   269,
     270,  1197,  1048,  -882,  4303,  4303,  4303,  4303,  4303,  4303,
    4303,  4303,  4303,  4303,  4303,  4303,  4401,   245,  4303,  4303,
    4303,  4303,  4303,  4303,  4303,  4303,  4303,  4303,  4303,  4303,
    4303,  -882,  -882,  -882,  -882,  -882,   275,   510,   -92,  -882,
    -882,  -882,   429,   291,   306,   -92,   313,   315,  -882,  -882,
    -882,  -882,  -882,   317,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,   614,   297,   299,
     314,   316,   318,   319,   325,  -882,   326,   327,   329,   330,
     331,  -882,   333,   334,   335,   336,   337,   339,   340,   341,
     342,   344,  -882,  -882,  -882,   361,  -882,   363,   364,   365,
     368,   338,  -882,   -57,   292,   387,   390,  -214,   371,   380,
     382,  4149,  4149,  4149,  -882,   383,  1080,  -882,  -882,  -882,
    -882,  -882,   394,  -882,   389,   449,  -882,   395,   -87,   554,
     423,  -234,   440,   448,   588,  -882,   592,   292,  1171,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,   752,   456,   -92,   457,
     458,   460,  -882,  -882,   596,   -75,   600,   298,   -78,   -92,
    -882,   602,   277,    85,  -882,   607,  -882,  -882,    94,  -882,
     469,   609,   470,   477,  -151,  -128,  4149,   478,  -882,  -882,
    -882,  -882,  -882,   420,   618,   480,   481,   482,   485,   488,
     489,   490,  -214,   493,   495,  4149,  4149,  -882,  1537,  -225,
     496,  -882,   292,   497,   498,  -214,  -882,  1707,  -882,  -882,
     -92,   612,   508,  -882,   491,   625,  1745,  -882,  -882,   513,
      -4,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,   660,   819,   515,  -882,   908,
     823,   -49,  -882,   824,    19,   -48,   664,   523,  4149,   524,
     525,   178,  4149,   397,  -882,  -882,   526,   312,  4149,   572,
    -882,  -882,  -882,  -882,  4149,  -882,   685,   547,   853,   549,
     292,   550,   551,  -214,   553,  -882,  1864,  -882,   556,   557,
     560,   561,  2119,  -882,   562,   563,   564,   587,  -882,   698,
     595,   384,   611,   613,  -214,   615,   617,   620,   623,   639,
    2163,   640,  -882,   661,  -882,   662,  4149,  4149,  4149,  4149,
    2187,   919,   -20,  4149,  -882,  -882,   760,  -882,  -882,  -882,
    4303,  4149,   642,   652,   653,   382,   383,   649,  -296,   649,
     -60,  1948,  2310,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,   256,   704,  1211,  2037,  -882,  -882,  -882,   965,   350,
     -83,  -882,  2136,  -882,  -882,  -882,  4149,   675,  -882,   691,
    2880,   682,  4149,  -882,   692,  -882,   835,   836,  4303,  4149,
    -882,  -882,  -882,  4303,  4303,  4303,  4303,  4303,  4303,  4303,
    4303,  4303,  4303,  4303,  4303,  4303,  -882,  4303,  4149,  4149,
    4149,  4149,  4149,  4149,  4149,  4149,  4303,  4303,  4303,  4303,
    4149,  4303,  -882,  4149,  2358,  2334,  6125,  2415,  6149,  2439,
    6173,  2578,  1144,  2720,  6197,  2608,  6221,  2756,  6245,  2780,
    6269,  2814,  6293,  2839,  6317,  2917,  6341,  2956,   -59,  6365,
    3140,   674,  6389,  3168,  6413,  3193,  6437,  3324,  6461,  3348,
    6485,  3504,  6509,  3659,  6533,  3683,  6557,  3718,  6581,  3839,
    6605,  3872,  6629,  4083,  6653,  4107,  6677,  4169,   837,  -882,
     697,  -882,   699,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    4149,  4149,  4149,  4149,  4149,  4149,  4149,  4149,  4149,  4149,
    4149,  4149,  4149,   700,  4149,  4149,  4149,  4149,  4149,  4149,
    4149,  4149,  4149,  4149,  4149,  4149,  4149,  -882,   703,  -882,
     705,  -882,   167,  -882,  -882,  -882,  -882,  -882,   397,   512,
    4149,   684,   684,  4202,  4149,  4149,  4149,  4149,  4149,  4149,
    4149,   707,   852,   856,  1020,   702,   858,   717,  -882,   718,
    -882,  -882,  -882,   719,   706,  4149,   722,  -882,  -882,  -882,
    -882,  -882,   723,   -56,  -882,   724,  2136,   726,   728,  4149,
    -882,  -882,  -882,   729,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,   -97,   730,  -882,  -882,   712,
     714,  -882,   734,   735,  4149,   736,   -84,  -129,  -214,   738,
     720,  -882,   530,  -882,   740,    54,   744,   745,  -882,   746,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  1052,  -882,
    6701,  6701,  4149,   748,   749,  -882,  -882,  -882,  -882,  -882,
    4149,  -882,  1006,   732,   894,   895,  4149,  -882,  -882,  -882,
    -882,  -882,  -882,  1059,  -882,  -882,  1061,  -882,  -882,  1065,
    -882,   741,  1067,  -882,  1084,   761,  -882,  4241,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  4265,   780,   764,
    -882,  1087,  4335,  -882,   783,  -334,  4359,   785,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  4149,  -882,  -882,  -882,  -882,
    4149,  -882,  -882,  -882,  -882,   786,  -882,  1096,   792,   793,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  4149,
    -882,  -882,  -882,  6701,  6701,  6701,  6701,  4149,  -882,  1099,
    1101,  4383,   779,  3549,  6725,  2540,  4303,   427,   512,  -882,
    -882,   799,   800,   820,  -882,  -882,   941,  -882,   821,   825,
     183,    99,   828,   -61,  -882,   -92,  4149,  4149,   -92,  -882,
     832,   842,  -882,   845,  3057,  3386,  4149,  4149,    83,  -882,
     846,   847,   -40,   804,  1211,  -882,  -882,  -882,   848,  1129,
     976,   849,   850,  1636,   870,   183,    -5,   871,    83,   873,
    1158,   -41,   874,   875,  3819,  3973,    -6,   877,   859,  2037,
    -882,  -882,  -882,   880,   881,   882,   101,    83,  -882,  -882,
     -50,  -882,   883,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  6701,  -882,  -882,   -78,   884,  4149,   885,
     886,   889,  4149,   890,  6701,  -882,  4421,  -882,   872,   876,
    3549,  6725,  3298,  2863,  6725,   967,  1605,  1605,  1605,  1605,
    1605,  1605,    44,  -121,    44,  -121,  -337,    16,  -337,    16,
    -337,    16,   925,   802,   673,  3040,  3040,  3040,  3040,  3040,
    3040,    44,  -121,    44,  -121,  -337,   114,  -337,   216,  -304,
    -337,   267,   763,  -882,   840,  -882,  -882,  -882,  -882,  -882,
    -882,  4303,  4149,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  4303,  4303,  -882,  -882,  -882,
     878,  -882,  -882,  -882,  -882,  -882,  -882,  4149,  4149,  4303,
    4149,  4303,  4149,  4303,  4149,  -882,  -882,  4303,  4303,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,   682,  -882,  -882,
    4445,  4469,  4493,  4517,  4541,  4565,  4589,  4613,  4637,  4661,
    4685,  4709,  4733,   887,  4757,  4781,  4805,  4829,  4853,  4877,
    4901,  4925,  4949,  4973,  4997,  5021,  5045,  -882,  -882,   879,
    -882,  -344,  5069,  5093,   309,   309,  -304,  -304,  -304,  5117,
     888,  1007,  -882,   892,  1198,    12,  -882,  -882,  -882,   896,
    6701,  -882,  -882,  -882,  -882,  -882,  2136,   893,  -882,  5141,
    -882,   898,   312,  4149,   397,  -882,  4149,   900,  -882,  -882,
    6701,   897,  2136,   903,   904,  -882,  -882,  -882,   901,  4149,
    4149,  -882,   -50,  -882,   906,  -882,  -882,  -882,  -882,  6701,
    -882,  -882,  5165,  -882,  4149,  1021,  -882,  6701,  -882,  -882,
    -882,  1063,   142,  -882,   922,  4149,  4149,  4149,   926,  -882,
    4149,  4149,  4303,    72,  4149,  4149,  -882,  6701,  5189,  -882,
    -882,  -882,   920,  5213,  5237,  1188,  -882,  4149,  4149,   940,
     932,  -344,  -882,  -882,  -882,   942,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,   934,  -882,  6701,
    6701,  -882,  -882,  -882,  -882,   954,   -92,   956,   957,   958,
    -882,  6701,   959,   -92,   961,   962,   963,  -882,  6701,  6701,
    6701,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  4149,  -882,  -882,  -882,   964,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,     1,  6701,     8,  6701,  -882,
    -882,  -882,  1270,  -882,  -882,  -882,  -882,   966,   968,  -882,
    -882,  -882,  -882,  -882,  -882,  6701,  -882,  -882,  -882,  6701,
    -882,  4149,   969,   984,   985,  4303,  4149,  6749,  5261,  6773,
    5285,  6797,  5309,  5333,  5357,  6821,  5381,  6845,  5405,  6869,
    5429,  6893,  5453,  6917,  5477,   987,  -882,  4149,  -882,  4149,
     397,  4149,  -882,  -882,  4149,   989,  1223,  1295,   972,  1298,
     977,   973,  -882,   997,  4149,   980,  -882,  1305,  5501,  -301,
    5525,  -882,  1001,  -882,  -882,  -882,  1003,  6701,  6701,  -882,
    -882,  4149,  5549,  1260,  1311,  -882,  -882,  -882,  1011,   986,
    5573,  6701,  6701,  -882,  6701,  6701,  3549,  6725,  -882,  -882,
    -882,  -882,  5597,  5621,  1009,  1010,  4149,  4149,  1317,  5645,
    5669,  4149,   397,  -882,  1023,  -882,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  6701,  -882,  -882,  -882,  -882,
    -882,  -882,  -882,  -882,  6701,   999,  1002,  4149,  3549,  6701,
    -882,  -882,  -882,  -882,  -882,  1012,  4149,  4149,  -882,  -882,
    -882,  -882,  -882,  -882,  -882,  -882,  -882,  1013,  -882,  5693,
    5717,   -72,  5741,  5765,  1015,  1008,  1018,  1334,  1019,  1336,
    -882,  -882,  6701,  1032,  -882,  4149,  4303,   287,  4149,  4149,
    -882,  -882,  5789,  4149,  1022,   142,  1024,  -882,  4149,  1038,
    4149,  -882,  -882,  5813,  5837,  -882,  1338,  4149,  6701,   -72,
    -882,  -882,  -882,  6701,  5861,  5885,  -882,  -882,  4149,  -882,
    4149,  1044,  1189,  1351,  1030,  1369,  1046,  -882,  6701,  3549,
    6725,  -882,  5909,  5933,  4149,  6701,  1209,    -8,  1068,  5957,
    4149,  5981,  4149,  4149,  -882,  6005,  -882,  -882,  6029,  6701,
    -882,  -882,  1050,  1375,  1055,  1380,  1075,  1078,  6701,  -882,
    1079,  1077,  -882,  4149,  6701,  4149,  6701,  6701,  1097,  -882,
    1100,  1083,  1408,  1085,  4149,  1086,  1088,  1106,  6701,   854,
    1089,  1095,  1115,  1102,  1421,  6701,  -882,  1118,  1324,  4149,
    -882,  1120,  1104,  1430,  1107,  1108,  1128,  6053,  -882,  1130,
    1110,  1132,  1133,  -882,  4149,  -882,  1134,  1116,  1117,  6077,
    1119,  1137,  1139,  4149,  1140,  -882,  1122,  6101,  1124,  1145,
    4149,  1160,  1366,  6701,  -882,  1165,  -882
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -882,  -882,  -882,  1149,  -882,  -882,  -882,  -882,  -102,  -882,
    -882,  -882,    51,  -882,  -882,  -882,  1488,  1504,   190,  1239,
    -882,   654,  -882,   630,  -882,  -882,  -882,  -882,   -54,  -768,
    -882,  -882,  -420,  -882,  -881,   645,    84,   -39,  -882,  -720,
    -882,  -882,  -354,    52,  -882,  -882,  -882,   212,     0
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -734
static const yytype_int16 yytable[] =
{
     172,   819,   425,   294,   705,   824,  1176,  1204,   205,  1265,
     353,   436,   968,  1181,  1173,   222,  1319,   841,   787,   765,
    1195,   324,   468,   332,   306,   344,   762,   709,   763,   355,
    1299,   799,   793,   345,   206,   377,   386,   325,   678,   766,
     767,   393,   818,   823,   202,   968,   307,   202,   334,  1192,
     304,   173,   167,   434,   444,   768,  1389,   296,   390,   441,
     441,   794,   209,  1232,  1329,   455,   500,  1316,    49,    50,
     318,  1457,  1233,   218,   565,  1381,   219,   685,  1459,   686,
      53,   437,  1233,  1182,   567,  1234,  1340,   391,   501,   378,
     387,   335,   497,   498,   499,   308,  1506,   396,   336,   853,
    1428,  1429,   548,  1390,   889,   640,   210,   699,   526,   218,
     545,   394,   219,   890,  -733,  1507,   298,   583,  1508,   456,
     873,   438,   462,  1183,   549,   583,   480,   356,   490,   769,
     969,   747,   748,   749,   750,   751,   752,   753,   754,  1193,
     211,   679,   174,   357,  1312,  1430,   175,   337,   207,   968,
     176,   508,   510,   512,   641,  1320,   680,   445,  1321,   821,
    1196,   645,   710,   969,   305,  1291,  1292,   397,   338,   442,
     442,   795,   521,  1194,   585,   587,   589,   591,   593,   595,
     597,   599,   601,   603,   605,   607,   610,   339,   613,   615,
     617,   619,   621,   623,   625,   627,   629,   631,   633,   635,
     637,  1266,  1341,   506,   506,   506,  1174,   203,   299,   970,
     203,   723,   171,   747,   748,   749,   750,   751,   752,   753,
     754,   218,  1300,  1309,   219,   395,   506,   506,   506,   506,
     506,   506,   506,   506,   506,   506,   506,   506,   506,   427,
     506,   506,   506,   506,   506,   506,   506,   506,   506,   506,
     506,   506,   506,   358,  1293,   218,  1330,  1313,   219,   346,
     347,   522,   441,  1458,   764,   177,   398,   969,   770,   399,
    1460,   309,  1258,   822,   718,   340,   310,   481,  1337,   400,
    1338,   178,   401,   730,  1314,   732,   341,   578,   579,   580,
     581,   691,   692,   693,   526,   218,   724,   300,   219,   731,
     583,   771,   218,   132,   133,   219,   342,   301,  1262,   778,
     179,   354,   311,   971,  1263,  1428,  1429,   439,   402,  1184,
     428,   429,   494,   430,  1294,   491,   359,   360,   457,   706,
     463,   495,   458,   482,   464,  1267,   801,   141,   142,   755,
     469,   180,   143,   144,  1507,   218,  1202,  1538,   219,   830,
    -731,   901,   403,   181,  1259,   734,   772,   898,   312,  1087,
     313,   583,   735,   507,   509,   511,   470,  1295,   459,   483,
     465,   182,   442,   523,   524,   790,   791,   808,   809,   810,
     811,   326,  -678,   840,  1580,  1581,   584,   586,   588,   590,
     592,   594,   596,   598,   600,   602,   604,   606,   404,   327,
     612,   614,   616,   618,   620,   622,   624,   626,   628,   630,
     632,   634,   636,  -702,  1197,   183,   460,   870,   466,   831,
     441,  1149,   736,   737,  1260,   492,   580,   581,   827,   328,
     431,   726,   837,   902,  1409,   903,   343,   583,   842,  1015,
     361,   362,   184,   471,   846,   185,   493,   484,   472,   432,
     473,   186,   563,   564,   474,   565,  1392,   218,   187,   314,
     219,   738,   739,   740,   188,   567,   189,   905,   485,   936,
     961,  -714,  1403,   747,   748,   749,   750,   751,   752,   753,
     754,   329,   190,   330,   191,   838,   883,   884,   885,   886,
     407,   845,   192,   891,   461,   193,   467,   194,   741,   218,
     894,   895,   219,   408,   409,   218,   195,   405,   219,   196,
     410,   201,  -715,   218,  -716,   774,   219,   832,   833,   834,
     835,   727,   475,   963,  -667,   581,  -703,   964,   208,   965,
     442,  1415,  1416,  1417,   197,   583,  1003,   198,   199,   200,
    1014,   223,  1016,   742,   743,   775,   297,   302,  1021,  1022,
     316,   867,   506,   303,   411,   412,   315,   868,   319,   333,
     433,  1033,  1035,  1037,  1039,  1041,   350,   351,  1043,  1044,
    1045,  1046,  1047,  1048,  1049,  1050,  1052,  1054,  1056,  1058,
    1059,  1061,   744,  1062,   728,   218,   745,   352,   219,   392,
     331,   224,   966,   212,   413,   320,   414,   776,   213,   777,
     506,   388,   389,   729,   417,   506,   506,   506,   506,   506,
     506,   506,   506,   506,   506,   506,   506,   506,   418,   506,
     426,   218,   435,   440,   219,   496,  -668,   581,   506,   506,
     506,   506,   321,   506,   214,   513,   967,   583,   514,   515,
     516,   611,   322,   517,   518,   519,   215,   520,   550,   415,
    1120,  1121,  1122,  1123,  1124,  1125,  1126,  1127,  1128,  1129,
    1130,  1131,  1132,   546,  1134,  1135,  1136,  1137,  1138,  1139,
    1140,  1141,  1142,  1143,  1144,  1145,  1146,  -670,   581,   547,
     216,  1250,   217,   551,   552,   843,   642,   643,   583,   869,
    1152,   638,   130,   131,  1153,  1154,  1155,  1156,  1157,  1158,
    1159,   419,   644,   218,   218,   649,   219,   219,   682,   646,
     225,   647,   893,   648,   650,  1170,   651,   697,   698,   580,
     699,  1203,   134,   135,   136,   137,   138,   139,   140,  1179,
     583,   652,   703,   653,   677,   654,   655,   420,   421,   224,
    1150,  1151,   656,   657,   658,   843,   659,   660,   661,   422,
     662,   663,   664,   665,   666,   323,   667,   668,   669,   670,
    1020,   671,  1397,   224,  1190,  1023,  1025,  1026,  1027,  1028,
    1029,  1030,  1031,  1032,  1034,  1036,  1038,  1040,   672,  1042,
     673,   674,   675,   683,   407,   676,   684,   688,  1051,  1053,
    1055,  1057,  1209,  1060,   906,   907,   689,   408,   409,   908,
    1212,   690,   694,  1200,   410,   702,  1217,   130,   131,   701,
     704,   707,   423,   226,   227,   228,   229,   230,   231,   708,
     232,   233,   234,   235,   236,   237,   238,   239,   240,   241,
     242,   243,   244,   245,   246,   909,   711,   134,   135,   136,
     137,   138,   139,   140,   712,   713,   714,   716,   411,   412,
     844,   722,   717,   719,   720,  1237,   721,   725,   225,   733,
    1238,  1268,   756,   759,  1271,   758,   760,   130,   131,   910,
    1280,  1287,   911,   761,   773,   779,   780,   781,   782,  1243,
     936,   783,   225,   912,   784,   785,   786,  1244,   413,   788,
     414,   789,   796,   797,   798,   802,   610,   134,   135,   136,
     137,   138,   139,   140,   803,   961,   639,   804,   805,   807,
     814,   815,  1339,   816,   817,   820,  1269,  1270,   825,   826,
     828,   829,   839,   913,  1281,  1288,  1289,  1290,   568,   569,
     570,   571,   572,   573,   574,   575,   695,   696,   697,   698,
     580,   699,   847,   848,   849,   850,   851,   852,   506,   854,
    1251,   583,   856,   857,  1326,  1328,   858,   859,   861,   862,
     863,   226,   227,   228,   229,   230,   231,  1343,   232,   233,
     234,   235,   236,   237,   238,   239,   240,   241,   242,   243,
     244,   245,   246,   864,   865,   226,   227,   228,   229,   230,
     231,   866,   232,   233,   234,   235,   236,   237,   238,   239,
     240,   241,   242,   243,   244,   245,   246,   871,  1345,   872,
     888,   874,  1349,   875,   892,  1375,   876,   248,   249,   877,
     250,   251,   252,   253,   254,   255,    10,    11,   256,   257,
     258,   259,   260,   261,   262,   878,   880,   263,   264,   265,
     266,   267,   268,   269,   270,   271,   424,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   962,   881,   882,   896,
      37,    38,    39,    40,    41,    42,    43,    44,   897,   898,
     567,  1004,  1358,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,  1360,  1362,  1005,  1017,  1018,
    1019,  1117,  1090,  1118,   583,  1119,  1133,  1363,  1364,  1147,
    1366,  1148,  1368,  1160,  1370,   583,  1161,  1372,  1374,   281,
    1162,  1163,  1165,  1166,  1167,  1168,   282,  1164,  1171,  1172,
    1175,  1169,  1177,   506,  1178,  1180,  1185,  1186,  1396,  1187,
    1188,  1189,  1191,  1609,  1198,  1199,  1201,   506,   506,   283,
    1205,  1206,  1207,  1208,  1210,  1211,  1213,  1214,  1215,  1216,
    1218,   506,  1219,   506,   284,   506,  1220,  1221,  1222,   506,
     506,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,  1223,  1224,  1227,  1229,  1228,
    1231,  1236,  1239,  1398,   583,  1356,  1400,  1240,  1241,  1242,
    1245,   285,  1246,   286,  1248,  1252,  1253,   287,  1255,  1407,
    1408,   569,   570,   571,   572,   573,   574,   575,   695,   696,
     697,   698,   580,   699,  1412,   218,  1254,  1256,   219,  1301,
    1304,  1257,  1446,   583,  1264,  1420,  1421,  1422,  1272,  1451,
    1424,  1425,  1427,  1305,  1432,  1433,  1399,   915,  1273,   916,
     917,  1274,  1297,  1298,  1303,  1306,  1307,  1439,  1440,  1318,
    -671,   288,   568,   569,   570,   571,   572,   573,   574,   575,
     695,   696,   697,   698,   580,   699,  1310,  1315,   289,  1317,
    1323,  1324,   918,  1331,  1332,   583,  1334,  1335,  1336,  1342,
    1344,  1346,  1347,  1357,   506,  1348,  1350,  1352,  -705,  1388,
    1386,  1353,  1391,   919,  1395,  1380,  1401,  1359,  1361,  1404,
    1405,  1455,  1410,  1385,  1413,  1378,   920,  1387,  1393,   130,
     131,  1365,  1402,  1367,   921,  1369,  1406,  1414,  1419,  1371,
    1373,  1438,  1423,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,  1435,   565,  1441,  1443,   134,
     135,   136,   137,   138,   139,   140,   567,  1355,  1442,  1444,
    1445,  1464,  1447,  1448,  1449,  1450,  1469,  1452,  1453,  1454,
    1456,  1461,  1462,  1495,  1463,  1465,   922,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,  1489,   565,  1490,
    1466,  1492,  1467,  1488,  1493,  1494,  1496,  1497,   567,  1498,
     923,  1500,  1499,  1501,  1502,  1503,  1504,  1510,   924,  1511,
    1514,   925,  1515,  1516,  1517,  1521,  1522,   506,  1525,   290,
     150,  1512,   218,   926,   927,   219,   928,  1531,   929,  1530,
    1532,   291,  -683,  -717,  1542,  1544,   292,  1546,  1547,  1564,
    1541,   293,  1491,  1543,  1545,  1560,  1523,  1524,  1556,  1558,
    1570,  1528,  1572,  1571,  1426,  1573,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
    1574,  1575,   582,  1579,  1582,  1590,  1591,  1533,   930,   583,
    1592,  1593,  1594,   931,  1595,  1596,  1534,  1535,   568,   569,
     570,   571,   572,   573,   574,   575,   695,   696,   697,   698,
     580,   699,  1597,  1600,  1529,   700,  1601,   932,  1602,  1603,
    1604,   583,  1608,  1607,  1606,  1548,  1550,  1610,  1552,  1553,
    1611,  1612,  1614,  1555,  1615,  1616,  1618,  1613,  1559,  1619,
    1561,  1620,  1621,  1622,  1623,  1626,  1625,  1565,  1627,  1628,
    1630,  1631,  1632,  1635,  1634,  1636,  1638,  1639,  1568,  1641,
    1569,  1642,   553,   554,   555,   556,   557,   558,   559,   560,
     561,   562,   563,   564,  1578,   565,  1644,  1645,   506,  1071,
    1584,  1646,  1586,  1587,   813,   567,  1557,  1468,   348,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,  1598,   349,  1599,   715,   757,  1302,  1333,
    1311,  1551,   583,     0,  1605,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,   933,   565,  1617,
       0,   566,     0,     0,     0,     0,     0,     0,   567,     0,
       0,     0,     0,     0,  1629,     0,     0,     0,     0,     0,
       0,     0,     0,  1637,     0,    -2,     1,     0,     2,     3,
    1643,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,     0,     0,    19,    20,
      21,    22,    23,    24,    25,    26,    27,     0,    28,    29,
      30,    31,    32,    33,    34,    35,    36,     0,     0,     0,
       0,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,     0,     0,    47,    48,     0,     0,     0,    49,    50,
       0,     0,     0,     0,     0,    51,    52,     0,     0,     0,
      53,    54,    55,     0,     0,    56,     0,     0,  1549,    57,
       0,     0,    58,     0,   527,   528,     0,    59,     0,    60,
      61,     0,    62,    63,    64,    65,     0,   973,    66,     0,
      67,    68,    69,     0,    70,     0,     0,     0,   529,     0,
      71,     0,     0,     0,   530,    72,     0,   974,    73,     0,
     531,     0,   532,    74,    75,     0,     0,     0,     0,    76,
      77,    78,     0,     0,     0,    79,    80,     0,     0,    81,
       0,   533,     0,     0,     0,    82,     0,     0,     0,    83,
       0,     0,     0,     0,     0,     0,    84,    85,    86,    87,
       0,    88,     0,     0,     0,     0,   534,    89,   535,     0,
       0,    90,    91,     0,    92,    93,    94,    95,     0,     0,
     536,     0,     0,     0,     0,     0,     0,     0,    96,    97,
     218,     0,     0,   219,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    98,     0,     0,     0,     0,     0,     0,
       0,     0,   537,     0,     0,     0,     0,     0,   975,     0,
      99,   100,     0,     0,     0,     0,   101,   102,   103,   104,
     105,   106,   107,     0,     0,     0,     0,     0,     0,   108,
     109,   110,   111,   112,     0,   113,   538,   539,     0,   114,
     115,     0,     0,     0,   116,     0,   117,     0,     0,     0,
     118,     0,     0,     0,     0,     0,     0,   540,     0,     0,
     119,   120,   121,   122,     0,   541,     0,   123,     0,     0,
       0,   124,     0,   125,     0,     0,   126,   127,   128,   129,
     130,   131,     0,   132,   133,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,   542,   792,     0,     0,     0,     0,     0,   583,     0,
     134,   135,   136,   137,   138,   139,   140,   141,   142,     0,
       0,     0,   143,   144,     0,     0,     0,     0,     0,   976,
     977,   978,   979,   980,   981,   982,   983,   984,   985,   986,
     987,   988,   989,   990,   991,   992,   993,   994,   995,   996,
     997,   998,   999,  1000,  1001,  -734,  -734,  -734,  -734,  -734,
    -734,   561,   562,   563,   564,   145,   565,     0,     0,     0,
     146,   147,   148,     0,     0,     0,   567,     0,     0,     0,
     149,   150,  1308,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   151,     0,     0,     0,     0,   152,     0,     0,
     248,   249,   153,   250,   251,   252,   253,   254,   255,    10,
      11,   256,   257,   258,   259,   260,   261,   262,     0,     0,
     263,   264,   265,   266,   267,   268,   269,   270,   271,     0,
     272,   273,   274,   275,   276,   277,   278,   279,   280,     0,
       0,     0,     0,    37,    38,    39,    40,    41,    42,    43,
      44,     0,   938,     0,   939,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,   800,     0,     0,     0,     0,   940,   583,     0,
       0,   363,   941,     0,     0,     0,     0,   364,     0,     0,
       0,     0,   365,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
     806,     0,     0,     0,     0,     0,   583,     0,   942,     0,
     943,     0,     0,     0,     0,     0,     0,     0,   366,     0,
       0,     0,     0,     0,   944,     0,     0,   317,     0,     0,
       0,     0,     0,     0,   367,     0,     0,     0,     0,     0,
       0,     0,   945,     0,     0,   946,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   368,   947,     0,     0,     0,
     369,     0,     0,     0,   370,     0,   371,     0,     0,     0,
     372,     0,     0,     0,     0,     0,     0,   973,   218,     0,
       0,   219,   948,     0,     0,     0,     0,     0,   218,     0,
       0,   219,     0,     0,     0,     0,     0,   974,     0,     0,
     949,   950,   568,   569,   570,   571,   572,   573,   574,   575,
     695,   696,   697,   698,   580,   699,     0,     0,     0,   855,
       0,     0,     0,     0,   288,   583,     0,     0,     0,     0,
       0,     0,     0,   951,     0,   952,     0,     0,     0,     0,
       0,   289,     0,     0,     0,     0,   953,   954,     0,     0,
       0,   955,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   373,   956,     0,     0,     0,   957,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   374,     0,     0,
       0,     0,   130,   131,   375,     0,   553,   554,   555,   556,
     557,   558,   559,   560,   561,   562,   563,   564,   975,   565,
       0,     0,     0,     0,     0,     0,   899,     0,     0,   567,
       0,     0,   134,   135,   376,   137,   138,   139,   140,   248,
     249,     0,   250,   251,   252,   253,   254,   255,    10,    11,
     256,   257,   258,   259,   260,   261,   262,     0,     0,   263,
     264,   265,   266,   267,   268,   269,   270,   271,     0,   272,
     273,   274,   275,   276,   277,   278,   279,   280,     0,     0,
       0,     0,    37,    38,    39,    40,    41,    42,    43,    44,
       0,     0,     0,   958,     0,     0,     0,     0,     0,     0,
       0,     0,   290,   150,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   291,     0,     0,     0,     0,   292,
       0,     0,     0,     0,   293,     0,   446,     0,     0,     0,
       0,   447,     0,     0,     0,     0,     0,     0,   448,   976,
     977,   978,   979,   980,   981,   982,   983,   984,   985,   986,
     987,   988,   989,   990,   991,   992,   993,   994,   995,   996,
     997,   998,   999,  1000,  1001,     0,     0,   449,     0,     0,
       0,     0,     0,     0,     0,     0,   450,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,   860,     0,     0,     0,     0,     0,
     583,   451,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   452,     0,     0,     0,   453,
       0,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,   218,   879,     0,
     219,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,   887,     0,     0,     0,     0,     0,   583,     0,
       0,     0,     0,   288,     0,     0,     0,     0,     0,   454,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     289,   248,   249,     0,   250,   251,   252,   253,   254,   255,
      10,    11,   256,   257,   258,   259,   260,   261,   262,     0,
       0,   263,   264,   265,   266,   267,   268,   269,   270,   271,
       0,   272,   273,   274,   275,   276,   277,   278,   279,   280,
       0,   130,   131,     0,    37,    38,    39,    40,    41,    42,
      43,    44,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   379,     0,     0,
       0,   134,   135,   136,   137,   138,   139,   140,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,     0,   380,     0,     0,     0,     0,   900,     0,
     381,   583,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,   578,   579,   580,   581,     0,     0,     0,     0,
       0,     0,  1064,     0,     0,   583,   553,   554,   555,   556,
     557,   558,   559,   560,   561,   562,   563,   564,   382,   565,
       0,   290,   150,     0,     0,     0,  1063,     0,     0,   567,
       0,     0,     0,   291,     0,     0,     0,     0,   292,     0,
       0,     0,     0,   293,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   383,     0,   384,     0,     0,
       0,   385,     0,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,     0,     0,   218,
       0,     0,   219,  1066,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   578,   579,   580,
     581,     0,     0,     0,     0,     0,     0,  1068,     0,     0,
     583,     0,     0,     0,     0,   288,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   289,   248,   249,     0,   250,   251,   252,   253,
     254,   255,    10,    11,   256,   257,   258,   259,   260,   261,
     262,     0,     0,   263,   264,   265,   266,   267,   268,   269,
     270,   271,     0,   272,   273,   274,   275,   276,   277,   278,
     279,   280,     0,   130,   131,     0,    37,    38,    39,    40,
      41,    42,    43,    44,     0,     0,     0,     0,   568,   569,
     570,   571,   572,   573,   574,   575,   695,   696,   697,   698,
     580,   699,     0,   134,   135,   136,   137,   138,   139,   140,
    1249,   583,     0,     0,     0,     0,     0,     0,     0,  1006,
    1007,     0,     0,     0,     0,  1008,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
       0,     0,     0,     0,     0,     0,  1070,     0,     0,   583,
       0,     0,     0,     0,     0,  1009,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
     317,     0,     0,   290,   150,     0,  1074,     0,     0,   583,
       0,     0,     0,     0,     0,   291,     0,     0,     0,     0,
     292,     0,     0,     0,     0,   293,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1010,     0,  1011,
     248,   249,     0,   250,   251,   252,   253,   254,   255,    10,
      11,   256,   257,   258,   259,   260,   261,   262,     0,     0,
     263,   264,   265,   266,   267,   268,   269,   270,   271,     0,
     272,   273,   274,   275,   276,   277,   278,   279,   280,     0,
       0,     0,     0,    37,    38,    39,    40,    41,    42,    43,
      44,     0,     0,     0,     0,     0,     0,   288,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,     0,     0,   289,  1072,     0,  1012,  1013,     0,
       0,   583,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1275,     0,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,   578,   579,   580,   581,     0,     0,
       0,     0,     0,     0,  1076,   130,   131,   583,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,     0,     0,     0,     0,     0,   317,  1078,     0,
    1276,   583,     0,     0,     0,   134,   135,   136,   137,   138,
     139,   140,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,   578,   579,   580,   581,     0,     0,     0,     0,
       0,     0,  1080,     0,  1277,   583,  1278,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   578,   579,   580,
     581,     0,     0,     0,     0,     0,     0,  1082,   218,     0,
     583,   219,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,   290,   150,     0,     0,     0,
       0,     0,     0,     0,   567,     0,     0,   291,     0,     0,
       0,     0,   292,     0,   288,     0,     0,   293,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   289,     0,     0,  1279,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,   578,   579,   580,   581,     0,
       0,     0,     0,     0,     0,  1084,     0,     0,   583,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   130,   131,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,   578,   579,   580,   581,     0,     0,
       0,  1086,     0,     0,     0,     0,     0,   583,     0,     0,
       0,     0,   134,   135,   136,   137,   138,   139,   140,   248,
     249,     0,   250,   251,   252,   253,   254,   255,    10,    11,
     256,   257,   258,   259,   260,   261,   262,     0,     0,   263,
     264,   265,   266,   267,   268,   269,   270,   271,     0,   272,
     273,   274,   275,   276,   277,   278,   279,   280,     0,     0,
       0,     0,    37,    38,    39,    40,    41,    42,    43,    44,
    -734,  -734,  -734,  -734,  -734,  -734,   695,   696,   697,   698,
     580,   699,   290,   150,     0,     0,     0,     0,     0,     0,
       0,   583,     0,     0,   291,     0,     0,     0,     0,   292,
       0,     0,     0,     0,   293,     0,     0,     0,     0,     0,
       0,  1282,     0,     0,     0,     0,     0,   248,   249,     0,
     250,   251,   252,   253,   254,   255,    10,    11,   256,   257,
     258,   259,   260,   261,   262,     0,     0,   263,   264,   265,
     266,   267,   268,   269,   270,   271,     0,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   317,     0,     0,  1283,
      37,    38,    39,    40,    41,    42,    43,    44,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,     0,     0,     0,     0,     0,     0,  1089,     0,
       0,   583,     0,  1284,     0,  1285,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
       0,     0,     0,     0,     0,     0,  1092,   218,     0,   583,
     219,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,   578,   579,   580,   581,     0,     0,     0,     0,     0,
       0,  1094,     0,     0,   583,     0,     0,     0,     0,     0,
       0,     0,     0,   288,   317,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     289,   248,   249,  1286,   250,   251,   252,   253,   254,   255,
      10,    11,   256,   257,   258,   259,   260,   261,   262,     0,
       0,   263,   264,   265,   266,   267,   268,   269,   270,   271,
       0,   272,   273,   274,   275,   276,   277,   278,   279,   280,
       0,   130,   131,     0,    37,    38,    39,    40,    41,    42,
      43,    44,     0,     0,     0,     0,   568,   569,   570,   571,
     572,   573,   574,   575,   695,   696,   697,   698,   580,   699,
       0,   134,   135,   136,   137,   138,   139,   140,  1354,   583,
       0,   288,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,   578,   579,   580,   581,     0,     0,   289,     0,
       0,     0,  1096,     0,     0,   583,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
       0,     0,     0,  1098,     0,     0,     0,     0,     0,   583,
       0,     0,     0,     0,     0,     0,     0,     0,   317,   130,
     131,   290,   150,     0,     0,   476,   477,     0,     0,     0,
       0,     0,     0,   291,     0,     0,     0,     0,   292,     0,
       0,     0,     0,   293,     0,     0,     0,     0,     0,   134,
     135,   136,   137,   138,   139,   140,     0,     0,     0,   478,
     479,     0,   248,   249,     0,   250,   251,   252,   253,   254,
     255,    10,    11,   256,   257,   258,   259,   260,   261,   262,
       0,     0,   263,   264,   265,   266,   267,   268,   269,   270,
     271,     0,   272,   273,   274,   275,   276,   277,   278,   279,
     280,     0,     0,     0,     0,    37,    38,    39,    40,    41,
      42,    43,    44,     0,     0,   288,     0,     0,     0,   290,
     150,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   291,   289,     0,     0,     0,   292,     0,     0,     0,
       0,   293,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,   578,   579,   580,   581,     0,     0,     0,  1100,
       0,     0,     0,     0,     0,   583,     0,     0,     0,     0,
       0,     0,     0,   130,   131,     0,     0,     0,     0,   486,
     487,     0,     0,     0,     0,     0,     0,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,   317,
     565,     0,     0,   134,   135,   136,   137,   138,   139,   140,
     567,     0,     0,   488,   489,     0,   248,   249,     0,   250,
     251,   252,   253,   254,   255,    10,    11,   256,   257,   258,
     259,   260,   261,   262,     0,     0,   263,   264,   265,   266,
     267,   268,   269,   270,   271,     0,   272,   273,   274,   275,
     276,   277,   278,   279,   280,     0,     0,     0,     0,    37,
      38,    39,    40,    41,    42,    43,    44,     0,     0,     0,
       0,     0,     0,   290,   150,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   291,     0,     0,     0,     0,
     292,     0,     0,     0,     0,   293,   288,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   578,   579,   580,
     581,     0,     0,   289,  1102,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,   578,   579,   580,   581,     0,     0,     0,  1104,     0,
       0,     0,     0,     0,   583,  1325,     0,     0,     0,     0,
       0,     0,     0,   317,   130,   131,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
       0,     0,     0,     0,     0,     0,  1106,     0,     0,   583,
       0,     0,     0,     0,   134,   135,   136,   137,   138,   139,
     140,     0,   248,   249,     0,   250,   251,   252,   253,   254,
     255,    10,    11,   256,   257,   258,   259,   260,   261,   262,
       0,     0,   263,   264,   265,   266,   267,   268,   269,   270,
     271,     0,   272,   273,   274,   275,   276,   277,   278,   279,
     280,     0,     0,     0,     0,    37,    38,    39,    40,    41,
      42,    43,    44,     0,     0,     0,     0,     0,     0,     0,
     288,     0,     0,     0,   290,   150,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   291,   289,     0,     0,
       0,   292,     0,     0,     0,     0,   293,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   578,   579,   580,
     581,     0,     0,     0,  1108,     0,     0,     0,     0,  1327,
     583,     0,     0,     0,     0,     0,     0,     0,   130,   131,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
     578,   579,   580,   581,     0,     0,     0,     0,     0,   317,
    1110,     0,     0,   583,     0,     0,     0,     0,   134,   135,
     136,   137,   138,   139,   140,     0,     2,     3,     0,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,   502,     0,     0,    19,    20,    21,    22,
      23,    24,    25,    26,    27,     0,    28,    29,    30,    31,
      32,    33,    34,    35,    36,     0,     0,     0,     0,    37,
      38,    39,    40,    41,    42,    43,    44,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   290,   150,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     291,     0,     0,     0,     0,   292,   288,     0,     0,     0,
     293,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   289,     2,     3,     0,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,   502,     0,     0,    19,    20,    21,    22,    23,    24,
      25,    26,    27,     0,    28,    29,    30,    31,    32,    33,
      34,    35,    36,   503,   130,   131,     0,    37,    38,    39,
      40,    41,    42,    43,    44,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   134,   135,   136,   137,   138,   139,
     140,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,   578,   579,   580,   581,     0,     0,     0,     0,     0,
       0,  1112,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,   578,   579,   580,   581,     0,
       0,     0,     0,     0,     0,  1114,     0,     0,   583,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     288,   503,     0,     0,   290,   150,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   291,   504,     0,     0,
       0,   292,     0,     0,     0,     0,   293,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   578,   579,   580,
     581,     0,     0,     0,     0,     0,     0,  1116,     0,     0,
     583,     0,     0,     0,     0,     0,     0,     0,   130,   131,
     568,   569,   570,   571,   572,   573,   574,   575,   695,   696,
     697,   698,   580,   699,     0,     0,     0,     0,     0,     0,
     900,     0,     0,   583,     0,     0,     0,     0,   134,   135,
     136,   137,   138,   139,   140,     0,     0,     0,   288,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,   608,  1225,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1226,     0,     0,     0,     0,     0,   583,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   130,   131,   505,   150,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     151,     0,     0,     0,     0,   152,     0,     0,     0,     0,
     153,     0,     0,     0,     0,     0,   134,   135,   136,   137,
     138,   139,   140,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1230,     0,     0,     0,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1235,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,   505,   150,  1247,     0,
       0,     0,     0,     0,   583,     0,     0,     0,   151,     0,
       0,     0,     0,   152,     0,     0,     0,     0,   153,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1351,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1376,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,     0,     0,     0,  1066,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,     0,     0,
       0,  1068,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,     0,     0,     0,  1070,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1072,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1074,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,     0,     0,     0,  1076,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,     0,     0,
       0,  1078,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,     0,     0,     0,  1080,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,  1082,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1084,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1377,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,     0,     0,
       0,  1089,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,     0,     0,     0,  1092,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,  1094,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1096,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1098,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1100,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1102,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1104,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1106,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1379,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,     0,     0,
       0,  1110,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,     0,     0,     0,  1112,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,  1114,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1116,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,     0,     0,     0,     0,     0,  1382,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,     0,     0,
       0,     0,     0,  1383,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1384,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1394,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1411,     0,     0,     0,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1434,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1436,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1437,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,  1471,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,     0,     0,     0,
       0,     0,     0,  1473,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   578,   579,   580,
     581,     0,     0,     0,     0,     0,     0,  1475,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1476,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1477,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,  1479,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1481,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,     0,     0,     0,  1483,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,   578,   579,   580,   581,     0,     0,     0,     0,     0,
       0,  1485,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,   578,   579,   580,   581,     0,
       0,     0,     0,     0,     0,  1487,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1505,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1509,     0,     0,     0,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1513,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1518,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,     0,     0,     0,     0,     0,  1519,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1520,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1526,     0,     0,     0,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1527,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,     0,     0,
       0,  1536,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,     0,     0,     0,  1537,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,     0,
       0,  1539,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1540,     0,     0,     0,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,  1554,     0,     0,     0,     0,     0,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1562,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1563,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,  1566,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
       0,     0,     0,  1567,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,     0,     0,     0,     0,     0,  1576,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1577,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1583,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1585,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   695,   696,   697,   698,   580,   699,     0,     0,     0,
    1588,     0,     0,     0,     0,     0,   583,   568,   569,   570,
     571,   572,   573,   574,   575,   695,   696,   697,   698,   580,
     699,     0,     0,     0,     0,     0,     0,     0,     0,  1589,
     583,   568,   569,   570,   571,   572,   573,   574,   575,   695,
     696,   697,   698,   580,   699,     0,     0,     0,  1624,     0,
       0,     0,     0,     0,   583,   568,   569,   570,   571,   572,
     573,   574,   575,   695,   696,   697,   698,   580,   699,     0,
       0,     0,  1633,     0,     0,     0,     0,     0,   583,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,  1640,     0,     0,     0,
       0,     0,   583,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,     0,   565,     0,     0,     0,
       0,     0,     0,  1065,     0,     0,   567,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,     0,     0,     0,  1067,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,     0,     0,
       0,  1069,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1073,     0,     0,   567,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,     0,   565,     0,     0,     0,     0,     0,     0,  1075,
       0,     0,   567,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,     0,   565,     0,     0,     0,
       0,     0,     0,  1077,     0,     0,   567,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,     0,     0,     0,  1079,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,     0,     0,
       0,  1081,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1083,     0,     0,   567,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,     0,   565,     0,     0,     0,  1085,     0,     0,     0,
       0,     0,   567,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,     0,   565,     0,     0,     0,
       0,     0,     0,  1088,     0,     0,   567,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,     0,     0,     0,  1091,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,     0,     0,
       0,  1093,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1095,     0,     0,   567,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,     0,   565,     0,     0,     0,  1097,     0,     0,     0,
       0,     0,   567,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,     0,   565,     0,     0,     0,
    1099,     0,     0,     0,     0,     0,   567,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,  1101,     0,     0,     0,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,  1103,     0,
       0,     0,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1105,     0,     0,   567,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,     0,   565,     0,     0,     0,  1107,     0,     0,     0,
       0,     0,   567,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,     0,   565,     0,     0,     0,
       0,     0,     0,  1109,     0,     0,   567,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,     0,     0,     0,  1111,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,     0,     0,
       0,  1113,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1115,     0,     0,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   695,   696,   697,
     698,   580,   699,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   583,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   583,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,     0,     0,     0,  1470,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,     0,     0,
       0,  1472,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1474,     0,     0,   567,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,     0,   565,     0,     0,     0,     0,     0,     0,  1478,
       0,     0,   567,   553,   554,   555,   556,   557,   558,   559,
     560,   561,   562,   563,   564,     0,   565,     0,     0,     0,
       0,     0,     0,  1480,     0,     0,   567,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   564,     0,
     565,     0,     0,     0,     0,     0,     0,  1482,     0,     0,
     567,   553,   554,   555,   556,   557,   558,   559,   560,   561,
     562,   563,   564,     0,   565,     0,     0,     0,     0,     0,
       0,  1484,     0,     0,   567,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,     0,   565,     0,
       0,     0,     0,     0,     0,  1486,     0,     0,   567
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-882)))

#define yytable_value_is_error(Yytable_value) \
  (!!((Yytable_value) == (-734)))

static const yytype_int16 yycheck[] =
{
       0,   421,   104,    57,    91,   425,   726,   775,    47,    70,
      95,   110,    95,   110,    70,    54,    57,   437,   372,   147,
     149,    75,   124,    77,    63,    79,   177,   261,   179,   172,
      70,   385,   257,    28,   140,    89,    90,    76,    95,   167,
     168,    95,    91,    91,   122,    95,    91,   122,    28,   133,
     140,   414,     0,   107,   293,   183,    44,    57,   208,   108,
     108,   286,   254,   397,    70,   119,   397,   948,    63,    64,
      70,    70,   416,   201,   411,   419,   204,   291,    70,   293,
      75,   180,   416,   180,   421,   419,   967,   237,   419,    89,
      90,    71,   146,   147,   148,   140,   397,    97,    78,   453,
      28,    29,   397,    91,   124,   207,   298,   411,   162,   201,
     164,   140,   204,   133,   410,   416,    84,   421,   419,   119,
     474,   220,   122,   220,   419,   421,   126,   270,   128,   257,
     213,    46,    47,    48,    49,    50,    51,    52,    53,   223,
     332,   198,   417,   286,   149,    73,   417,   127,   254,    95,
     417,   151,   152,   153,   208,   196,   213,   396,   199,   140,
     289,   215,   396,   213,   254,    82,    83,    56,   148,   218,
     218,   396,    78,   257,   174,   175,   176,   177,   178,   179,
     180,   181,   182,   183,   184,   185,   186,   167,   188,   189,
     190,   191,   192,   193,   194,   195,   196,   197,   198,   199,
     200,   262,   970,   151,   152,   153,   262,   285,   176,   292,
     285,   286,     0,    46,    47,    48,    49,    50,    51,    52,
      53,   201,   262,   943,   204,   254,   174,   175,   176,   177,
     178,   179,   180,   181,   182,   183,   184,   185,   186,    95,
     188,   189,   190,   191,   192,   193,   194,   195,   196,   197,
     198,   199,   200,   396,   171,   201,   262,   262,   204,   254,
     255,   167,   108,   262,   415,   417,   155,   213,   396,   158,
     262,    90,    89,   254,   328,   255,    95,    84,   177,   168,
     179,   417,   171,   337,   289,   339,   266,   408,   409,   410,
     411,   291,   292,   293,   348,   201,   335,   265,   204,   338,
     421,   355,   201,   298,   299,   204,   286,   275,   209,   363,
     417,   396,   131,   396,   215,    28,    29,   416,   207,   416,
     176,   177,   132,   179,   241,   154,   140,   141,    91,   416,
      91,   141,    95,   140,    95,   396,   390,   332,   333,   254,
      69,   417,   337,   338,   416,   201,   292,   419,   204,   171,
     410,    95,   241,   417,   171,    78,   356,   416,   177,   418,
     179,   421,    85,   151,   152,   153,    95,   284,   131,   176,
     131,   417,   218,   279,   280,   375,   376,   381,   382,   383,
     384,    77,   410,   437,   392,   393,   174,   175,   176,   177,
     178,   179,   180,   181,   182,   183,   184,   185,   287,    95,
     188,   189,   190,   191,   192,   193,   194,   195,   196,   197,
     198,   199,   200,   410,   768,   417,   179,   471,   179,   241,
     108,   254,   145,   146,   241,   254,   410,   411,   428,   125,
     286,   133,   432,   177,  1202,   179,   416,   421,   438,   541,
     254,   255,   417,   172,   444,   417,   275,   254,   177,   305,
     179,   417,   408,   409,   183,   411,  1176,   201,   417,   278,
     204,   184,   185,   186,   417,   421,   417,   521,   275,   523,
     524,   410,  1192,    46,    47,    48,    49,    50,    51,    52,
      53,   177,   417,   179,   417,   433,   486,   487,   488,   489,
     124,   439,   417,   493,   257,   417,   257,   417,   221,   201,
     500,   501,   204,   137,   138,   201,   417,   396,   204,   417,
     144,   396,   410,   201,   410,    95,   204,   339,   340,   341,
     342,   223,   251,   173,   410,   411,   410,   177,   235,   179,
     218,   389,   390,   391,   417,   421,   536,   417,   417,   417,
     540,    95,   542,   266,   267,   125,    91,   396,   548,   549,
     417,   167,   500,   396,   188,   189,   396,   173,    91,    91,
     416,   561,   562,   563,   564,   565,   105,   417,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,   305,   583,   286,   201,   309,   305,   204,   395,
     286,    81,   242,    90,   228,    93,   230,   177,    95,   179,
     548,   396,   396,   305,   179,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,   565,    91,   567,
     417,   201,    91,   396,   204,   417,   410,   411,   576,   577,
     578,   579,   130,   581,   131,     0,   286,   421,   414,   414,
     414,   396,   140,   414,   414,   414,   143,   414,   414,   283,
     650,   651,   652,   653,   654,   655,   656,   657,   658,   659,
     660,   661,   662,   417,   664,   665,   666,   667,   668,   669,
     670,   671,   672,   673,   674,   675,   676,   410,   411,   417,
     177,   254,   179,   414,   414,   173,   257,   396,   421,   305,
     690,   416,   295,   296,   694,   695,   696,   697,   698,   699,
     700,    43,   396,   201,   201,    91,   204,   204,   416,   396,
     200,   396,   500,   396,   417,   715,   417,   408,   409,   410,
     411,   775,   325,   326,   327,   328,   329,   330,   331,   729,
     421,   417,   283,   417,   396,   417,   417,    79,    80,    81,
     688,   689,   417,   417,   417,   173,   417,   417,   417,    91,
     417,   417,   417,   417,   417,   253,   417,   417,   417,   417,
     548,   417,  1182,    81,   764,   553,   554,   555,   556,   557,
     558,   559,   560,   561,   562,   563,   564,   565,   417,   567,
     417,   417,   417,   396,   124,   417,   396,   416,   576,   577,
     578,   579,   792,   581,    90,    91,   416,   137,   138,    95,
     800,   419,   419,   273,   144,   416,   806,   295,   296,   415,
     415,   257,   154,   303,   304,   305,   306,   307,   308,   396,
     310,   311,   312,   313,   314,   315,   316,   317,   318,   319,
     320,   321,   322,   323,   324,   131,   396,   325,   326,   327,
     328,   329,   330,   331,   396,   257,   254,    95,   188,   189,
     278,   255,   396,   396,   396,   855,   396,   257,   200,   257,
     860,   915,   255,   254,   918,   396,   396,   295,   296,   165,
     924,   925,   168,   396,   396,   257,   396,   396,   396,   879,
     934,   396,   200,   179,   396,   396,   396,   887,   228,   396,
     230,   396,   396,   396,   396,   283,   896,   325,   326,   327,
     328,   329,   330,   331,   396,   959,   396,   416,   283,   396,
      91,   396,   966,     5,    91,    91,   916,   917,   254,   396,
     396,   396,   396,   219,   924,   925,   926,   927,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,   257,   396,    91,   396,   396,   396,   896,   396,
     898,   421,   396,   396,   954,   955,   396,   396,   396,   396,
     396,   303,   304,   305,   306,   307,   308,  1006,   310,   311,
     312,   313,   314,   315,   316,   317,   318,   319,   320,   321,
     322,   323,   324,   396,   286,   303,   304,   305,   306,   307,
     308,   396,   310,   311,   312,   313,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   396,  1008,   396,
      91,   396,  1012,   396,   254,  1117,   396,     3,     4,   396,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,   396,   396,    23,    24,    25,
      26,    27,    28,    29,    30,    31,   388,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    91,   396,   396,   417,
      46,    47,    48,    49,    50,    51,    52,    53,   416,   416,
     421,   396,  1072,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,  1085,  1086,   396,   396,   254,
     254,   254,   418,   396,   421,   396,   396,  1097,  1098,   396,
    1100,   396,  1102,   396,  1104,   421,   254,  1107,  1108,    95,
     254,    91,   254,   396,   396,   396,   102,   415,   396,   396,
     396,   415,   396,  1071,   396,   396,   396,   415,  1182,   415,
     396,   396,   396,   279,   396,   415,   396,  1085,  1086,   125,
     396,   396,   396,    91,   396,   396,   140,   415,   254,   254,
      91,  1099,    91,  1101,   140,  1103,    91,   416,    91,  1107,
    1108,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    91,   415,   397,    91,   415,
     397,   396,   396,  1183,   421,   422,  1186,    91,   396,   396,
      91,   177,    91,   179,   415,   396,   396,   183,   257,  1199,
    1200,   399,   400,   401,   402,   403,   404,   405,   406,   407,
     408,   409,   410,   411,  1214,   201,   396,   396,   204,   415,
      91,   396,  1276,   421,   396,  1225,  1226,  1227,   396,  1283,
    1230,  1231,  1232,   257,  1234,  1235,  1184,    26,   396,    28,
      29,   396,   396,   396,   396,   396,   396,  1247,  1248,    91,
     410,   237,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,   396,   396,   254,   396,
     396,   396,    61,   396,   415,   421,   396,   396,   396,   396,
     396,   396,   396,  1071,  1232,   396,   396,   415,   410,    91,
     283,   415,   396,    82,   396,   416,   396,  1085,  1086,   396,
     396,  1301,   396,   415,   283,   418,    95,   415,   415,   295,
     296,  1099,   415,  1101,   103,  1103,   415,   254,   396,  1107,
    1108,   133,   396,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   415,   411,   397,   396,   325,
     326,   327,   328,   329,   330,   331,   421,   422,   416,   415,
     396,  1351,   396,   396,   396,   396,  1356,   396,   396,   396,
     396,    91,   396,   140,   396,   396,   155,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,  1377,   411,  1379,
     396,  1381,   397,   396,  1384,   396,    91,   415,   421,    91,
     179,   418,   415,   396,  1394,   415,    91,   396,   187,   396,
     140,   190,    91,   392,   418,   396,   396,  1355,    91,   395,
     396,  1411,   201,   202,   203,   204,   205,   418,   207,   396,
     418,   407,   410,   410,   416,    91,   412,    91,   396,    91,
     415,   417,  1380,   415,   415,   397,  1436,  1437,   416,   415,
     396,  1441,    91,   254,  1232,   415,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
      91,   415,   414,   254,   396,   415,    91,  1467,   257,   421,
     415,    91,   397,   262,   396,   396,  1476,  1477,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,   415,   396,  1442,   415,   396,   286,   415,    91,
     415,   421,   396,   415,   418,  1505,  1506,   418,  1508,  1509,
     415,   396,    91,  1513,   396,   191,   396,   415,  1518,   415,
    1520,    91,   415,   415,   396,   415,   396,  1527,   396,   396,
     396,   415,   415,   396,   415,   396,   396,   415,  1538,   415,
    1540,   396,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,  1554,   411,   396,   191,  1506,   415,
    1560,   396,  1562,  1563,   415,   421,  1515,  1355,    80,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,  1583,    80,  1585,   415,   348,   934,   959,
     945,  1507,   421,    -1,  1594,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   396,   411,  1609,
      -1,   414,    -1,    -1,    -1,    -1,    -1,    -1,   421,    -1,
      -1,    -1,    -1,    -1,  1624,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1633,    -1,     0,     1,    -1,     3,     4,
    1640,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    -1,    -1,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    -1,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    -1,    -1,    -1,
      -1,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    -1,    -1,    58,    59,    -1,    -1,    -1,    63,    64,
      -1,    -1,    -1,    -1,    -1,    70,    71,    -1,    -1,    -1,
      75,    76,    77,    -1,    -1,    80,    -1,    -1,  1506,    84,
      -1,    -1,    87,    -1,    95,    96,    -1,    92,    -1,    94,
      95,    -1,    97,    98,    99,   100,    -1,   101,   103,    -1,
     105,   106,   107,    -1,   109,    -1,    -1,    -1,   119,    -1,
     115,    -1,    -1,    -1,   125,   120,    -1,   121,   123,    -1,
     131,    -1,   133,   128,   129,    -1,    -1,    -1,    -1,   134,
     135,   136,    -1,    -1,    -1,   140,   141,    -1,    -1,   144,
      -1,   152,    -1,    -1,    -1,   150,    -1,    -1,    -1,   154,
      -1,    -1,    -1,    -1,    -1,    -1,   161,   162,   163,   164,
      -1,   166,    -1,    -1,    -1,    -1,   177,   172,   179,    -1,
      -1,   176,   177,    -1,   179,   180,   181,   182,    -1,    -1,
     191,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   193,   194,
     201,    -1,    -1,   204,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   208,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   223,    -1,    -1,    -1,    -1,    -1,   222,    -1,
     225,   226,    -1,    -1,    -1,    -1,   231,   232,   233,   234,
     235,   236,   237,    -1,    -1,    -1,    -1,    -1,    -1,   244,
     245,   246,   247,   248,    -1,   250,   257,   258,    -1,   254,
     255,    -1,    -1,    -1,   259,    -1,   261,    -1,    -1,    -1,
     265,    -1,    -1,    -1,    -1,    -1,    -1,   278,    -1,    -1,
     275,   276,   277,   278,    -1,   286,    -1,   282,    -1,    -1,
      -1,   286,    -1,   288,    -1,    -1,   291,   292,   293,   294,
     295,   296,    -1,   298,   299,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,   322,   415,    -1,    -1,    -1,    -1,    -1,   421,    -1,
     325,   326,   327,   328,   329,   330,   331,   332,   333,    -1,
      -1,    -1,   337,   338,    -1,    -1,    -1,    -1,    -1,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   352,   353,
     354,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   380,   411,    -1,    -1,    -1,
     385,   386,   387,    -1,    -1,    -1,   421,    -1,    -1,    -1,
     395,   396,   396,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   407,    -1,    -1,    -1,    -1,   412,    -1,    -1,
       3,     4,   417,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    -1,    -1,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    -1,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    -1,
      -1,    -1,    -1,    46,    47,    48,    49,    50,    51,    52,
      53,    -1,    65,    -1,    67,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    90,   421,    -1,
      -1,    84,    95,    -1,    -1,    -1,    -1,    90,    -1,    -1,
      -1,    -1,    95,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,    -1,   131,    -1,
     133,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   131,    -1,
      -1,    -1,    -1,    -1,   147,    -1,    -1,   140,    -1,    -1,
      -1,    -1,    -1,    -1,   147,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   165,    -1,    -1,   168,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   168,   179,    -1,    -1,    -1,
     173,    -1,    -1,    -1,   177,    -1,   179,    -1,    -1,    -1,
     183,    -1,    -1,    -1,    -1,    -1,    -1,   101,   201,    -1,
      -1,   204,   205,    -1,    -1,    -1,    -1,    -1,   201,    -1,
      -1,   204,    -1,    -1,    -1,    -1,    -1,   121,    -1,    -1,
     223,   224,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,    -1,    -1,    -1,   415,
      -1,    -1,    -1,    -1,   237,   421,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   256,    -1,   258,    -1,    -1,    -1,    -1,
      -1,   254,    -1,    -1,    -1,    -1,   269,   270,    -1,    -1,
      -1,   274,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   275,   286,    -1,    -1,    -1,   290,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   290,    -1,    -1,
      -1,    -1,   295,   296,   297,    -1,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   222,   411,
      -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,
      -1,    -1,   325,   326,   327,   328,   329,   330,   331,     3,
       4,    -1,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    -1,    -1,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    -1,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    -1,    -1,
      -1,    -1,    46,    47,    48,    49,    50,    51,    52,    53,
      -1,    -1,    -1,   396,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   395,   396,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   407,    -1,    -1,    -1,    -1,   412,
      -1,    -1,    -1,    -1,   417,    -1,    90,    -1,    -1,    -1,
      -1,    95,    -1,    -1,    -1,    -1,    -1,    -1,   102,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   352,   353,
     354,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,    -1,    -1,   131,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   140,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   165,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   179,    -1,    -1,    -1,   183,
      -1,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,   201,   415,    -1,
     204,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,    -1,
      -1,    -1,    -1,   237,    -1,    -1,    -1,    -1,    -1,   243,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     254,     3,     4,    -1,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    -1,
      -1,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      -1,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      -1,   295,   296,    -1,    46,    47,    48,    49,    50,    51,
      52,    53,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,    -1,    -1,
      -1,   325,   326,   327,   328,   329,   330,   331,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,    -1,    95,    -1,    -1,    -1,    -1,   418,    -1,
     102,   421,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,
      -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   140,   411,
      -1,   395,   396,    -1,    -1,    -1,   418,    -1,    -1,   421,
      -1,    -1,    -1,   407,    -1,    -1,    -1,    -1,   412,    -1,
      -1,    -1,    -1,   417,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   177,    -1,   179,    -1,    -1,
      -1,   183,    -1,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,   201,
      -1,    -1,   204,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,    -1,    -1,    -1,    -1,   237,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   254,     3,     4,    -1,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    -1,    -1,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    -1,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    -1,   295,   296,    -1,    46,    47,    48,    49,
      50,    51,    52,    53,    -1,    -1,    -1,    -1,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,    -1,   325,   326,   327,   328,   329,   330,   331,
     420,   421,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    89,
      90,    -1,    -1,    -1,    -1,    95,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
      -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,
      -1,    -1,    -1,    -1,    -1,   125,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
     140,    -1,    -1,   395,   396,    -1,   418,    -1,    -1,   421,
      -1,    -1,    -1,    -1,    -1,   407,    -1,    -1,    -1,    -1,
     412,    -1,    -1,    -1,    -1,   417,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   177,    -1,   179,
       3,     4,    -1,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    -1,    -1,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    -1,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    -1,
      -1,    -1,    -1,    46,    47,    48,    49,    50,    51,    52,
      53,    -1,    -1,    -1,    -1,    -1,    -1,   237,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,    -1,    -1,   254,   415,    -1,   257,   258,    -1,
      -1,   421,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    95,    -1,   398,   399,   400,   401,   402,   403,
     404,   405,   406,   407,   408,   409,   410,   411,    -1,    -1,
      -1,    -1,    -1,    -1,   418,   295,   296,   421,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,    -1,    -1,    -1,    -1,    -1,   140,   418,    -1,
     143,   421,    -1,    -1,    -1,   325,   326,   327,   328,   329,
     330,   331,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,
      -1,    -1,   418,    -1,   177,   421,   179,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,   201,    -1,
     421,   204,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,   395,   396,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   421,    -1,    -1,   407,    -1,    -1,
      -1,    -1,   412,    -1,   237,    -1,    -1,   417,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   254,    -1,    -1,   257,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   295,   296,   398,   399,   400,   401,   402,   403,
     404,   405,   406,   407,   408,   409,   410,   411,    -1,    -1,
      -1,   415,    -1,    -1,    -1,    -1,    -1,   421,    -1,    -1,
      -1,    -1,   325,   326,   327,   328,   329,   330,   331,     3,
       4,    -1,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    -1,    -1,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    -1,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    -1,    -1,
      -1,    -1,    46,    47,    48,    49,    50,    51,    52,    53,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,   395,   396,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   421,    -1,    -1,   407,    -1,    -1,    -1,    -1,   412,
      -1,    -1,    -1,    -1,   417,    -1,    -1,    -1,    -1,    -1,
      -1,    95,    -1,    -1,    -1,    -1,    -1,     3,     4,    -1,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    -1,    -1,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    -1,    33,    34,    35,
      36,    37,    38,    39,    40,    41,   140,    -1,    -1,   143,
      46,    47,    48,    49,    50,    51,    52,    53,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,
      -1,   421,    -1,   177,    -1,   179,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
      -1,    -1,    -1,    -1,    -1,    -1,   418,   201,    -1,   421,
     204,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   237,   140,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     254,     3,     4,   257,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    -1,
      -1,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      -1,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      -1,   295,   296,    -1,    46,    47,    48,    49,    50,    51,
      52,    53,    -1,    -1,    -1,    -1,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
      -1,   325,   326,   327,   328,   329,   330,   331,   420,   421,
      -1,   237,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,    -1,    -1,   254,    -1,
      -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
      -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   140,   295,
     296,   395,   396,    -1,    -1,   301,   302,    -1,    -1,    -1,
      -1,    -1,    -1,   407,    -1,    -1,    -1,    -1,   412,    -1,
      -1,    -1,    -1,   417,    -1,    -1,    -1,    -1,    -1,   325,
     326,   327,   328,   329,   330,   331,    -1,    -1,    -1,   335,
     336,    -1,     3,     4,    -1,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      -1,    -1,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    -1,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    -1,    -1,    -1,    -1,    46,    47,    48,    49,    50,
      51,    52,    53,    -1,    -1,   237,    -1,    -1,    -1,   395,
     396,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   407,   254,    -1,    -1,    -1,   412,    -1,    -1,    -1,
      -1,   417,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,    -1,    -1,    -1,   415,
      -1,    -1,    -1,    -1,    -1,   421,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   295,   296,    -1,    -1,    -1,    -1,   301,
     302,    -1,    -1,    -1,    -1,    -1,    -1,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   140,
     411,    -1,    -1,   325,   326,   327,   328,   329,   330,   331,
     421,    -1,    -1,   335,   336,    -1,     3,     4,    -1,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    -1,    -1,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    -1,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    -1,    -1,    -1,    -1,    46,
      47,    48,    49,    50,    51,    52,    53,    -1,    -1,    -1,
      -1,    -1,    -1,   395,   396,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   407,    -1,    -1,    -1,    -1,
     412,    -1,    -1,    -1,    -1,   417,   237,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,   254,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   286,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   140,   295,   296,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
      -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,
      -1,    -1,    -1,    -1,   325,   326,   327,   328,   329,   330,
     331,    -1,     3,     4,    -1,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      -1,    -1,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    -1,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    -1,    -1,    -1,    -1,    46,    47,    48,    49,    50,
      51,    52,    53,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     237,    -1,    -1,    -1,   395,   396,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   407,   254,    -1,    -1,
      -1,   412,    -1,    -1,    -1,    -1,   417,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,   286,
     421,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   295,   296,
     398,   399,   400,   401,   402,   403,   404,   405,   406,   407,
     408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,   140,
     418,    -1,    -1,   421,    -1,    -1,    -1,    -1,   325,   326,
     327,   328,   329,   330,   331,    -1,     3,     4,    -1,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    -1,    -1,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    -1,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    -1,    -1,    -1,    -1,    46,
      47,    48,    49,    50,    51,    52,    53,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   395,   396,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     407,    -1,    -1,    -1,    -1,   412,   237,    -1,    -1,    -1,
     417,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   254,     3,     4,    -1,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    -1,    -1,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    -1,    33,    34,    35,    36,    37,    38,
      39,    40,    41,   140,   295,   296,    -1,    46,    47,    48,
      49,    50,    51,    52,    53,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   325,   326,   327,   328,   329,   330,
     331,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     237,   140,    -1,    -1,   395,   396,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   407,   254,    -1,    -1,
      -1,   412,    -1,    -1,    -1,    -1,   417,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   295,   296,
     398,   399,   400,   401,   402,   403,   404,   405,   406,   407,
     408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,
     418,    -1,    -1,   421,    -1,    -1,    -1,    -1,   325,   326,
     327,   328,   329,   330,   331,    -1,    -1,    -1,   237,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,   254,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   295,   296,   395,   396,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     407,    -1,    -1,    -1,    -1,   412,    -1,    -1,    -1,    -1,
     417,    -1,    -1,    -1,    -1,    -1,   325,   326,   327,   328,
     329,   330,   331,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,   395,   396,   415,    -1,
      -1,    -1,    -1,    -1,   421,    -1,    -1,    -1,   407,    -1,
      -1,    -1,    -1,   412,    -1,    -1,    -1,    -1,   417,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   420,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   420,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   420,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   420,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   420,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   420,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,    -1,
      -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,    -1,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,    -1,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,    -1,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,    -1,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,    -1,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,    -1,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,    -1,   411,    -1,    -1,    -1,
     415,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,   415,    -1,
      -1,    -1,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,    -1,   411,    -1,    -1,    -1,   415,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,    -1,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,    -1,   411,    -1,    -1,    -1,    -1,    -1,    -1,   418,
      -1,    -1,   421,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,    -1,   411,    -1,    -1,    -1,
      -1,    -1,    -1,   418,    -1,    -1,   421,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,    -1,
     411,    -1,    -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,
     421,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,    -1,   411,    -1,    -1,    -1,    -1,    -1,
      -1,   418,    -1,    -1,   421,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,    -1,   411,    -1,
      -1,    -1,    -1,    -1,    -1,   418,    -1,    -1,   421
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     1,     3,     4,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    58,    59,    63,
      64,    70,    71,    75,    76,    77,    80,    84,    87,    92,
      94,    95,    97,    98,    99,   100,   103,   105,   106,   107,
     109,   115,   120,   123,   128,   129,   134,   135,   136,   140,
     141,   144,   150,   154,   161,   162,   163,   164,   166,   172,
     176,   177,   179,   180,   181,   182,   193,   194,   208,   225,
     226,   231,   232,   233,   234,   235,   236,   237,   244,   245,
     246,   247,   248,   250,   254,   255,   259,   261,   265,   275,
     276,   277,   278,   282,   286,   288,   291,   292,   293,   294,
     295,   296,   298,   299,   325,   326,   327,   328,   329,   330,
     331,   332,   333,   337,   338,   380,   385,   386,   387,   395,
     396,   407,   412,   417,   424,   425,   427,   428,   429,   430,
     436,   438,   439,   440,   449,   453,   454,   466,   467,   468,
     469,   470,   471,   414,   417,   417,   417,   417,   417,   417,
     417,   417,   417,   417,   417,   417,   417,   417,   417,   417,
     417,   417,   417,   417,   417,   417,   417,   417,   417,   417,
     417,   396,   122,   285,   460,   460,   140,   254,   235,   254,
     298,   332,    90,    95,   131,   143,   177,   179,   201,   204,
     441,   451,   460,    95,    81,   200,   303,   304,   305,   306,
     307,   308,   310,   311,   312,   313,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   431,     3,     4,
       6,     7,     8,     9,    10,    11,    14,    15,    16,    17,
      18,    19,    20,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    95,   102,   125,   140,   177,   179,   183,   237,   254,
     395,   407,   412,   417,   451,   466,   471,    91,    84,   176,
     265,   275,   396,   396,   140,   254,   460,    91,   140,    90,
      95,   131,   177,   179,   278,   396,   417,   140,   471,    91,
      93,   130,   140,   253,   451,   460,    77,    95,   125,   177,
     179,   286,   451,    91,    28,    71,    78,   127,   148,   167,
     255,   266,   286,   416,   451,    28,   254,   255,   439,   440,
     105,   417,   305,    95,   396,   172,   270,   286,   396,   140,
     141,   254,   255,    84,    90,    95,   131,   147,   168,   173,
     177,   179,   183,   275,   290,   297,   327,   451,   471,    69,
      95,   102,   140,   177,   179,   183,   451,   471,   396,   396,
     208,   237,   395,   451,   140,   254,   471,    56,   155,   158,
     168,   171,   207,   241,   287,   396,   464,   124,   137,   138,
     144,   188,   189,   228,   230,   283,   426,   179,    91,    43,
      79,    80,    91,   154,   388,   431,   417,    95,   176,   177,
     179,   286,   305,   416,   451,    91,   110,   180,   220,   416,
     396,   108,   218,   455,   293,   396,    90,    95,   102,   131,
     140,   165,   179,   183,   243,   451,   471,    91,    95,   131,
     179,   257,   471,    91,    95,   131,   179,   257,   431,    69,
      95,   172,   177,   179,   183,   251,   301,   302,   335,   336,
     471,    84,   140,   176,   254,   275,   301,   302,   335,   336,
     471,   154,   254,   275,   441,   441,   417,   451,   451,   451,
     397,   419,    20,   140,   254,   395,   466,   470,   471,   470,
     471,   470,   471,     0,   414,   414,   414,   414,   414,   414,
     414,    78,   167,   279,   280,   442,   451,    95,    96,   119,
     125,   131,   133,   152,   177,   179,   191,   223,   257,   258,
     278,   286,   322,   437,   450,   451,   417,   417,   397,   419,
     414,   414,   414,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   411,   414,   421,   398,   399,
     400,   401,   402,   403,   404,   405,   406,   407,   408,   409,
     410,   411,   414,   421,   470,   471,   470,   471,   470,   471,
     470,   471,   470,   471,   470,   471,   470,   471,   470,   471,
     470,   471,   470,   471,   470,   471,   470,   471,   254,   470,
     471,   396,   470,   471,   470,   471,   470,   471,   470,   471,
     470,   471,   470,   471,   470,   471,   470,   471,   470,   471,
     470,   471,   470,   471,   470,   471,   470,   471,   416,   396,
     431,   451,   257,   396,   396,   451,   396,   396,   396,    91,
     417,   417,   417,   417,   417,   417,   417,   417,   417,   417,
     417,   417,   417,   417,   417,   417,   417,   417,   417,   417,
     417,   417,   417,   417,   417,   417,   417,   396,    95,   198,
     213,   456,   416,   396,   396,   291,   293,   465,   416,   416,
     419,   471,   471,   471,   419,   406,   407,   408,   409,   411,
     415,   415,   416,   283,   415,    91,   416,   257,   396,   261,
     396,   396,   396,   257,   254,   415,    95,   396,   451,   396,
     396,   396,   255,   286,   460,   257,   133,   223,   286,   305,
     451,   460,   451,   257,    78,    85,   145,   146,   184,   185,
     186,   221,   266,   267,   305,   309,   432,    46,    47,    48,
      49,    50,    51,    52,    53,   254,   255,   442,   396,   254,
     396,   396,   177,   179,   415,   147,   167,   168,   183,   257,
     396,   451,   471,   396,    95,   125,   177,   179,   451,   257,
     396,   396,   396,   396,   396,   396,   396,   465,   396,   396,
     471,   471,   415,   257,   286,   396,   396,   396,   396,   465,
     415,   451,   283,   396,   416,   283,   415,   396,   381,   382,
     383,   384,   433,   426,    91,   396,     5,    91,    91,   455,
      91,   140,   254,    91,   455,   254,   396,   471,   396,   396,
     171,   241,   339,   340,   341,   342,   434,   471,   466,   396,
     451,   455,   471,   173,   278,   466,   471,   257,   396,    91,
     396,   396,   396,   465,   396,   415,   396,   396,   396,   396,
     415,   396,   396,   396,   396,   286,   396,   167,   173,   305,
     451,   396,   396,   465,   396,   396,   396,   396,   396,   415,
     396,   396,   396,   471,   471,   471,   471,   415,    91,   124,
     133,   471,   254,   470,   471,   471,   417,   416,   416,   418,
     418,    95,   177,   179,   448,   451,    90,    91,    95,   131,
     165,   168,   179,   219,   447,    26,    28,    29,    61,    82,
      95,   103,   155,   179,   187,   190,   202,   203,   205,   207,
     257,   262,   286,   396,   443,   444,   451,   461,    65,    67,
      90,    95,   131,   133,   147,   165,   168,   179,   205,   223,
     224,   256,   258,   269,   270,   274,   286,   290,   396,   445,
     446,   451,    91,   173,   177,   179,   242,   286,    95,   213,
     292,   396,   452,   101,   121,   222,   343,   344,   345,   346,
     347,   348,   349,   350,   351,   352,   353,   354,   355,   356,
     357,   358,   359,   360,   361,   362,   363,   364,   365,   366,
     367,   368,   462,   471,   396,   396,    89,    90,    95,   125,
     177,   179,   257,   258,   471,   431,   471,   396,   254,   254,
     470,   471,   471,   470,   471,   470,   470,   470,   470,   470,
     470,   470,   470,   471,   470,   471,   470,   471,   470,   471,
     470,   471,   470,   471,   471,   471,   471,   471,   471,   471,
     471,   470,   471,   470,   471,   470,   471,   470,   471,   471,
     470,   471,   471,   418,   418,   418,   418,   418,   418,   418,
     418,   415,   415,   418,   418,   418,   418,   418,   418,   418,
     418,   418,   418,   418,   418,   415,   415,   418,   418,   418,
     418,   418,   418,   418,   418,   418,   418,   415,   415,   415,
     415,   415,   415,   415,   415,   418,   418,   415,   415,   418,
     418,   418,   418,   418,   418,   418,   418,   254,   396,   396,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   396,   471,   471,   471,   471,   471,   471,
     471,   471,   471,   471,   471,   471,   471,   396,   396,   254,
     466,   466,   471,   471,   471,   471,   471,   471,   471,   471,
     396,   254,   254,    91,   415,   254,   396,   396,   396,   415,
     471,   396,   396,    70,   262,   396,   462,   396,   396,   471,
     396,   110,   180,   220,   416,   396,   415,   415,   396,   396,
     471,   396,   133,   223,   257,   149,   289,   465,   396,   415,
     273,   396,   292,   451,   452,   396,   396,   396,    91,   471,
     396,   396,   471,   140,   415,   254,   254,   471,    91,    91,
      91,   416,    91,    91,   415,   415,   415,   397,   415,    91,
     415,   397,   397,   416,   419,   415,   396,   471,   471,   396,
      91,   396,   396,   471,   471,    91,    91,   415,   415,   420,
     254,   466,   396,   396,   396,   257,   396,   396,    89,   171,
     241,   458,   209,   215,   396,    70,   262,   396,   451,   471,
     471,   451,   396,   396,   396,    95,   143,   177,   179,   257,
     451,   471,    95,   143,   177,   179,   257,   451,   471,   471,
     471,    82,    83,   171,   241,   284,   457,   396,   396,    70,
     262,   415,   444,   396,    91,   257,   396,   396,   396,   462,
     396,   458,   149,   262,   289,   396,   457,   396,    91,    57,
     196,   199,   463,   396,   396,   286,   471,   286,   471,    70,
     262,   396,   415,   446,   396,   396,   396,   177,   179,   451,
     457,   452,   396,   460,   396,   471,   396,   396,   396,   471,
     396,   415,   415,   415,   420,   422,   422,   470,   471,   470,
     471,   470,   471,   471,   471,   470,   471,   470,   471,   470,
     471,   470,   471,   470,   471,   431,   418,   415,   418,   415,
     416,   419,   420,   420,   415,   415,   283,   415,    91,    44,
      91,   396,   462,   415,   415,   396,   451,   455,   471,   466,
     471,   396,   415,   462,   396,   396,   415,   471,   471,   452,
     396,   415,   471,   283,   254,   389,   390,   391,   435,   396,
     471,   471,   471,   396,   471,   471,   470,   471,    28,    29,
      73,   459,   471,   471,   415,   415,   415,   415,   133,   471,
     471,   397,   416,   396,   415,   396,   451,   396,   396,   396,
     396,   451,   396,   396,   396,   471,   396,    70,   262,    70,
     262,    91,   396,   396,   471,   396,   396,   397,   470,   471,
     418,   418,   418,   418,   418,   418,   415,   415,   418,   418,
     418,   418,   418,   418,   418,   418,   418,   418,   396,   471,
     471,   466,   471,   471,   396,   140,    91,   415,    91,   415,
     418,   396,   471,   415,    91,   415,   397,   416,   419,   415,
     396,   396,   471,   415,   140,    91,   392,   418,   415,   420,
     415,   396,   396,   471,   471,    91,   415,   415,   471,   466,
     396,   418,   418,   471,   471,   471,   418,   418,   419,   420,
     415,   415,   416,   415,    91,   415,    91,   396,   471,   470,
     471,   459,   471,   471,   415,   471,   416,   435,   415,   471,
     397,   471,   415,   415,    91,   471,   418,   418,   471,   471,
     396,   254,    91,   415,    91,   415,   420,   415,   471,   254,
     392,   393,   396,   415,   471,   415,   471,   471,   415,   420,
     415,    91,   415,    91,   397,   396,   396,   415,   471,   471,
     396,   396,   415,    91,   415,   471,   418,   415,   396,   279,
     418,   415,   396,   415,    91,   396,   191,   471,   396,   415,
      91,   415,   415,   396,   415,   396,   415,   396,   396,   471,
     396,   415,   415,   415,   415,   396,   396,   471,   396,   415,
     415,   415,   396,   471,   396,   191,   396
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
        break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
        break;
    }
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
/* Line 1792 of yacc.c  */
#line 532 "pars.y"
    {}
    break;

  case 4:
/* Line 1792 of yacc.c  */
#line 533 "pars.y"
    {}
    break;

  case 5:
/* Line 1792 of yacc.c  */
#line 534 "pars.y"
    {}
    break;

  case 6:
/* Line 1792 of yacc.c  */
#line 535 "pars.y"
    {
	    result = (yyvsp[(1) - (2)].val);
	}
    break;

  case 7:
/* Line 1792 of yacc.c  */
#line 538 "pars.y"
    {
	    result = *(yyvsp[(1) - (2)].ptr);
	}
    break;

  case 8:
/* Line 1792 of yacc.c  */
#line 541 "pars.y"
    {}
    break;

  case 9:
/* Line 1792 of yacc.c  */
#line 542 "pars.y"
    {}
    break;

  case 10:
/* Line 1792 of yacc.c  */
#line 543 "pars.y"
    {}
    break;

  case 11:
/* Line 1792 of yacc.c  */
#line 544 "pars.y"
    {}
    break;

  case 12:
/* Line 1792 of yacc.c  */
#line 545 "pars.y"
    {}
    break;

  case 13:
/* Line 1792 of yacc.c  */
#line 546 "pars.y"
    {}
    break;

  case 14:
/* Line 1792 of yacc.c  */
#line 547 "pars.y"
    {}
    break;

  case 15:
/* Line 1792 of yacc.c  */
#line 548 "pars.y"
    {
	    return 1;
	}
    break;

  case 16:
/* Line 1792 of yacc.c  */
#line 554 "pars.y"
    {
	    if ((yyvsp[(2) - (3)].pset) == FILEP) {
		set_printer(FILEP, (yyvsp[(3) - (3)].pset));
	    }
	    else {
		set_printer((yyvsp[(2) - (3)].pset), (yyvsp[(3) - (3)].pset));
	    }
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 17:
/* Line 1792 of yacc.c  */
#line 563 "pars.y"
    {
	    if ((yyvsp[(3) - (4)].pset) == FILEP) {
		set_printer(FILEP, (yyvsp[(4) - (4)].pset));
	    }
	    else {
		set_printer((yyvsp[(3) - (4)].pset), (yyvsp[(4) - (4)].pset));
	    }
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 18:
/* Line 1792 of yacc.c  */
#line 572 "pars.y"
    {
	    if ((yyvsp[(3) - (3)].pset) == FILEP) {
		set_printer(FILEP, NULL);
	    }
	    else {
		set_printer((yyvsp[(3) - (3)].pset), NULL);
	    }
	}
    break;

  case 19:
/* Line 1792 of yacc.c  */
#line 580 "pars.y"
    {
	    tdevice = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 20:
/* Line 1792 of yacc.c  */
#line 583 "pars.y"
    {
	    hdevice = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 21:
/* Line 1792 of yacc.c  */
#line 586 "pars.y"
    {
	    do_hardcopy();
	}
    break;

  case 22:
/* Line 1792 of yacc.c  */
#line 590 "pars.y"
    { 
	    extern int pslwfactor; /* fudge factor for linewidths in the PS driver */
	    pslwfactor = (yyvsp[(3) - (3)].val); 
	}
    break;

  case 23:
/* Line 1792 of yacc.c  */
#line 597 "pars.y"
    { (yyval.pset) = GR_PS_P; }
    break;

  case 24:
/* Line 1792 of yacc.c  */
#line 598 "pars.y"
    { (yyval.pset) = GR_PS_L; }
    break;

  case 25:
/* Line 1792 of yacc.c  */
#line 599 "pars.y"
    { (yyval.pset) = GR_MIF_P; }
    break;

  case 26:
/* Line 1792 of yacc.c  */
#line 600 "pars.y"
    { (yyval.pset) = GR_MIF_L; }
    break;

  case 27:
/* Line 1792 of yacc.c  */
#line 601 "pars.y"
    { (yyval.pset) = GR_GIFP; }
    break;

  case 28:
/* Line 1792 of yacc.c  */
#line 602 "pars.y"
    { (yyval.pset) = GR_GIFL; }
    break;

  case 29:
/* Line 1792 of yacc.c  */
#line 603 "pars.y"
    { (yyval.pset) = hdevice; }
    break;

  case 30:
/* Line 1792 of yacc.c  */
#line 604 "pars.y"
    { (yyval.pset) = FILEP; }
    break;

  case 31:
/* Line 1792 of yacc.c  */
#line 608 "pars.y"
    {
	    rg[(yyvsp[(1) - (2)].pset)].active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 32:
/* Line 1792 of yacc.c  */
#line 611 "pars.y"
    {
	    rg[(yyvsp[(1) - (3)].pset)].type = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 33:
/* Line 1792 of yacc.c  */
#line 614 "pars.y"
    {
	    rg[(yyvsp[(1) - (3)].pset)].color = checkon(COLOR, rg[(yyvsp[(1) - (3)].pset)].color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 34:
/* Line 1792 of yacc.c  */
#line 617 "pars.y"
    {
	    rg[(yyvsp[(1) - (3)].pset)].lines = checkon(LINESTYLE, rg[(yyvsp[(1) - (3)].pset)].lines, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 35:
/* Line 1792 of yacc.c  */
#line 620 "pars.y"
    {
	    rg[(yyvsp[(1) - (3)].pset)].linew = checkon(LINEWIDTH, rg[(yyvsp[(1) - (3)].pset)].linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 36:
/* Line 1792 of yacc.c  */
#line 624 "pars.y"
    {
	    rg[(yyvsp[(1) - (9)].pset)].x1 = (yyvsp[(3) - (9)].val);
	    rg[(yyvsp[(1) - (9)].pset)].y1 = (yyvsp[(5) - (9)].val);
	    rg[(yyvsp[(1) - (9)].pset)].x2 = (yyvsp[(7) - (9)].val);
	    rg[(yyvsp[(1) - (9)].pset)].y2 = (yyvsp[(9) - (9)].val);
	}
    break;

  case 37:
/* Line 1792 of yacc.c  */
#line 631 "pars.y"
    {
	    if (rg[(yyvsp[(1) - (5)].pset)].x == NULL || rg[(yyvsp[(1) - (5)].pset)].n == 0) {
		rg[(yyvsp[(1) - (5)].pset)].n = 0;
		rg[(yyvsp[(1) - (5)].pset)].x = (double *) calloc(1, sizeof(double));
		rg[(yyvsp[(1) - (5)].pset)].y = (double *) calloc(1, sizeof(double));
	    } else {
		rg[(yyvsp[(1) - (5)].pset)].x = (double *) realloc(rg[(yyvsp[(1) - (5)].pset)].x, (rg[(yyvsp[(1) - (5)].pset)].n + 1) * sizeof(double));
		rg[(yyvsp[(1) - (5)].pset)].y = (double *) realloc(rg[(yyvsp[(1) - (5)].pset)].y, (rg[(yyvsp[(1) - (5)].pset)].n + 1) * sizeof(double));
	    }
	    rg[(yyvsp[(1) - (5)].pset)].x[rg[(yyvsp[(1) - (5)].pset)].n] = (yyvsp[(3) - (5)].val);
	    rg[(yyvsp[(1) - (5)].pset)].y[rg[(yyvsp[(1) - (5)].pset)].n] = (yyvsp[(5) - (5)].val);
	    rg[(yyvsp[(1) - (5)].pset)].n++;
	}
    break;

  case 38:
/* Line 1792 of yacc.c  */
#line 644 "pars.y"
    {
	    rg[(yyvsp[(2) - (4)].pset)].linkto[(yyvsp[(4) - (4)].pset)] = TRUE;
	}
    break;

  case 39:
/* Line 1792 of yacc.c  */
#line 650 "pars.y"
    {
	    drawgraph();
	}
    break;

  case 40:
/* Line 1792 of yacc.c  */
#line 653 "pars.y"
    {
	    auto_redraw = ((yyvsp[(3) - (3)].pset) == ON);
	}
    break;

  case 41:
/* Line 1792 of yacc.c  */
#line 656 "pars.y"
    {
	    char buf[MAXPATHLEN];
	    strcpy(buf, (char *) (yyvsp[(2) - (2)].pset));
	    expand_tilde(buf); 
	    if (chdir(buf) >= 0) {
		strcpy(workingdir, buf);
	    	if (inwin) {
		    set_title(workingdir);
	        }
	    }
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 42:
/* Line 1792 of yacc.c  */
#line 668 "pars.y"
    {
	    if (inwin) {
		set_left_footer((yyvsp[(2) - (2)].pset));
	    }
	    else {
		printf("%s\n", (yyvsp[(2) - (2)].pset));
	    }
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 43:
/* Line 1792 of yacc.c  */
#line 677 "pars.y"
    {
	    setbgcolor((int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 44:
/* Line 1792 of yacc.c  */
#line 680 "pars.y"
    {
	    xlibsetcmap((int) (yyvsp[(2) - (8)].val), (int) (yyvsp[(4) - (8)].val), (int) (yyvsp[(6) - (8)].val), (int) (yyvsp[(8) - (8)].val));
	}
    break;

  case 45:
/* Line 1792 of yacc.c  */
#line 683 "pars.y"
    {
	    SetCorieTime((int) (yyvsp[(2) - (2)].pset) == TRUEP);
	}
    break;

  case 46:
/* Line 1792 of yacc.c  */
#line 686 "pars.y"
    {
	    exit(0);
	}
    break;

  case 47:
/* Line 1792 of yacc.c  */
#line 690 "pars.y"
    {
            set_pagelayout((yyvsp[(3) - (3)].pset));
        }
    break;

  case 48:
/* Line 1792 of yacc.c  */
#line 694 "pars.y"
    {
            set_toolbars((yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].pset) == ON);
        }
    break;

  case 49:
/* Line 1792 of yacc.c  */
#line 698 "pars.y"
    {
            set_toolbars((yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].pset) == ON);
        }
    break;

  case 50:
/* Line 1792 of yacc.c  */
#line 702 "pars.y"
    {
            set_toolbars((yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].pset) == ON);
        }
    break;

  case 51:
/* Line 1792 of yacc.c  */
#line 706 "pars.y"
    {
	    if (inwin) {
		my_draw2((double) (yyvsp[(2) - (4)].val), (double) (yyvsp[(4) - (4)].val));
		flush_pending();
	    }
	}
    break;

  case 52:
/* Line 1792 of yacc.c  */
#line 712 "pars.y"
    {
	    if (inwin) {
		my_move2((double) (yyvsp[(2) - (4)].val), (double) (yyvsp[(4) - (4)].val));
	    }
	}
    break;

  case 53:
/* Line 1792 of yacc.c  */
#line 717 "pars.y"
    {
	    if (inwin) {
		double x = (double) (yyvsp[(2) - (6)].val);
		double y = (double) (yyvsp[(4) - (6)].val);
		drawpolysym(&x, &y, 1, (int) (yyvsp[(6) - (6)].val), 0, 0, 1.0);
		flush_pending();
	    }
	}
    break;

  case 54:
/* Line 1792 of yacc.c  */
#line 725 "pars.y"
    {
	    if (inwin) {
		setcolor((int) (yyvsp[(2) - (2)].val));
	    }
	}
    break;

  case 55:
/* Line 1792 of yacc.c  */
#line 730 "pars.y"
    {
	    if (inwin) {
		setlinewidth((int) (yyvsp[(2) - (2)].val));
	    }
	}
    break;

  case 56:
/* Line 1792 of yacc.c  */
#line 735 "pars.y"
    {
	    if (inwin) {
		setlinestyle((int) (yyvsp[(2) - (2)].val));
	    }
	}
    break;

  case 57:
/* Line 1792 of yacc.c  */
#line 740 "pars.y"
    {
	}
    break;

  case 58:
/* Line 1792 of yacc.c  */
#line 744 "pars.y"
    {
	    switch ((yyvsp[(2) - (2)].pset)) {
	    case UP:
		gwindup_proc();
		break;
	    case DOWN:
		gwinddown_proc();
		break;
	    case RIGHT:
		gwindright_proc();
		break;
	    case LEFT:
		gwindleft_proc();
		break;
	    case IN:
		gwindshrink_proc();
		break;
	    case OUT:
		gwindexpand_proc();
		break;
	    }
	}
    break;

  case 59:
/* Line 1792 of yacc.c  */
#line 766 "pars.y"
    {
	    scroll_proc((int) (yyvsp[(2) - (2)].val));
	}
    break;

  case 60:
/* Line 1792 of yacc.c  */
#line 769 "pars.y"
    {
	    scrollinout_proc((int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 61:
/* Line 1792 of yacc.c  */
#line 772 "pars.y"
    {
	    scrolling_islinked = (yyvsp[(3) - (3)].pset) == ON;
	}
    break;

  case 62:
/* Line 1792 of yacc.c  */
#line 775 "pars.y"
    {
	    my_doublebuffer((yyvsp[(2) - (2)].pset) == TRUEP);
	}
    break;

  case 63:
/* Line 1792 of yacc.c  */
#line 778 "pars.y"
    {
	    my_frontbuffer((yyvsp[(2) - (2)].pset) == TRUEP);
	}
    break;

  case 64:
/* Line 1792 of yacc.c  */
#line 781 "pars.y"
    {
	    my_backbuffer((yyvsp[(2) - (2)].pset) == TRUEP);
	}
    break;

  case 65:
/* Line 1792 of yacc.c  */
#line 784 "pars.y"
    {
	    my_swapbuffer();
	}
    break;

  case 66:
/* Line 1792 of yacc.c  */
#line 787 "pars.y"
    {
	    sleep((int) (yyvsp[(2) - (2)].val));
	}
    break;

  case 67:
/* Line 1792 of yacc.c  */
#line 790 "pars.y"
    {	/* TODO add delay function */
	}
    break;

  case 68:
/* Line 1792 of yacc.c  */
#line 792 "pars.y"
    {		/* TODO add abort flag and function */
	}
    break;

  case 69:
/* Line 1792 of yacc.c  */
#line 795 "pars.y"
    {
	    gotparams = TRUE;
	    strcpy(paramfile, (char *) (yyvsp[(2) - (2)].pset));
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 70:
/* Line 1792 of yacc.c  */
#line 801 "pars.y"
    {
	    if (!fexists((char *) (yyvsp[(2) - (2)].pset))) {
		FILE *pp = fopen((char *) (yyvsp[(2) - (2)].pset), "w");
		if (pp != NULL) {
		    putparms(cg, pp, 0);
		    fclose(pp);
		} else {
		    errwin("Unable to write parameter file");
		}
	    }
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 71:
/* Line 1792 of yacc.c  */
#line 813 "pars.y"
    {
	    cg = (int) (yyvsp[(2) - (2)].pset);
	    g[cg].parmsread = TRUE;
	    change_gno = cg;
	}
    break;

  case 72:
/* Line 1792 of yacc.c  */
#line 818 "pars.y"
    {
	    curset = (int) (yyvsp[(2) - (2)].pset);
	}
    break;

  case 73:
/* Line 1792 of yacc.c  */
#line 823 "pars.y"
    {
	    set_hotlink(cg, (yyvsp[(1) - (4)].pset), 1, (char *) (yyvsp[(4) - (4)].pset), (yyvsp[(3) - (4)].pset));
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 74:
/* Line 1792 of yacc.c  */
#line 827 "pars.y"
    {
	    set_hotlink((yyvsp[(1) - (6)].pset), (yyvsp[(3) - (6)].pset), 1, (char *) (yyvsp[(6) - (6)].pset), (yyvsp[(5) - (6)].pset));
	    free((char *) (yyvsp[(6) - (6)].pset));
	}
    break;

  case 75:
/* Line 1792 of yacc.c  */
#line 831 "pars.y"
    {
	    set_hotlink(cg, (yyvsp[(1) - (3)].pset), (yyvsp[(3) - (3)].pset) == ON, NULL, 0);
	}
    break;

  case 76:
/* Line 1792 of yacc.c  */
#line 834 "pars.y"
    {
	    set_hotlink((yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].pset), (yyvsp[(5) - (5)].pset) == ON, NULL, 0);
	}
    break;

  case 77:
/* Line 1792 of yacc.c  */
#line 838 "pars.y"
    {
	    do_activateset(cg, (yyvsp[(2) - (3)].pset), (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 78:
/* Line 1792 of yacc.c  */
#line 841 "pars.y"
    {
            activateset(cg, (yyvsp[(2) - (4)].pset));
            settype(cg, (yyvsp[(2) - (4)].pset), (yyvsp[(3) - (4)].pset));
            setlength(cg, (yyvsp[(2) - (4)].pset), (yyvsp[(4) - (4)].val));
            setcomment(cg, (yyvsp[(2) - (4)].pset), "Generated set");
            updatesetminmax(cg, (yyvsp[(2) - (4)].pset));
	}
    break;

  case 79:
/* Line 1792 of yacc.c  */
#line 848 "pars.y"
    {
            activateset((yyvsp[(2) - (6)].pset), (yyvsp[(4) - (6)].pset));
            settype((yyvsp[(2) - (6)].pset), (yyvsp[(4) - (6)].pset), (yyvsp[(5) - (6)].pset));
            setlength((yyvsp[(2) - (6)].pset), (yyvsp[(4) - (6)].pset), (yyvsp[(6) - (6)].val));
            setcomment((yyvsp[(2) - (6)].pset), (yyvsp[(4) - (6)].pset), "Generated set");
            updatesetminmax((yyvsp[(2) - (6)].pset), (yyvsp[(4) - (6)].pset));
	}
    break;

  case 80:
/* Line 1792 of yacc.c  */
#line 855 "pars.y"
    {
	    add_point(cg, (yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].val), (yyvsp[(5) - (5)].val), 0.0, 0.0, XY);
	}
    break;

  case 81:
/* Line 1792 of yacc.c  */
#line 858 "pars.y"
    {
	    add_point((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), (yyvsp[(5) - (7)].val), (yyvsp[(7) - (7)].val), 0.0, 0.0, XY);
	}
    break;

  case 82:
/* Line 1792 of yacc.c  */
#line 862 "pars.y"
    {
	    int start = (int) (yyvsp[(3) - (5)].val) - 1;
	    int stop = (int) (yyvsp[(5) - (5)].val) - 1;
	    int dist = stop - start + 1;
	    if (dist > 0 && start >= 0) {
	        droppoints(cg, (yyvsp[(1) - (5)].pset), start, stop, dist);
	    }
	}
    break;

  case 83:
/* Line 1792 of yacc.c  */
#line 870 "pars.y"
    {
	    int start = (int) (yyvsp[(5) - (7)].val) - 1;
	    int stop = (int) (yyvsp[(7) - (7)].val) - 1;
	    int dist = stop - start + 1;
	    if (dist > 0 && start >= 0) {
	        droppoints((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), start, stop, dist);
	    }
	}
    break;

  case 84:
/* Line 1792 of yacc.c  */
#line 878 "pars.y"
    {
	    do_copyset(cg, (yyvsp[(2) - (4)].pset), cg, (yyvsp[(4) - (4)].pset));
	}
    break;

  case 85:
/* Line 1792 of yacc.c  */
#line 881 "pars.y"
    {
	    do_copyset((yyvsp[(2) - (8)].pset), (yyvsp[(4) - (8)].pset), (yyvsp[(6) - (8)].pset), (yyvsp[(8) - (8)].pset));
	}
    break;

  case 86:
/* Line 1792 of yacc.c  */
#line 884 "pars.y"
    {
	    do_moveset(cg, (yyvsp[(2) - (4)].pset), cg, (yyvsp[(4) - (4)].pset));
	}
    break;

  case 87:
/* Line 1792 of yacc.c  */
#line 887 "pars.y"
    {
	    do_moveset((yyvsp[(2) - (8)].pset), (yyvsp[(4) - (8)].pset), (yyvsp[(6) - (8)].pset), (yyvsp[(8) - (8)].pset));
	}
    break;

  case 88:
/* Line 1792 of yacc.c  */
#line 891 "pars.y"
    {
	    killset(cg, (yyvsp[(2) - (2)].pset));
	}
    break;

  case 89:
/* Line 1792 of yacc.c  */
#line 895 "pars.y"
    {
	    int i;
	    for (i = 0; i < g[cg].maxplot; i++) {
		killset(cg, i);
	    }
	}
    break;

  case 90:
/* Line 1792 of yacc.c  */
#line 902 "pars.y"
    {
	    kill_graph((yyvsp[(2) - (2)].pset));
	}
    break;

  case 91:
/* Line 1792 of yacc.c  */
#line 906 "pars.y"
    {
	    kill_graph(maxgraph);
	}
    break;

  case 92:
/* Line 1792 of yacc.c  */
#line 910 "pars.y"
    {
	    wipeout(0);
	}
    break;

  case 93:
/* Line 1792 of yacc.c  */
#line 914 "pars.y"
    {
	    int i;
	    for (i = 0; i < (int) (yyvsp[(3) - (7)].val); i++) {
		(yyvsp[(2) - (7)].ptr)[i] = (yyvsp[(5) - (7)].val) + (yyvsp[(7) - (7)].val) * i;
	    }
	}
    break;

  case 94:
/* Line 1792 of yacc.c  */
#line 921 "pars.y"
    {
	    int setno = (yyvsp[(3) - (6)].pset), ideg = (int) (yyvsp[(5) - (6)].val);
	    do_regress(setno, ideg, 0, -1, 0);
	}
    break;

  case 95:
/* Line 1792 of yacc.c  */
#line 926 "pars.y"
    {
	    do_running_command((yyvsp[(1) - (6)].pset), (yyvsp[(3) - (6)].pset), (int) (yyvsp[(5) - (6)].val));
	}
    break;

  case 96:
/* Line 1792 of yacc.c  */
#line 930 "pars.y"
    {
	    do_fourier_command((yyvsp[(1) - (6)].pset), (yyvsp[(3) - (6)].pset), (int) (yyvsp[(5) - (6)].val));
	}
    break;

  case 97:
/* Line 1792 of yacc.c  */
#line 934 "pars.y"
    {
	    do_spline((yyvsp[(3) - (10)].pset), (yyvsp[(5) - (10)].val), (yyvsp[(7) - (10)].val), (int) (yyvsp[(9) - (10)].val));
	}
    break;

  case 98:
/* Line 1792 of yacc.c  */
#line 938 "pars.y"
    {
	    do_histo_command((yyvsp[(3) - (10)].pset), -1, -1, (yyvsp[(5) - (10)].val), (yyvsp[(7) - (10)].val), (int) (yyvsp[(9) - (10)].val));
	}
    break;

  case 99:
/* Line 1792 of yacc.c  */
#line 942 "pars.y"
    {
	    do_differ((yyvsp[(3) - (6)].pset), (int) (yyvsp[(5) - (6)].val));
	}
    break;

  case 100:
/* Line 1792 of yacc.c  */
#line 946 "pars.y"
    {
	    do_int((yyvsp[(3) - (4)].pset), 0);
	}
    break;

  case 101:
/* Line 1792 of yacc.c  */
#line 950 "pars.y"
    {
	    if (activeset(cg)) {
		defaultgraph(cg);
		default_axis(cg, g[cg].auto_type, X_AXIS);
		default_axis(cg, g[cg].auto_type, ZX_AXIS);
		default_axis(cg, g[cg].auto_type, Y_AXIS);
		default_axis(cg, g[cg].auto_type, ZY_AXIS);
		update_world(cg);
		drawgraph();
	    } else {
		errwin("No active sets!");
	    }
	}
    break;

  case 102:
/* Line 1792 of yacc.c  */
#line 964 "pars.y"
    {
	    if (activeset(cg)) {
		defaultx(cg, -1);
		default_axis(cg, g[cg].auto_type, X_AXIS);
		default_axis(cg, g[cg].auto_type, ZX_AXIS);
		update_world(cg);
		drawgraph();
	    } else {
		errwin("No active sets!");
	    }
	}
    break;

  case 103:
/* Line 1792 of yacc.c  */
#line 976 "pars.y"
    {
	    if (activeset(cg)) {
		defaulty(cg, -1);
		default_axis(cg, g[cg].auto_type, Y_AXIS);
		default_axis(cg, g[cg].auto_type, ZY_AXIS);
		update_world(cg);
		drawgraph();
	    } else {
		errwin("No active sets!");
	    }
	}
    break;

  case 104:
/* Line 1792 of yacc.c  */
#line 988 "pars.y"
    {
	    if (isactive_set(cg, (yyvsp[(2) - (2)].pset))) {
		defaultsetgraph(cg, (yyvsp[(2) - (2)].pset));
		default_axis(cg, g[cg].auto_type, X_AXIS);
		default_axis(cg, g[cg].auto_type, ZX_AXIS);
		default_axis(cg, g[cg].auto_type, Y_AXIS);
		default_axis(cg, g[cg].auto_type, ZY_AXIS);
		update_world(cg);
		drawgraph();
	    } else {
		errwin("Set not active");
	    }
	}
    break;

  case 105:
/* Line 1792 of yacc.c  */
#line 1002 "pars.y"
    {
	    extern int go_locateflag;
	    go_locateflag = ((yyvsp[(2) - (2)].pset) == ON);
	}
    break;

  case 106:
/* Line 1792 of yacc.c  */
#line 1007 "pars.y"
    {
	    draw_focus(cg);
	    cg = (int) (yyvsp[(2) - (2)].pset);
	    defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, 
			islogx(cg), islogy(cg));
	    viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
	    draw_focus(cg);
	    update_all(cg);
	}
    break;

  case 107:
/* Line 1792 of yacc.c  */
#line 1016 "pars.y"
    {
	    draw_focus_flag = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 108:
/* Line 1792 of yacc.c  */
#line 1019 "pars.y"
    {
	    focus_policy = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 109:
/* Line 1792 of yacc.c  */
#line 1022 "pars.y"
    {
	    focus_policy = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 110:
/* Line 1792 of yacc.c  */
#line 1025 "pars.y"
    {
	    focus_policy = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 111:
/* Line 1792 of yacc.c  */
#line 1028 "pars.y"
    {
	    cursource = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 112:
/* Line 1792 of yacc.c  */
#line 1031 "pars.y"
    {
	    curtype = (yyvsp[(2) - (2)].pset);
	    change_type = curtype;
	}
    break;

  case 113:
/* Line 1792 of yacc.c  */
#line 1036 "pars.y"
    {
	    gotread = TRUE;
	    readtype = curtype;
	    readsrc = cursource;
	    strcpy(readfile, (char *) (yyvsp[(2) - (2)].pset));
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 114:
/* Line 1792 of yacc.c  */
#line 1044 "pars.y"
    {
	    gotbatch = TRUE;
	    strcpy(batchfile, (char *) (yyvsp[(3) - (3)].pset));
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 115:
/* Line 1792 of yacc.c  */
#line 1050 "pars.y"
    {
	    getdata(cg, (char *) (yyvsp[(3) - (3)].pset), DISK, BLOCK);
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 116:
/* Line 1792 of yacc.c  */
#line 1055 "pars.y"
    {
	    getdata(cg, (char *) (yyvsp[(4) - (4)].pset), (yyvsp[(3) - (4)].pset), BLOCK);
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 117:
/* Line 1792 of yacc.c  */
#line 1060 "pars.y"
    {
	    create_set_fromblock(cg, (yyvsp[(2) - (3)].pset), (char *) (yyvsp[(3) - (3)].pset));
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 118:
/* Line 1792 of yacc.c  */
#line 1065 "pars.y"
    {
	    gotread = TRUE;
	    readtype = (yyvsp[(2) - (3)].pset);
	    readsrc = cursource;
	    strcpy(readfile, (char *) (yyvsp[(3) - (3)].pset));
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 119:
/* Line 1792 of yacc.c  */
#line 1073 "pars.y"
    {
	    gotread = TRUE;
	    strcpy(readfile, (char *) (yyvsp[(4) - (4)].pset));
	    readtype = (yyvsp[(2) - (4)].pset);
	    readsrc = (yyvsp[(3) - (4)].pset);
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 120:
/* Line 1792 of yacc.c  */
#line 1081 "pars.y"
    {
	    read_image((char *) (yyvsp[(3) - (3)].pset));
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 121:
/* Line 1792 of yacc.c  */
#line 1086 "pars.y"
    {
	    getdata(cg, (char  *) (yyvsp[(4) - (4)].pset), DISK, CTD);
	}
    break;

  case 122:
/* Line 1792 of yacc.c  */
#line 1090 "pars.y"
    { /* set, file, type, TIME node, level, start, stop, skip, missing val*/
	}
    break;

  case 123:
/* Line 1792 of yacc.c  */
#line 1093 "pars.y"
    { /* set, file, type, TIME level, start, stop, skip*/
	}
    break;

  case 124:
/* Line 1792 of yacc.c  */
#line 1096 "pars.y"
    { /* set, file, type, PROFILE node*/
	}
    break;

  case 125:
/* Line 1792 of yacc.c  */
#line 1099 "pars.y"
    { /* set, file, type, TIME level, start, stop, skip*/
	}
    break;

  case 126:
/* Line 1792 of yacc.c  */
#line 1102 "pars.y"
    {
	    write_image((char *) (yyvsp[(3) - (3)].pset));
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 127:
/* Line 1792 of yacc.c  */
#line 1107 "pars.y"
    {
	    imagex = (int) (yyvsp[(3) - (5)].val);
	    imagey = (int) (yyvsp[(5) - (5)].val);
	}
    break;

  case 128:
/* Line 1792 of yacc.c  */
#line 1112 "pars.y"
    {
	    outputset(cg, (yyvsp[(2) - (2)].pset), (char *) NULL, (char *) NULL);
	}
    break;

  case 129:
/* Line 1792 of yacc.c  */
#line 1116 "pars.y"
    {
	    outputset(cg, (yyvsp[(2) - (4)].pset), (char *) NULL, (char *) (yyvsp[(4) - (4)].pset));
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 130:
/* Line 1792 of yacc.c  */
#line 1121 "pars.y"
    {
	    outputset(cg, (yyvsp[(2) - (4)].pset), (char *) (yyvsp[(4) - (4)].pset), (char *) NULL);
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 131:
/* Line 1792 of yacc.c  */
#line 1126 "pars.y"
    {
	    outputset(cg, (yyvsp[(2) - (6)].pset), (char *) (yyvsp[(4) - (6)].pset), (char *) (yyvsp[(6) - (6)].pset));
	    free((char *) (yyvsp[(4) - (6)].pset));
	    free((char *) (yyvsp[(6) - (6)].pset));
	}
    break;

  case 132:
/* Line 1792 of yacc.c  */
#line 1132 "pars.y"
    {
            extern char sformat[];
            do_writesets(maxgraph, -1, 1, (char *) (yyvsp[(2) - (2)].pset), sformat);
            free((char *) (yyvsp[(2) - (2)].pset));
        }
    break;

  case 133:
/* Line 1792 of yacc.c  */
#line 1137 "pars.y"
    {
	    push_world();
	}
    break;

  case 134:
/* Line 1792 of yacc.c  */
#line 1140 "pars.y"
    {
	    pop_world();
	}
    break;

  case 135:
/* Line 1792 of yacc.c  */
#line 1143 "pars.y"
    {
	    cycle_world_stack();
	}
    break;

  case 136:
/* Line 1792 of yacc.c  */
#line 1146 "pars.y"
    {
	    if ((int) (yyvsp[(2) - (2)].val) > 0)
		show_world_stack((int) (yyvsp[(2) - (2)].val) - 1);
	}
    break;

  case 137:
/* Line 1792 of yacc.c  */
#line 1151 "pars.y"
    {
	    add_world(cg, (yyvsp[(3) - (17)].val), (yyvsp[(5) - (17)].val), (yyvsp[(7) - (17)].val), (yyvsp[(9) - (17)].val), (yyvsp[(11) - (17)].val), (yyvsp[(13) - (17)].val), (yyvsp[(15) - (17)].val), (yyvsp[(17) - (17)].val));
	}
    break;

  case 138:
/* Line 1792 of yacc.c  */
#line 1154 "pars.y"
    {
	    clear_world_stack();
	}
    break;

  case 139:
/* Line 1792 of yacc.c  */
#line 1157 "pars.y"
    {
	    do_clear_boxes();
	}
    break;

  case 140:
/* Line 1792 of yacc.c  */
#line 1160 "pars.y"
    {
	    curbox = next_box();
	}
    break;

  case 141:
/* Line 1792 of yacc.c  */
#line 1163 "pars.y"
    {
	    curbox = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 142:
/* Line 1792 of yacc.c  */
#line 1166 "pars.y"
    {
	    boxes[curbox].active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 143:
/* Line 1792 of yacc.c  */
#line 1169 "pars.y"
    {
	    boxes[curbox].gno = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 144:
/* Line 1792 of yacc.c  */
#line 1173 "pars.y"
    {
	    if (curbox >= 0 && curbox < MAXBOXES) {
		boxes[curbox].x1 = (yyvsp[(2) - (8)].val);
		boxes[curbox].y1 = (yyvsp[(4) - (8)].val);
		boxes[curbox].x2 = (yyvsp[(6) - (8)].val);
		boxes[curbox].y2 = (yyvsp[(8) - (8)].val);
	    }
	}
    break;

  case 145:
/* Line 1792 of yacc.c  */
#line 1181 "pars.y"
    {
	    box_loctype = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 146:
/* Line 1792 of yacc.c  */
#line 1184 "pars.y"
    {
	    box_lines = checkon(LINESTYLE, box_lines, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 147:
/* Line 1792 of yacc.c  */
#line 1187 "pars.y"
    {
	    box_linew = checkon(LINEWIDTH, box_linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 148:
/* Line 1792 of yacc.c  */
#line 1190 "pars.y"
    {
	    box_color = checkon(COLOR, box_color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 149:
/* Line 1792 of yacc.c  */
#line 1193 "pars.y"
    {
	    box_fill = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 150:
/* Line 1792 of yacc.c  */
#line 1196 "pars.y"
    {
	    box_fillcolor = checkon(COLOR, box_fillcolor, (int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 151:
/* Line 1792 of yacc.c  */
#line 1199 "pars.y"
    {
	    box_fillpat = checkon(PATTERN, box_fillpat, (int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 152:
/* Line 1792 of yacc.c  */
#line 1203 "pars.y"
    {
	    if (curbox >= 0 && curbox < MAXBOXES) {
		boxes[curbox].lines = box_lines;
		boxes[curbox].linew = box_linew;
		boxes[curbox].color = box_color;
		boxes[curbox].fill = box_fill;
		boxes[curbox].fillcolor = box_fillcolor;
		boxes[curbox].fillpattern = box_fillpat;
		boxes[curbox].loctype = box_loctype;
	    }
	}
    break;

  case 153:
/* Line 1792 of yacc.c  */
#line 1214 "pars.y"
    {
	    curline = next_line();
	}
    break;

  case 154:
/* Line 1792 of yacc.c  */
#line 1217 "pars.y"
    {
	    curline = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 155:
/* Line 1792 of yacc.c  */
#line 1220 "pars.y"
    {
	    do_clear_lines();
	}
    break;

  case 156:
/* Line 1792 of yacc.c  */
#line 1223 "pars.y"
    {
	    lines[curline].active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 157:
/* Line 1792 of yacc.c  */
#line 1226 "pars.y"
    {
	    lines[curline].gno = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 158:
/* Line 1792 of yacc.c  */
#line 1230 "pars.y"
    {
	    lines[curline].x1 = (yyvsp[(2) - (8)].val);
	    lines[curline].y1 = (yyvsp[(4) - (8)].val);
	    lines[curline].x2 = (yyvsp[(6) - (8)].val);
	    lines[curline].y2 = (yyvsp[(8) - (8)].val);
	}
    break;

  case 159:
/* Line 1792 of yacc.c  */
#line 1236 "pars.y"
    {
	    line_loctype = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 160:
/* Line 1792 of yacc.c  */
#line 1239 "pars.y"
    {
	    line_linew = checkon(LINEWIDTH, line_linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 161:
/* Line 1792 of yacc.c  */
#line 1242 "pars.y"
    {
	    line_lines = checkon(LINESTYLE, line_lines, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 162:
/* Line 1792 of yacc.c  */
#line 1245 "pars.y"
    {
	    line_color = checkon(COLOR, line_color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 163:
/* Line 1792 of yacc.c  */
#line 1248 "pars.y"
    {
	    line_arrow = checkon(ARROW, line_arrow, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 164:
/* Line 1792 of yacc.c  */
#line 1251 "pars.y"
    {
	    line_asize = (yyvsp[(4) - (4)].val);
	}
    break;

  case 165:
/* Line 1792 of yacc.c  */
#line 1254 "pars.y"
    {
	    line_atype = (int) (yyvsp[(4) - (4)].val);
	}
    break;

  case 166:
/* Line 1792 of yacc.c  */
#line 1258 "pars.y"
    {
	    if (curline >= 0 && curline < MAXLINES) {
		lines[curline].lines = line_lines;
		lines[curline].linew = line_linew;
		lines[curline].color = line_color;
		lines[curline].arrow = line_arrow;
		lines[curline].asize = line_asize;
		lines[curline].atype = line_atype;
		lines[curline].loctype = line_loctype;
	    }
	}
    break;

  case 167:
/* Line 1792 of yacc.c  */
#line 1269 "pars.y"
    {
	    do_clear_text();
	}
    break;

  case 168:
/* Line 1792 of yacc.c  */
#line 1272 "pars.y"
    {
	    curstring = next_string();
	}
    break;

  case 169:
/* Line 1792 of yacc.c  */
#line 1275 "pars.y"
    {
	    curstring = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 170:
/* Line 1792 of yacc.c  */
#line 1278 "pars.y"
    {
	    pstr[curstring].active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 171:
/* Line 1792 of yacc.c  */
#line 1281 "pars.y"
    {
	    pstr[curstring].gno = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 172:
/* Line 1792 of yacc.c  */
#line 1285 "pars.y"
    {
	    pstr[curstring].x = (yyvsp[(2) - (4)].val);
	    pstr[curstring].y = (yyvsp[(4) - (4)].val);
	}
    break;

  case 173:
/* Line 1792 of yacc.c  */
#line 1289 "pars.y"
    {
	    string_loctype = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 174:
/* Line 1792 of yacc.c  */
#line 1292 "pars.y"
    {
	    string_linew = checkon(LINEWIDTH, string_linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 175:
/* Line 1792 of yacc.c  */
#line 1295 "pars.y"
    {
	    string_color = checkon(COLOR, string_color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 176:
/* Line 1792 of yacc.c  */
#line 1298 "pars.y"
    {
	    string_rot = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 177:
/* Line 1792 of yacc.c  */
#line 1301 "pars.y"
    {
	    string_font = checkon(FONTP, string_font, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 178:
/* Line 1792 of yacc.c  */
#line 1304 "pars.y"
    {
	    string_just = checkon(JUST, string_just, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 179:
/* Line 1792 of yacc.c  */
#line 1307 "pars.y"
    {
	    string_size = (yyvsp[(4) - (4)].val);
	}
    break;

  case 180:
/* Line 1792 of yacc.c  */
#line 1311 "pars.y"
    {
	    set_plotstr_string(&pstr[curstring], (char *) (yyvsp[(3) - (3)].pset));
	    pstr[curstring].linew = string_linew;
	    pstr[curstring].color = string_color;
	    pstr[curstring].font = string_font;
	    pstr[curstring].just = string_just;
	    pstr[curstring].loctype = string_loctype;
	    pstr[curstring].rot = string_rot;
	    pstr[curstring].charsize = string_size;
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 181:
/* Line 1792 of yacc.c  */
#line 1322 "pars.y"
    {
	    grdefaults.lines = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 182:
/* Line 1792 of yacc.c  */
#line 1325 "pars.y"
    {
	    grdefaults.linew = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 183:
/* Line 1792 of yacc.c  */
#line 1328 "pars.y"
    {
	    grdefaults.color = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 184:
/* Line 1792 of yacc.c  */
#line 1331 "pars.y"
    {
	    grdefaults.charsize = (yyvsp[(4) - (4)].val);
	}
    break;

  case 185:
/* Line 1792 of yacc.c  */
#line 1334 "pars.y"
    {
	    grdefaults.font = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 186:
/* Line 1792 of yacc.c  */
#line 1337 "pars.y"
    {
	    grdefaults.fontsrc = (int) (yyvsp[(4) - (4)].val);
	}
    break;

  case 187:
/* Line 1792 of yacc.c  */
#line 1340 "pars.y"
    {
	    grdefaults.symsize = (yyvsp[(4) - (4)].val);
	}
    break;

  case 188:
/* Line 1792 of yacc.c  */
#line 1344 "pars.y"
    {
	    g[cg].w.xg1 = (yyvsp[(2) - (8)].val);
	    g[cg].w.yg1 = (yyvsp[(4) - (8)].val);
	    g[cg].w.xg2 = (yyvsp[(6) - (8)].val);
	    g[cg].w.yg2 = (yyvsp[(8) - (8)].val);
	}
    break;

  case 189:
/* Line 1792 of yacc.c  */
#line 1350 "pars.y"
    {
	    g[cg].w.xg1 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 190:
/* Line 1792 of yacc.c  */
#line 1353 "pars.y"
    {
	    g[cg].w.xg2 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 191:
/* Line 1792 of yacc.c  */
#line 1356 "pars.y"
    {
	    g[cg].w.yg1 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 192:
/* Line 1792 of yacc.c  */
#line 1359 "pars.y"
    {
	    g[cg].w.yg2 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 193:
/* Line 1792 of yacc.c  */
#line 1363 "pars.y"
    {
	    g[cg].v.xv1 = (yyvsp[(2) - (8)].val);
	    g[cg].v.yv1 = (yyvsp[(4) - (8)].val);
	    g[cg].v.xv2 = (yyvsp[(6) - (8)].val);
	    g[cg].v.yv2 = (yyvsp[(8) - (8)].val);
	}
    break;

  case 194:
/* Line 1792 of yacc.c  */
#line 1369 "pars.y"
    {
	    g[cg].v.xv1 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 195:
/* Line 1792 of yacc.c  */
#line 1372 "pars.y"
    {
	    g[cg].v.xv2 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 196:
/* Line 1792 of yacc.c  */
#line 1375 "pars.y"
    {
	    g[cg].v.yv1 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 197:
/* Line 1792 of yacc.c  */
#line 1378 "pars.y"
    {
	    g[cg].v.yv2 = (yyvsp[(3) - (3)].val);
	}
    break;

  case 198:
/* Line 1792 of yacc.c  */
#line 1381 "pars.y"
    {
	    set_plotstr_string(&g[cg].labs.title, (char *) (yyvsp[(2) - (2)].pset));
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 199:
/* Line 1792 of yacc.c  */
#line 1385 "pars.y"
    {
	    g[cg].labs.title.font = checkon(FONTP, g[cg].labs.title.font, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 200:
/* Line 1792 of yacc.c  */
#line 1388 "pars.y"
    {
	    g[cg].labs.title.charsize = (yyvsp[(3) - (3)].val);
	}
    break;

  case 201:
/* Line 1792 of yacc.c  */
#line 1391 "pars.y"
    {
	    g[cg].labs.title.color = checkon(COLOR, g[cg].labs.title.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 202:
/* Line 1792 of yacc.c  */
#line 1395 "pars.y"
    {
	    g[cg].labs.title.linew = checkon(LINEWIDTH, g[cg].labs.title.linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 203:
/* Line 1792 of yacc.c  */
#line 1398 "pars.y"
    {
	    set_plotstr_string(&g[cg].labs.stitle, (char *) (yyvsp[(2) - (2)].pset));
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 204:
/* Line 1792 of yacc.c  */
#line 1403 "pars.y"
    {
	    g[cg].labs.stitle.font = checkon(FONTP, g[cg].labs.stitle.font, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 205:
/* Line 1792 of yacc.c  */
#line 1406 "pars.y"
    {
	    g[cg].labs.stitle.charsize = (yyvsp[(3) - (3)].val);
	}
    break;

  case 206:
/* Line 1792 of yacc.c  */
#line 1410 "pars.y"
    {
	    g[cg].labs.stitle.color = checkon(COLOR, g[cg].labs.stitle.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 207:
/* Line 1792 of yacc.c  */
#line 1414 "pars.y"
    {
	    g[cg].labs.stitle.linew = checkon(LINEWIDTH, g[cg].labs.stitle.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 208:
/* Line 1792 of yacc.c  */
#line 1417 "pars.y"
    {
	    realloc_plots((int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 209:
/* Line 1792 of yacc.c  */
#line 1420 "pars.y"
    {
	    g[cg].l.active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 210:
/* Line 1792 of yacc.c  */
#line 1423 "pars.y"
    {
	    g[cg].l.loctype = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 211:
/* Line 1792 of yacc.c  */
#line 1426 "pars.y"
    {
	    g[cg].l.layout = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 212:
/* Line 1792 of yacc.c  */
#line 1429 "pars.y"
    {
	    g[cg].l.vgap = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 213:
/* Line 1792 of yacc.c  */
#line 1432 "pars.y"
    {
	    g[cg].l.hgap = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 214:
/* Line 1792 of yacc.c  */
#line 1435 "pars.y"
    {
	    g[cg].l.len = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 215:
/* Line 1792 of yacc.c  */
#line 1438 "pars.y"
    {
	    g[cg].l.box = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 216:
/* Line 1792 of yacc.c  */
#line 1441 "pars.y"
    {
	    g[cg].l.boxfill = (yyvsp[(4) - (4)].pset);
	}
    break;

  case 217:
/* Line 1792 of yacc.c  */
#line 1444 "pars.y"
    {
	    g[cg].l.boxfillusing = (yyvsp[(5) - (5)].pset);
	}
    break;

  case 218:
/* Line 1792 of yacc.c  */
#line 1448 "pars.y"
    {
	    if ((yyvsp[(4) - (5)].pset) == COLOR) {
		g[cg].l.boxfillcolor = (int) (yyvsp[(5) - (5)].val);
	    } else {
		g[cg].l.boxfillpat = (int) (yyvsp[(5) - (5)].val);
	    }
	}
    break;

  case 219:
/* Line 1792 of yacc.c  */
#line 1455 "pars.y"
    {
	    g[cg].l.boxlcolor = checkon(COLOR, g[cg].l.boxlcolor, (int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 220:
/* Line 1792 of yacc.c  */
#line 1458 "pars.y"
    {
	    g[cg].l.boxlines = checkon(LINESTYLE, g[cg].l.boxlines, (int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 221:
/* Line 1792 of yacc.c  */
#line 1461 "pars.y"
    {
	    g[cg].l.boxlinew = checkon(LINEWIDTH, g[cg].l.boxlinew, (int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 222:
/* Line 1792 of yacc.c  */
#line 1464 "pars.y"
    {
	    g[cg].l.legx = (yyvsp[(2) - (4)].val);
	    g[cg].l.legy = (yyvsp[(4) - (4)].val);
	}
    break;

  case 223:
/* Line 1792 of yacc.c  */
#line 1468 "pars.y"
    {
	    g[cg].l.legx = (yyvsp[(3) - (3)].val);
	}
    break;

  case 224:
/* Line 1792 of yacc.c  */
#line 1471 "pars.y"
    {
	    g[cg].l.legy = (yyvsp[(3) - (3)].val);
	}
    break;

  case 225:
/* Line 1792 of yacc.c  */
#line 1474 "pars.y"
    {
	    g[cg].l.charsize = (yyvsp[(4) - (4)].val);
	}
    break;

  case 226:
/* Line 1792 of yacc.c  */
#line 1477 "pars.y"
    {
	    g[cg].l.font = checkon(FONTP, g[cg].l.font, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 227:
/* Line 1792 of yacc.c  */
#line 1480 "pars.y"
    {
	    g[cg].l.lines = checkon(LINESTYLE, g[cg].l.lines, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 228:
/* Line 1792 of yacc.c  */
#line 1483 "pars.y"
    {
	    g[cg].l.linew = checkon(LINEWIDTH, g[cg].l.linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 229:
/* Line 1792 of yacc.c  */
#line 1486 "pars.y"
    {
	    g[cg].l.color = checkon(COLOR, g[cg].l.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 230:
/* Line 1792 of yacc.c  */
#line 1489 "pars.y"
    {
	    set_plotstr_string(&g[cg].l.str[(int) (yyvsp[(3) - (4)].val)], (char *) (yyvsp[(4) - (4)].pset));
	    free((char *) (yyvsp[(4) - (4)].pset));
	}
    break;

  case 231:
/* Line 1792 of yacc.c  */
#line 1493 "pars.y"
    {
	    g[cg].isol.nisol = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 232:
/* Line 1792 of yacc.c  */
#line 1496 "pars.y"
    {
	    g[cg].isol.isoltype = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 233:
/* Line 1792 of yacc.c  */
#line 1499 "pars.y"
    {
	    g[cg].isol.cis[(int) (yyvsp[(2) - (4)].val)] = (yyvsp[(4) - (4)].val);
	}
    break;

  case 234:
/* Line 1792 of yacc.c  */
#line 1502 "pars.y"
    {
	}
    break;

  case 235:
/* Line 1792 of yacc.c  */
#line 1504 "pars.y"
    {
	}
    break;

  case 236:
/* Line 1792 of yacc.c  */
#line 1506 "pars.y"
    {
	}
    break;

  case 237:
/* Line 1792 of yacc.c  */
#line 1508 "pars.y"
    {
	    g[cg].isol.cis[0] = (double) (yyvsp[(3) - (5)].val);
	    g[cg].isol.cint = (double) (yyvsp[(5) - (5)].val);
	}
    break;

  case 238:
/* Line 1792 of yacc.c  */
#line 1512 "pars.y"
    {
	    g[cg].isol.lactive = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 239:
/* Line 1792 of yacc.c  */
#line 1515 "pars.y"
    {
	    g[cg].isol.loctype = (yyvsp[(4) - (4)].pset);
	}
    break;

  case 240:
/* Line 1792 of yacc.c  */
#line 1518 "pars.y"
    {
	    g[cg].isol.x = (double) (yyvsp[(3) - (5)].val);
	    g[cg].isol.y = (double) (yyvsp[(5) - (5)].val);
	}
    break;

  case 241:
/* Line 1792 of yacc.c  */
#line 1522 "pars.y"
    {
	    g[cg].isol.layout = (yyvsp[(4) - (4)].pset);
	}
    break;

  case 242:
/* Line 1792 of yacc.c  */
#line 1525 "pars.y"
    {
	    g[cg].isol.layout = (yyvsp[(4) - (4)].pset);
	}
    break;

  case 243:
/* Line 1792 of yacc.c  */
#line 1528 "pars.y"
    {
	    g[cg].isol.p.format = (yyvsp[(5) - (5)].pset);
	}
    break;

  case 244:
/* Line 1792 of yacc.c  */
#line 1531 "pars.y"
    {
	    g[cg].isol.p.prec = (int) (yyvsp[(5) - (5)].val);
	}
    break;

  case 245:
/* Line 1792 of yacc.c  */
#line 1534 "pars.y"
    {
	    g[cg].isol.p.charsize = (yyvsp[(5) - (5)].val);
	}
    break;

  case 246:
/* Line 1792 of yacc.c  */
#line 1537 "pars.y"
    {
	    g[cg].isol.xlen = (yyvsp[(4) - (6)].val);
	    g[cg].isol.ylen = (yyvsp[(6) - (6)].val);
	}
    break;

  case 247:
/* Line 1792 of yacc.c  */
#line 1541 "pars.y"
    {
	    g[cg].isol.xgap = (yyvsp[(4) - (6)].val);
	    g[cg].isol.ygap = (yyvsp[(6) - (6)].val);
	}
    break;

  case 248:
/* Line 1792 of yacc.c  */
#line 1545 "pars.y"
    {
	    g[cg].f.active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 249:
/* Line 1792 of yacc.c  */
#line 1548 "pars.y"
    {
	    g[cg].f.type = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 250:
/* Line 1792 of yacc.c  */
#line 1551 "pars.y"
    {
	    g[cg].f.lines = checkon(LINESTYLE, g[cg].f.lines, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 251:
/* Line 1792 of yacc.c  */
#line 1554 "pars.y"
    {
	    g[cg].f.linew = checkon(LINEWIDTH, g[cg].f.linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 252:
/* Line 1792 of yacc.c  */
#line 1557 "pars.y"
    {
	    g[cg].f.color = checkon(COLOR, g[cg].f.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 253:
/* Line 1792 of yacc.c  */
#line 1560 "pars.y"
    {
	    g[cg].f.fillbg = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 254:
/* Line 1792 of yacc.c  */
#line 1563 "pars.y"
    {
	    g[cg].f.bgcolor = (int) (yyvsp[(4) - (4)].val);
	}
    break;

  case 255:
/* Line 1792 of yacc.c  */
#line 1566 "pars.y"
    {
	    g[(yyvsp[(1) - (2)].pset)].active = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 256:
/* Line 1792 of yacc.c  */
#line 1569 "pars.y"
    {
	    g[(yyvsp[(1) - (3)].pset)].label = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 257:
/* Line 1792 of yacc.c  */
#line 1572 "pars.y"
    {
	    g[(yyvsp[(1) - (4)].pset)].auto_type = (yyvsp[(4) - (4)].pset);
	}
    break;

  case 258:
/* Line 1792 of yacc.c  */
#line 1575 "pars.y"
    {
	    g[(yyvsp[(1) - (4)].pset)].auto_type = (yyvsp[(4) - (4)].pset);
	}
    break;

  case 259:
/* Line 1792 of yacc.c  */
#line 1578 "pars.y"
    {
	    g[(yyvsp[(1) - (3)].pset)].parmsread = ((yyvsp[(3) - (3)].pset) == FALSEP);
	}
    break;

  case 260:
/* Line 1792 of yacc.c  */
#line 1581 "pars.y"
    {
	    g[(yyvsp[(1) - (3)].pset)].hidden = ((yyvsp[(3) - (3)].pset) == TRUEP);
	}
    break;

  case 261:
/* Line 1792 of yacc.c  */
#line 1584 "pars.y"
    {
	    g[(yyvsp[(1) - (3)].pset)].type = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 262:
/* Line 1792 of yacc.c  */
#line 1587 "pars.y"
    {
	    g[(yyvsp[(1) - (4)].pset)].barwid = (yyvsp[(4) - (4)].val);
	}
    break;

  case 263:
/* Line 1792 of yacc.c  */
#line 1590 "pars.y"
    {
	    g[(yyvsp[(1) - (4)].pset)].sbarwid = (yyvsp[(4) - (4)].val);
	}
    break;

  case 264:
/* Line 1792 of yacc.c  */
#line 1593 "pars.y"
    {
	    g[(yyvsp[(1) - (3)].pset)].pointset = ((yyvsp[(3) - (3)].pset) == ON);
	}
    break;

  case 265:
/* Line 1792 of yacc.c  */
#line 1597 "pars.y"
    {
	    g[(yyvsp[(1) - (5)].pset)].fx = (yyvsp[(4) - (5)].pset);
	    g[(yyvsp[(1) - (5)].pset)].fy = (yyvsp[(5) - (5)].pset);
	}
    break;

  case 266:
/* Line 1792 of yacc.c  */
#line 1602 "pars.y"
    {
	    g[(yyvsp[(1) - (6)].pset)].px = (yyvsp[(4) - (6)].val);
	    g[(yyvsp[(1) - (6)].pset)].py = (yyvsp[(6) - (6)].val);
	}
    break;

  case 267:
/* Line 1792 of yacc.c  */
#line 1607 "pars.y"
    {
	    g[(yyvsp[(1) - (6)].pset)].dsx = (yyvsp[(4) - (6)].val);
	    g[(yyvsp[(1) - (6)].pset)].dsy = (yyvsp[(6) - (6)].val);
	}
    break;

  case 268:
/* Line 1792 of yacc.c  */
#line 1611 "pars.y"
    {
	    g[(yyvsp[(1) - (4)].pset)].pt_type = (int) (yyvsp[(4) - (4)].val);
	}
    break;

  case 269:
/* Line 1792 of yacc.c  */
#line 1614 "pars.y"
    {
	    realloc_graph_plots((yyvsp[(1) - (4)].pset), (int) (yyvsp[(4) - (4)].val));
	}
    break;

  case 270:
/* Line 1792 of yacc.c  */
#line 1621 "pars.y"
    { /* set the DB host name, database */
#ifdef PGSQL
		SetDBHost((char *) (yyvsp[(2) - (4)].pset), (char *) (yyvsp[(4) - (4)].pset));
#endif
	}
    break;

  case 271:
/* Line 1792 of yacc.c  */
#line 1627 "pars.y"
    { /* graph set site, iid, x, y, cday1, cday2 */
#ifdef PGSQL
		ReadDB((yyvsp[(2) - (15)].pset), (yyvsp[(4) - (15)].pset), (char *) (yyvsp[(5) - (15)].pset), (char *) (yyvsp[(7) - (15)].pset), (char *) (yyvsp[(9) - (15)].pset), 
			(char *) (yyvsp[(11) - (15)].pset), (double) (yyvsp[(13) - (15)].val), (double) (yyvsp[(15) - (15)].val));
#endif
	}
    break;

  case 272:
/* Line 1792 of yacc.c  */
#line 1634 "pars.y"
    { /* default graph, next set site, iid, x, y, cday1, cday2 */
#ifdef PGSQL
		ReadDB(cg, -1, (char *) (yyvsp[(2) - (12)].pset), (char *) (yyvsp[(4) - (12)].pset), (char *) (yyvsp[(6) - (12)].pset), 
			(char *) (yyvsp[(8) - (12)].pset), (double) (yyvsp[(10) - (12)].val), (double) (yyvsp[(12) - (12)].val));
#endif
	}
    break;

  case 273:
/* Line 1792 of yacc.c  */
#line 1641 "pars.y"
    { /* default graph, next set site, iid, x, y, cday1, cday2 */
#ifdef PGSQL
		ReadDB((yyvsp[(2) - (13)].pset), -1, (char *) (yyvsp[(3) - (13)].pset), (char *) (yyvsp[(5) - (13)].pset), (char *) (yyvsp[(7) - (13)].pset), 
			(char *) (yyvsp[(9) - (13)].pset), (double) (yyvsp[(11) - (13)].val), (double) (yyvsp[(13) - (13)].val));
#endif
	}
    break;

  case 274:
/* Line 1792 of yacc.c  */
#line 1648 "pars.y"
    { /* graph set site, iid, x, y, bin#, cday1, cday2 */
#ifdef PGSQL
		ReadDBADP((yyvsp[(2) - (18)].pset), (yyvsp[(4) - (18)].pset), (char *) (yyvsp[(6) - (18)].pset), (char *) (yyvsp[(8) - (18)].pset), (char *) (yyvsp[(10) - (18)].pset), 
			(char *) (yyvsp[(12) - (18)].pset), (int) (yyvsp[(14) - (18)].val), (double) (yyvsp[(16) - (18)].val), (double) (yyvsp[(18) - (18)].val));
#endif
	}
    break;

  case 275:
/* Line 1792 of yacc.c  */
#line 1655 "pars.y"
    { /* graph set SQL */
#ifdef PGSQL
		ReadDBSQL((yyvsp[(2) - (5)].pset), (yyvsp[(4) - (5)].pset), (char *) (yyvsp[(5) - (5)].pset));
#endif
	}
    break;

  case 276:
/* Line 1792 of yacc.c  */
#line 1663 "pars.y"
    { g[cg].vp.vscale = (yyvsp[(3) - (3)].val); }
    break;

  case 277:
/* Line 1792 of yacc.c  */
#line 1664 "pars.y"
    { g[cg].vp.active = (yyvsp[(3) - (3)].pset); }
    break;

  case 278:
/* Line 1792 of yacc.c  */
#line 1665 "pars.y"
    { 
		    g[cg].vp.velx = (yyvsp[(4) - (6)].val); 
		    g[cg].vp.vely = (yyvsp[(6) - (6)].val); 
	}
    break;

  case 279:
/* Line 1792 of yacc.c  */
#line 1669 "pars.y"
    { g[cg].vp.userlength = (double) (yyvsp[(4) - (4)].val); }
    break;

  case 280:
/* Line 1792 of yacc.c  */
#line 1670 "pars.y"
    { set_plotstr_string(&g[cg].vp.vstr, (char *) (yyvsp[(4) - (4)].pset)); }
    break;

  case 281:
/* Line 1792 of yacc.c  */
#line 1671 "pars.y"
    { g[cg].vp.linew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 282:
/* Line 1792 of yacc.c  */
#line 1672 "pars.y"
    { g[cg].vp.lines = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 283:
/* Line 1792 of yacc.c  */
#line 1673 "pars.y"
    { g[cg].vp.color = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 284:
/* Line 1792 of yacc.c  */
#line 1674 "pars.y"
    { g[cg].vp.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 285:
/* Line 1792 of yacc.c  */
#line 1675 "pars.y"
    { g[cg].vp.arrowtype = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 286:
/* Line 1792 of yacc.c  */
#line 1679 "pars.y"
    {
	    (yyval.pset) = XY;
	}
    break;

  case 287:
/* Line 1792 of yacc.c  */
#line 1682 "pars.y"
    {
	    (yyval.pset) = XYARC;
	}
    break;

  case 288:
/* Line 1792 of yacc.c  */
#line 1685 "pars.y"
    {
	    (yyval.pset) = XYBOX;
	}
    break;

  case 289:
/* Line 1792 of yacc.c  */
#line 1688 "pars.y"
    {
	    (yyval.pset) = XYBOXPLOT;
	}
    break;

  case 290:
/* Line 1792 of yacc.c  */
#line 1691 "pars.y"
    {
	    (yyval.pset) = XYHILO;
	}
    break;

  case 291:
/* Line 1792 of yacc.c  */
#line 1694 "pars.y"
    {
	    (yyval.pset) = XYRT;
	}
    break;

  case 292:
/* Line 1792 of yacc.c  */
#line 1697 "pars.y"
    {
	    (yyval.pset) = XYSEG;
	}
    break;

  case 293:
/* Line 1792 of yacc.c  */
#line 1700 "pars.y"
    {
	    (yyval.pset) = XYSTRING;
	}
    break;

  case 294:
/* Line 1792 of yacc.c  */
#line 1703 "pars.y"
    {
	    (yyval.pset) = XYDX;
	}
    break;

  case 295:
/* Line 1792 of yacc.c  */
#line 1706 "pars.y"
    {
	    (yyval.pset) = XYDY;
	}
    break;

  case 296:
/* Line 1792 of yacc.c  */
#line 1709 "pars.y"
    {
	    (yyval.pset) = XYDXDX;
	}
    break;

  case 297:
/* Line 1792 of yacc.c  */
#line 1712 "pars.y"
    {
	    (yyval.pset) = XYDYDY;
	}
    break;

  case 298:
/* Line 1792 of yacc.c  */
#line 1715 "pars.y"
    {
	    (yyval.pset) = XYDXDY;
	}
    break;

  case 299:
/* Line 1792 of yacc.c  */
#line 1718 "pars.y"
    {
	    (yyval.pset) = XYX2Y2;
	}
    break;

  case 300:
/* Line 1792 of yacc.c  */
#line 1721 "pars.y"
    {
	    (yyval.pset) = XYXX;
	}
    break;

  case 301:
/* Line 1792 of yacc.c  */
#line 1724 "pars.y"
    {
	    (yyval.pset) = XYYY;
	}
    break;

  case 302:
/* Line 1792 of yacc.c  */
#line 1727 "pars.y"
    {
	    (yyval.pset) = XYZ;
	}
    break;

  case 303:
/* Line 1792 of yacc.c  */
#line 1730 "pars.y"
    {
	    (yyval.pset) = XYZW;
	}
    break;

  case 304:
/* Line 1792 of yacc.c  */
#line 1733 "pars.y"
    {
	    (yyval.pset) = XYUV;
	}
    break;

  case 305:
/* Line 1792 of yacc.c  */
#line 1736 "pars.y"
    {
	    (yyval.pset) = RECTGRID;
	}
    break;

  case 306:
/* Line 1792 of yacc.c  */
#line 1739 "pars.y"
    {
	    (yyval.pset) = FEGRID;
	}
    break;

  case 307:
/* Line 1792 of yacc.c  */
#line 1742 "pars.y"
    {
	    (yyval.pset) = NXY;
	}
    break;

  case 308:
/* Line 1792 of yacc.c  */
#line 1745 "pars.y"
    {
	    (yyval.pset) = BIN;
	}
    break;

  case 309:
/* Line 1792 of yacc.c  */
#line 1751 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 310:
/* Line 1792 of yacc.c  */
#line 1754 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 311:
/* Line 1792 of yacc.c  */
#line 1757 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 312:
/* Line 1792 of yacc.c  */
#line 1760 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 313:
/* Line 1792 of yacc.c  */
#line 1763 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 314:
/* Line 1792 of yacc.c  */
#line 1766 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 315:
/* Line 1792 of yacc.c  */
#line 1769 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 316:
/* Line 1792 of yacc.c  */
#line 1772 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 317:
/* Line 1792 of yacc.c  */
#line 1775 "pars.y"
    {
	    (yyval.pset) = XY;		/* not active */
	}
    break;

  case 318:
/* Line 1792 of yacc.c  */
#line 1778 "pars.y"
    {
	    (yyval.pset) = XY;		/* not active */
	}
    break;

  case 319:
/* Line 1792 of yacc.c  */
#line 1781 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 320:
/* Line 1792 of yacc.c  */
#line 1784 "pars.y"
    {
	    (yyval.pset) = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 321:
/* Line 1792 of yacc.c  */
#line 1790 "pars.y"
    {
            (yyval.pset) = FREE;
        }
    break;

  case 322:
/* Line 1792 of yacc.c  */
#line 1793 "pars.y"
    {
            (yyval.pset) = LANDSCAPE;
        }
    break;

  case 323:
/* Line 1792 of yacc.c  */
#line 1796 "pars.y"
    {
            (yyval.pset) = PORTRAIT;
        }
    break;

  case 324:
/* Line 1792 of yacc.c  */
#line 1799 "pars.y"
    {
            (yyval.pset) = FIXED;
        }
    break;

  case 325:
/* Line 1792 of yacc.c  */
#line 1805 "pars.y"
    {}
    break;

  case 326:
/* Line 1792 of yacc.c  */
#line 1806 "pars.y"
    {}
    break;

  case 327:
/* Line 1792 of yacc.c  */
#line 1807 "pars.y"
    {}
    break;

  case 328:
/* Line 1792 of yacc.c  */
#line 1808 "pars.y"
    {}
    break;

  case 329:
/* Line 1792 of yacc.c  */
#line 1809 "pars.y"
    {}
    break;

  case 330:
/* Line 1792 of yacc.c  */
#line 1810 "pars.y"
    {}
    break;

  case 331:
/* Line 1792 of yacc.c  */
#line 1814 "pars.y"
    {}
    break;

  case 332:
/* Line 1792 of yacc.c  */
#line 1815 "pars.y"
    {}
    break;

  case 333:
/* Line 1792 of yacc.c  */
#line 1816 "pars.y"
    {}
    break;

  case 334:
/* Line 1792 of yacc.c  */
#line 1820 "pars.y"
    {}
    break;

  case 335:
/* Line 1792 of yacc.c  */
#line 1835 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ACTIVE, (yyvsp[(1) - (1)].pset), 0);
	}
    break;

  case 336:
/* Line 1792 of yacc.c  */
#line 1838 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ACTIVE, (yyvsp[(1) - (1)].pset), 0);
	}
    break;

  case 337:
/* Line 1792 of yacc.c  */
#line 1841 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, TYPE, (yyvsp[(2) - (2)].pset), 0);
	}
    break;

  case 338:
/* Line 1792 of yacc.c  */
#line 1844 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, MISSINGP, (yyvsp[(2) - (2)].val), 0);
	}
    break;

  case 339:
/* Line 1792 of yacc.c  */
#line 1847 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, PREC, (int) (yyvsp[(2) - (2)].val), 0);
	}
    break;

  case 340:
/* Line 1792 of yacc.c  */
#line 1850 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, FORMAT, (yyvsp[(2) - (2)].pset), 0);
	}
    break;

  case 341:
/* Line 1792 of yacc.c  */
#line 1853 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, TYPE, (int) (yyvsp[(2) - (2)].val), 0);
	}
    break;

  case 342:
/* Line 1792 of yacc.c  */
#line 1856 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, FILL, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 343:
/* Line 1792 of yacc.c  */
#line 1859 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, CENTER, ((yyvsp[(3) - (3)].pset) == TRUEP), 0);
	}
    break;

  case 344:
/* Line 1792 of yacc.c  */
#line 1862 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, SIZE, (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 345:
/* Line 1792 of yacc.c  */
#line 1865 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, CHAR, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 346:
/* Line 1792 of yacc.c  */
#line 1868 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, SKIP, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 347:
/* Line 1792 of yacc.c  */
#line 1871 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, COLOR, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 348:
/* Line 1792 of yacc.c  */
#line 1874 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, LINEWIDTH, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 349:
/* Line 1792 of yacc.c  */
#line 1877 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SYMBOL, LINESTYLE, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 350:
/* Line 1792 of yacc.c  */
#line 1880 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, (yyvsp[(1) - (2)].pset), (int) (yyvsp[(2) - (2)].val), 0);
	}
    break;

  case 351:
/* Line 1792 of yacc.c  */
#line 1883 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, FILL, TYPE, (int) (yyvsp[(2) - (2)].val), 0);
	}
    break;

  case 352:
/* Line 1792 of yacc.c  */
#line 1886 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, FILL, WITH, (yyvsp[(3) - (3)].pset), 0);
	}
    break;

  case 353:
/* Line 1792 of yacc.c  */
#line 1889 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, FILL, (yyvsp[(2) - (3)].pset), (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 354:
/* Line 1792 of yacc.c  */
#line 1892 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, SKIP, (int) (yyvsp[(2) - (2)].val), 0);
	}
    break;

  case 355:
/* Line 1792 of yacc.c  */
#line 1895 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, LENGTH, (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 356:
/* Line 1792 of yacc.c  */
#line 1898 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, TYPE, (yyvsp[(3) - (3)].pset), 0);
	}
    break;

  case 357:
/* Line 1792 of yacc.c  */
#line 1901 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, LINEWIDTH, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 358:
/* Line 1792 of yacc.c  */
#line 1904 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, LINESTYLE, (int) (yyvsp[(3) - (3)].val), 0);
	}
    break;

  case 359:
/* Line 1792 of yacc.c  */
#line 1907 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, RISER, ACTIVE, (yyvsp[(3) - (3)].pset), 0);
	}
    break;

  case 360:
/* Line 1792 of yacc.c  */
#line 1910 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, RISER, LINEWIDTH, (int) (yyvsp[(4) - (4)].val), 0);
	}
    break;

  case 361:
/* Line 1792 of yacc.c  */
#line 1913 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, ERRORBAR, RISER, LINESTYLE, (int) (yyvsp[(4) - (4)].val), 0);
	}
    break;

  case 362:
/* Line 1792 of yacc.c  */
#line 1916 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, XYZ, (yyvsp[(2) - (4)].val), (yyvsp[(4) - (4)].val), 0);
	}
    break;

  case 363:
/* Line 1792 of yacc.c  */
#line 1919 "pars.y"
    {
	    set_prop(whichgraph, SET, SETNUM, whichset, COMMENT, (char *) (yyvsp[(2) - (2)].pset), 0);
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 364:
/* Line 1792 of yacc.c  */
#line 1926 "pars.y"
    {}
    break;

  case 365:
/* Line 1792 of yacc.c  */
#line 1927 "pars.y"
    {}
    break;

  case 366:
/* Line 1792 of yacc.c  */
#line 1928 "pars.y"
    {}
    break;

  case 367:
/* Line 1792 of yacc.c  */
#line 1929 "pars.y"
    {}
    break;

  case 368:
/* Line 1792 of yacc.c  */
#line 1930 "pars.y"
    {}
    break;

  case 369:
/* Line 1792 of yacc.c  */
#line 1931 "pars.y"
    {}
    break;

  case 370:
/* Line 1792 of yacc.c  */
#line 1935 "pars.y"
    {}
    break;

  case 371:
/* Line 1792 of yacc.c  */
#line 1936 "pars.y"
    {}
    break;

  case 372:
/* Line 1792 of yacc.c  */
#line 1937 "pars.y"
    {}
    break;

  case 373:
/* Line 1792 of yacc.c  */
#line 1938 "pars.y"
    {}
    break;

  case 374:
/* Line 1792 of yacc.c  */
#line 1939 "pars.y"
    {}
    break;

  case 375:
/* Line 1792 of yacc.c  */
#line 1940 "pars.y"
    {}
    break;

  case 376:
/* Line 1792 of yacc.c  */
#line 1944 "pars.y"
    {}
    break;

  case 377:
/* Line 1792 of yacc.c  */
#line 1945 "pars.y"
    {}
    break;

  case 378:
/* Line 1792 of yacc.c  */
#line 1946 "pars.y"
    {}
    break;

  case 379:
/* Line 1792 of yacc.c  */
#line 1950 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (1)].pset), 0.0);
	}
    break;

  case 380:
/* Line 1792 of yacc.c  */
#line 1953 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val));
	}
    break;

  case 381:
/* Line 1792 of yacc.c  */
#line 1956 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val));
	}
    break;

  case 382:
/* Line 1792 of yacc.c  */
#line 1959 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val));
	}
    break;

  case 383:
/* Line 1792 of yacc.c  */
#line 1962 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val));
	}
    break;

  case 384:
/* Line 1792 of yacc.c  */
#line 1965 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (3)].pset), (yyvsp[(3) - (3)].val));
	}
    break;

  case 385:
/* Line 1792 of yacc.c  */
#line 1968 "pars.y"
    {
	    set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].pset));
	}
    break;

  case 386:
/* Line 1792 of yacc.c  */
#line 1974 "pars.y"
    {}
    break;

  case 387:
/* Line 1792 of yacc.c  */
#line 1975 "pars.y"
    {}
    break;

  case 388:
/* Line 1792 of yacc.c  */
#line 1976 "pars.y"
    {}
    break;

  case 389:
/* Line 1792 of yacc.c  */
#line 1977 "pars.y"
    {}
    break;

  case 390:
/* Line 1792 of yacc.c  */
#line 1978 "pars.y"
    {
	    g[cg].t[naxis].active = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 391:
/* Line 1792 of yacc.c  */
#line 1984 "pars.y"
    {}
    break;

  case 392:
/* Line 1792 of yacc.c  */
#line 1985 "pars.y"
    {}
    break;

  case 393:
/* Line 1792 of yacc.c  */
#line 1990 "pars.y"
    {
	    g[cg].t[naxis].t_flag = (yyvsp[(1) - (1)].pset);
	    g[cg].t[naxis].t_mflag = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 394:
/* Line 1792 of yacc.c  */
#line 1994 "pars.y"
    {
	    g[cg].t[naxis].t_flag = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 395:
/* Line 1792 of yacc.c  */
#line 1997 "pars.y"
    {
	    g[cg].t[naxis].t_mflag = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 396:
/* Line 1792 of yacc.c  */
#line 2000 "pars.y"
    {
	    g[cg].t[naxis].tmajor = (yyvsp[(2) - (2)].val);
	}
    break;

  case 397:
/* Line 1792 of yacc.c  */
#line 2003 "pars.y"
    {
	    g[cg].t[naxis].tminor = (yyvsp[(2) - (2)].val);
	}
    break;

  case 398:
/* Line 1792 of yacc.c  */
#line 2006 "pars.y"
    {
	    g[cg].t[naxis].offsx = (yyvsp[(2) - (2)].val);
	}
    break;

  case 399:
/* Line 1792 of yacc.c  */
#line 2009 "pars.y"
    {
	    g[cg].t[naxis].offsy = (yyvsp[(2) - (2)].val);
	}
    break;

  case 400:
/* Line 1792 of yacc.c  */
#line 2012 "pars.y"
    {
	    g[cg].t[naxis].alt = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 401:
/* Line 1792 of yacc.c  */
#line 2015 "pars.y"
    {
	    g[cg].t[naxis].tmin = (yyvsp[(2) - (2)].val);
	}
    break;

  case 402:
/* Line 1792 of yacc.c  */
#line 2018 "pars.y"
    {
	    g[cg].t[naxis].tmax = (yyvsp[(2) - (2)].val);
	}
    break;

  case 403:
/* Line 1792 of yacc.c  */
#line 2021 "pars.y"
    {
	    g[cg].t[naxis].t_num = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 404:
/* Line 1792 of yacc.c  */
#line 2024 "pars.y"
    {
	    g[cg].t[naxis].t_inout = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 405:
/* Line 1792 of yacc.c  */
#line 2027 "pars.y"
    {
	    g[cg].t[naxis].t_log = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 406:
/* Line 1792 of yacc.c  */
#line 2030 "pars.y"
    {
	    g[cg].t[naxis].t_size = (yyvsp[(2) - (2)].val);
	}
    break;

  case 407:
/* Line 1792 of yacc.c  */
#line 2033 "pars.y"
    {
	    g[cg].t[naxis].t_size = (yyvsp[(3) - (3)].val);
	}
    break;

  case 408:
/* Line 1792 of yacc.c  */
#line 2036 "pars.y"
    {
	    g[cg].t[naxis].t_msize = (yyvsp[(3) - (3)].val);
	}
    break;

  case 409:
/* Line 1792 of yacc.c  */
#line 2039 "pars.y"
    {
	    g[cg].t[naxis].t_color = g[cg].t[naxis].t_mcolor = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 410:
/* Line 1792 of yacc.c  */
#line 2042 "pars.y"
    {
	    g[cg].t[naxis].t_linew = g[cg].t[naxis].t_mlinew = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 411:
/* Line 1792 of yacc.c  */
#line 2045 "pars.y"
    {
	    g[cg].t[naxis].t_color = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 412:
/* Line 1792 of yacc.c  */
#line 2048 "pars.y"
    {
	    g[cg].t[naxis].t_mcolor = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 413:
/* Line 1792 of yacc.c  */
#line 2051 "pars.y"
    {
	    g[cg].t[naxis].t_linew = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 414:
/* Line 1792 of yacc.c  */
#line 2054 "pars.y"
    {
	    g[cg].t[naxis].t_mlinew = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 415:
/* Line 1792 of yacc.c  */
#line 2057 "pars.y"
    {
	    g[cg].t[naxis].t_lines = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 416:
/* Line 1792 of yacc.c  */
#line 2060 "pars.y"
    {
	    g[cg].t[naxis].t_mlines = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 417:
/* Line 1792 of yacc.c  */
#line 2063 "pars.y"
    {
	    g[cg].t[naxis].t_gridflag = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 418:
/* Line 1792 of yacc.c  */
#line 2066 "pars.y"
    {
	    g[cg].t[naxis].t_mgridflag = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 419:
/* Line 1792 of yacc.c  */
#line 2069 "pars.y"
    {
	    g[cg].t[naxis].t_op = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 420:
/* Line 1792 of yacc.c  */
#line 2072 "pars.y"
    {
	    g[cg].t[naxis].t_type = AUTO;
	}
    break;

  case 421:
/* Line 1792 of yacc.c  */
#line 2075 "pars.y"
    {
	    g[cg].t[naxis].t_type = SPEC;
	}
    break;

  case 422:
/* Line 1792 of yacc.c  */
#line 2078 "pars.y"
    {
	    g[cg].t[naxis].t_spec = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 423:
/* Line 1792 of yacc.c  */
#line 2081 "pars.y"
    {
	    g[cg].t[naxis].t_specloc[(int) (yyvsp[(1) - (3)].val)] = (yyvsp[(3) - (3)].val);
	}
    break;

  case 426:
/* Line 1792 of yacc.c  */
#line 2092 "pars.y"
    {
	    g[cg].t[naxis].tl_flag = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 427:
/* Line 1792 of yacc.c  */
#line 2095 "pars.y"
    {
	    g[cg].t[naxis].tl_type = AUTO;
	}
    break;

  case 428:
/* Line 1792 of yacc.c  */
#line 2098 "pars.y"
    {
	    g[cg].t[naxis].tl_type = SPEC;
	}
    break;

  case 429:
/* Line 1792 of yacc.c  */
#line 2101 "pars.y"
    {
	    g[cg].t[naxis].tl_prec = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 430:
/* Line 1792 of yacc.c  */
#line 2104 "pars.y"
    {
	    g[cg].t[naxis].tl_format = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 431:
/* Line 1792 of yacc.c  */
#line 2107 "pars.y"
    {
	    g[cg].t[naxis].tl_format = (yyvsp[(2) - (2)].val);
	}
    break;

  case 432:
/* Line 1792 of yacc.c  */
#line 2110 "pars.y"
    {
	    strcpy(g[cg].t[naxis].tl_appstr, (char *) (yyvsp[(2) - (2)].pset));
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 433:
/* Line 1792 of yacc.c  */
#line 2114 "pars.y"
    {
	    strcpy(g[cg].t[naxis].tl_prestr, (char *) (yyvsp[(2) - (2)].pset));
	    free((char *) (yyvsp[(2) - (2)].pset));
	}
    break;

  case 434:
/* Line 1792 of yacc.c  */
#line 2118 "pars.y"
    {
	    g[cg].t[naxis].tl_layout = HORIZONTAL;
	}
    break;

  case 435:
/* Line 1792 of yacc.c  */
#line 2121 "pars.y"
    {
	    g[cg].t[naxis].tl_layout = VERTICAL;
	}
    break;

  case 436:
/* Line 1792 of yacc.c  */
#line 2124 "pars.y"
    {
	    g[cg].t[naxis].tl_layout = SPEC;
	}
    break;

  case 437:
/* Line 1792 of yacc.c  */
#line 2127 "pars.y"
    {
	    g[cg].t[naxis].tl_angle = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 438:
/* Line 1792 of yacc.c  */
#line 2130 "pars.y"
    {
	    g[cg].t[naxis].tl_just = (int) (yyvsp[(2) - (2)].pset);
	}
    break;

  case 439:
/* Line 1792 of yacc.c  */
#line 2133 "pars.y"
    {
	    g[cg].t[naxis].tl_skip = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 440:
/* Line 1792 of yacc.c  */
#line 2136 "pars.y"
    {
	    g[cg].t[naxis].tl_staggered = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 441:
/* Line 1792 of yacc.c  */
#line 2139 "pars.y"
    {
	    g[cg].t[naxis].tl_op = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 442:
/* Line 1792 of yacc.c  */
#line 2142 "pars.y"
    {
	    g[cg].t[naxis].tl_sign = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 443:
/* Line 1792 of yacc.c  */
#line 2145 "pars.y"
    {
	    g[cg].t[naxis].tl_start = (yyvsp[(2) - (2)].val);
	}
    break;

  case 444:
/* Line 1792 of yacc.c  */
#line 2148 "pars.y"
    {
	    g[cg].t[naxis].tl_stop = (yyvsp[(2) - (2)].val);
	}
    break;

  case 445:
/* Line 1792 of yacc.c  */
#line 2151 "pars.y"
    {
	    g[cg].t[naxis].tl_starttype = (int) (yyvsp[(3) - (3)].pset);
	}
    break;

  case 446:
/* Line 1792 of yacc.c  */
#line 2154 "pars.y"
    {
	    g[cg].t[naxis].tl_starttype = (int) (yyvsp[(3) - (3)].pset);
	}
    break;

  case 447:
/* Line 1792 of yacc.c  */
#line 2157 "pars.y"
    {
	    g[cg].t[naxis].tl_stoptype = (int) (yyvsp[(3) - (3)].pset);
	}
    break;

  case 448:
/* Line 1792 of yacc.c  */
#line 2160 "pars.y"
    {
	    g[cg].t[naxis].tl_stoptype = (int) (yyvsp[(3) - (3)].pset);
	}
    break;

  case 449:
/* Line 1792 of yacc.c  */
#line 2163 "pars.y"
    {
	    g[cg].t[naxis].tl_vgap = (yyvsp[(2) - (2)].val);
	}
    break;

  case 450:
/* Line 1792 of yacc.c  */
#line 2166 "pars.y"
    {
	    g[cg].t[naxis].tl_hgap = (yyvsp[(2) - (2)].val);
	}
    break;

  case 451:
/* Line 1792 of yacc.c  */
#line 2169 "pars.y"
    {
	    g[cg].t[naxis].tl_charsize = (yyvsp[(3) - (3)].val);
	}
    break;

  case 452:
/* Line 1792 of yacc.c  */
#line 2172 "pars.y"
    {
	    g[cg].t[naxis].tl_font = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 453:
/* Line 1792 of yacc.c  */
#line 2175 "pars.y"
    {
	    g[cg].t[naxis].tl_color = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 454:
/* Line 1792 of yacc.c  */
#line 2178 "pars.y"
    {
	    g[cg].t[naxis].tl_linew = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 455:
/* Line 1792 of yacc.c  */
#line 2181 "pars.y"
    {
	    set_plotstr_string(&g[cg].t[naxis].t_speclab[(int) (yyvsp[(1) - (3)].val)], (char *) (yyvsp[(3) - (3)].pset));
	    free((char *) (yyvsp[(3) - (3)].pset));
	}
    break;

  case 456:
/* Line 1792 of yacc.c  */
#line 2188 "pars.y"
    {
	    set_plotstr_string(&g[cg].t[naxis].label, (char *) (yyvsp[(1) - (1)].pset));
	    free((char *) (yyvsp[(1) - (1)].pset));
	}
    break;

  case 457:
/* Line 1792 of yacc.c  */
#line 2192 "pars.y"
    {
	    g[cg].t[naxis].label_layout = PERP;
	}
    break;

  case 458:
/* Line 1792 of yacc.c  */
#line 2195 "pars.y"
    {
	    g[cg].t[naxis].label_layout = PARA;
	}
    break;

  case 459:
/* Line 1792 of yacc.c  */
#line 2198 "pars.y"
    {
	    g[cg].t[naxis].label_place = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 460:
/* Line 1792 of yacc.c  */
#line 2201 "pars.y"
    {
	    g[cg].t[naxis].label_place = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 461:
/* Line 1792 of yacc.c  */
#line 2204 "pars.y"
    {
	    g[cg].t[naxis].label.x = (yyvsp[(2) - (4)].val);
	    g[cg].t[naxis].label.y = (yyvsp[(4) - (4)].val);
	}
    break;

  case 462:
/* Line 1792 of yacc.c  */
#line 2208 "pars.y"
    {
	    g[cg].t[naxis].label.just = (int) (yyvsp[(2) - (2)].pset);
	}
    break;

  case 463:
/* Line 1792 of yacc.c  */
#line 2211 "pars.y"
    {
	    g[cg].t[naxis].label.charsize = (yyvsp[(3) - (3)].val);
	}
    break;

  case 464:
/* Line 1792 of yacc.c  */
#line 2214 "pars.y"
    {
	    g[cg].t[naxis].label.font = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 465:
/* Line 1792 of yacc.c  */
#line 2217 "pars.y"
    {
	    g[cg].t[naxis].label.color = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 466:
/* Line 1792 of yacc.c  */
#line 2220 "pars.y"
    {
	    g[cg].t[naxis].label.linew = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 467:
/* Line 1792 of yacc.c  */
#line 2226 "pars.y"
    {
	    g[cg].t[naxis].t_drawbar = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 468:
/* Line 1792 of yacc.c  */
#line 2229 "pars.y"
    {
	    g[cg].t[naxis].t_drawbarcolor = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 469:
/* Line 1792 of yacc.c  */
#line 2232 "pars.y"
    {
	    g[cg].t[naxis].t_drawbarlines = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 470:
/* Line 1792 of yacc.c  */
#line 2235 "pars.y"
    {
	    g[cg].t[naxis].t_drawbarlinew = (int) (yyvsp[(2) - (2)].val);
	}
    break;

  case 471:
/* Line 1792 of yacc.c  */
#line 2242 "pars.y"
    {
	    whichgraph = (yyvsp[(1) - (3)].pset);
	    whichset = (yyvsp[(3) - (3)].pset);
	}
    break;

  case 472:
/* Line 1792 of yacc.c  */
#line 2247 "pars.y"
    {
	    whichgraph = cg;
	    whichset = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 473:
/* Line 1792 of yacc.c  */
#line 2252 "pars.y"
    {
	    whichgraph = cg;
	    whichset = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 474:
/* Line 1792 of yacc.c  */
#line 2257 "pars.y"
    {
	    whichgraph = (yyvsp[(1) - (2)].pset);
	    whichset = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 475:
/* Line 1792 of yacc.c  */
#line 2262 "pars.y"
    {
	    whichgraph = (yyvsp[(1) - (2)].pset);
	    whichset = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 476:
/* Line 1792 of yacc.c  */
#line 2267 "pars.y"
    {
	    whichgraph = (yyvsp[(1) - (2)].pset);
	    whichset = (yyvsp[(2) - (2)].pset);
	}
    break;

  case 477:
/* Line 1792 of yacc.c  */
#line 2274 "pars.y"
    {
            (yyval.pset) = (yyvsp[(1) - (1)].pset);
        }
    break;

  case 478:
/* Line 1792 of yacc.c  */
#line 2277 "pars.y"
    {
            (yyval.pset) = (yyvsp[(1) - (1)].pset);
        }
    break;

  case 479:
/* Line 1792 of yacc.c  */
#line 2280 "pars.y"
    {
            (yyval.pset) = (yyvsp[(1) - (1)].pset);
        }
    break;

  case 480:
/* Line 1792 of yacc.c  */
#line 2283 "pars.y"
    {
            (yyval.pset) = (yyvsp[(1) - (1)].pset);
        }
    break;

  case 481:
/* Line 1792 of yacc.c  */
#line 2286 "pars.y"
    {
            (yyval.pset) = (yyvsp[(1) - (1)].pset);
        }
    break;

  case 482:
/* Line 1792 of yacc.c  */
#line 2292 "pars.y"
    {
	    (yyval.pset) = ON;
	}
    break;

  case 483:
/* Line 1792 of yacc.c  */
#line 2295 "pars.y"
    {
	    (yyval.pset) = OFF;
	}
    break;

  case 484:
/* Line 1792 of yacc.c  */
#line 2301 "pars.y"
    {
	    (yyval.pset) = COLOR;
	}
    break;

  case 485:
/* Line 1792 of yacc.c  */
#line 2304 "pars.y"
    {
	    (yyval.pset) = PATTERN;
	}
    break;

  case 486:
/* Line 1792 of yacc.c  */
#line 2310 "pars.y"
    {
	    (yyval.pset) = RUNAVG;
	}
    break;

  case 487:
/* Line 1792 of yacc.c  */
#line 2313 "pars.y"
    {
	    (yyval.pset) = RUNSTD;
	}
    break;

  case 488:
/* Line 1792 of yacc.c  */
#line 2316 "pars.y"
    {
	    (yyval.pset) = RUNMED;
	}
    break;

  case 489:
/* Line 1792 of yacc.c  */
#line 2319 "pars.y"
    {
	    (yyval.pset) = RUNMAX;
	}
    break;

  case 490:
/* Line 1792 of yacc.c  */
#line 2322 "pars.y"
    {
	    (yyval.pset) = RUNMIN;
	}
    break;

  case 491:
/* Line 1792 of yacc.c  */
#line 2328 "pars.y"
    {
	    (yyval.pset) = DFT;
	}
    break;

  case 492:
/* Line 1792 of yacc.c  */
#line 2331 "pars.y"
    {
	    (yyval.pset) = FFT;
	}
    break;

  case 493:
/* Line 1792 of yacc.c  */
#line 2334 "pars.y"
    {
	    (yyval.pset) = INVDFT;
	}
    break;

  case 494:
/* Line 1792 of yacc.c  */
#line 2337 "pars.y"
    {
	    (yyval.pset) = INVFFT;
	}
    break;

  case 495:
/* Line 1792 of yacc.c  */
#line 2343 "pars.y"
    {
	    (yyval.pset) = DISK;
	}
    break;

  case 496:
/* Line 1792 of yacc.c  */
#line 2346 "pars.y"
    {
	    (yyval.pset) = PIPE;
	}
    break;

  case 497:
/* Line 1792 of yacc.c  */
#line 2352 "pars.y"
    {
	    (yyval.pset) = PATTERN;
	}
    break;

  case 498:
/* Line 1792 of yacc.c  */
#line 2355 "pars.y"
    {
	    (yyval.pset) = COLOR;
	}
    break;

  case 499:
/* Line 1792 of yacc.c  */
#line 2358 "pars.y"
    {
	    (yyval.pset) = NONE;
	}
    break;

  case 500:
/* Line 1792 of yacc.c  */
#line 2364 "pars.y"
    {
	    (yyval.pset) = TOP;
	}
    break;

  case 501:
/* Line 1792 of yacc.c  */
#line 2367 "pars.y"
    {
	    (yyval.pset) = BOTTOM;
	}
    break;

  case 502:
/* Line 1792 of yacc.c  */
#line 2370 "pars.y"
    {
	    (yyval.pset) = LEFT;
	}
    break;

  case 503:
/* Line 1792 of yacc.c  */
#line 2373 "pars.y"
    {
	    (yyval.pset) = RIGHT;
	}
    break;

  case 504:
/* Line 1792 of yacc.c  */
#line 2376 "pars.y"
    {
	    (yyval.pset) = BOTH;
	}
    break;

  case 505:
/* Line 1792 of yacc.c  */
#line 2382 "pars.y"
    {
	    (yyval.pset) = RIGHT;
	}
    break;

  case 506:
/* Line 1792 of yacc.c  */
#line 2385 "pars.y"
    {
	    (yyval.pset) = LEFT;
	}
    break;

  case 507:
/* Line 1792 of yacc.c  */
#line 2388 "pars.y"
    {
	    (yyval.pset) = CENTER;
	}
    break;

  case 508:
/* Line 1792 of yacc.c  */
#line 2394 "pars.y"
    {
	    (yyval.pset) = MINP;
	}
    break;

  case 509:
/* Line 1792 of yacc.c  */
#line 2397 "pars.y"
    {
	    (yyval.pset) = MAXP;
	}
    break;

  case 510:
/* Line 1792 of yacc.c  */
#line 2403 "pars.y"
    {
	    (yyval.pset) = TRUEP;
	}
    break;

  case 511:
/* Line 1792 of yacc.c  */
#line 2406 "pars.y"
    {
	    (yyval.pset) = FALSEP;
	}
    break;

  case 512:
/* Line 1792 of yacc.c  */
#line 2412 "pars.y"
    {
	    (yyval.pset) = IN;
	}
    break;

  case 513:
/* Line 1792 of yacc.c  */
#line 2415 "pars.y"
    {
	    (yyval.pset) = OUT;
	}
    break;

  case 514:
/* Line 1792 of yacc.c  */
#line 2418 "pars.y"
    {
	    (yyval.pset) = BOTH;
	}
    break;

  case 515:
/* Line 1792 of yacc.c  */
#line 2424 "pars.y"
    {
	    (yyval.pset) = DECIMAL;
	}
    break;

  case 516:
/* Line 1792 of yacc.c  */
#line 2427 "pars.y"
    {
	    (yyval.pset) = EXPONENTIAL;
	}
    break;

  case 517:
/* Line 1792 of yacc.c  */
#line 2430 "pars.y"
    {
	    (yyval.pset) = POWER;
	}
    break;

  case 518:
/* Line 1792 of yacc.c  */
#line 2433 "pars.y"
    {
	    (yyval.pset) = GENERAL;
	}
    break;

  case 519:
/* Line 1792 of yacc.c  */
#line 2436 "pars.y"
    {
	    (yyval.pset) = DDMMYY;
	}
    break;

  case 520:
/* Line 1792 of yacc.c  */
#line 2439 "pars.y"
    {
	    (yyval.pset) = MMDDYY;
	}
    break;

  case 521:
/* Line 1792 of yacc.c  */
#line 2442 "pars.y"
    {
	    (yyval.pset) = MMYY;
	}
    break;

  case 522:
/* Line 1792 of yacc.c  */
#line 2445 "pars.y"
    {
	    (yyval.pset) = MMDD;
	}
    break;

  case 523:
/* Line 1792 of yacc.c  */
#line 2448 "pars.y"
    {
	    (yyval.pset) = MONTHDAY;
	}
    break;

  case 524:
/* Line 1792 of yacc.c  */
#line 2451 "pars.y"
    {
	    (yyval.pset) = DAYMONTH;
	}
    break;

  case 525:
/* Line 1792 of yacc.c  */
#line 2454 "pars.y"
    {
	    (yyval.pset) = DDMONTHSYYHHMMSS;
	}
    break;

  case 526:
/* Line 1792 of yacc.c  */
#line 2457 "pars.y"
    {
	    (yyval.pset) = DDMONTHSYY;
	}
    break;

  case 527:
/* Line 1792 of yacc.c  */
#line 2460 "pars.y"
    {
	    (yyval.pset) = MONTHS;
	}
    break;

  case 528:
/* Line 1792 of yacc.c  */
#line 2463 "pars.y"
    {
	    (yyval.pset) = MONTHL;
	}
    break;

  case 529:
/* Line 1792 of yacc.c  */
#line 2466 "pars.y"
    {
	    (yyval.pset) = DAYOFWEEKS;
	}
    break;

  case 530:
/* Line 1792 of yacc.c  */
#line 2469 "pars.y"
    {
	    (yyval.pset) = DAYOFWEEKL;
	}
    break;

  case 531:
/* Line 1792 of yacc.c  */
#line 2472 "pars.y"
    {
	    (yyval.pset) = DAYOFYEAR;
	}
    break;

  case 532:
/* Line 1792 of yacc.c  */
#line 2475 "pars.y"
    {
	    (yyval.pset) = HMS;
	}
    break;

  case 533:
/* Line 1792 of yacc.c  */
#line 2478 "pars.y"
    {
	    (yyval.pset) = HH;
	}
    break;

  case 534:
/* Line 1792 of yacc.c  */
#line 2481 "pars.y"
    {
	    (yyval.pset) = MMDDHMS;
	}
    break;

  case 535:
/* Line 1792 of yacc.c  */
#line 2484 "pars.y"
    {
	    (yyval.pset) = MMDDYYHMS;
	}
    break;

  case 536:
/* Line 1792 of yacc.c  */
#line 2487 "pars.y"
    {
	    (yyval.pset) = DEGREESLON;
	}
    break;

  case 537:
/* Line 1792 of yacc.c  */
#line 2490 "pars.y"
    {
	    (yyval.pset) = DEGREESMMLON;
	}
    break;

  case 538:
/* Line 1792 of yacc.c  */
#line 2493 "pars.y"
    {
	    (yyval.pset) = DEGREESMMSSLON;
	}
    break;

  case 539:
/* Line 1792 of yacc.c  */
#line 2496 "pars.y"
    {
	    (yyval.pset) = MMSSLON;
	}
    break;

  case 540:
/* Line 1792 of yacc.c  */
#line 2499 "pars.y"
    {
	    (yyval.pset) = DEGREESLAT;
	}
    break;

  case 541:
/* Line 1792 of yacc.c  */
#line 2502 "pars.y"
    {
	    (yyval.pset) = DEGREESMMLAT;
	}
    break;

  case 542:
/* Line 1792 of yacc.c  */
#line 2505 "pars.y"
    {
	    (yyval.pset) = DEGREESMMSSLAT;
	}
    break;

  case 543:
/* Line 1792 of yacc.c  */
#line 2508 "pars.y"
    {
	    (yyval.pset) = MMSSLAT;
	}
    break;

  case 544:
/* Line 1792 of yacc.c  */
#line 2514 "pars.y"
    {
	    (yyval.pset) = NORMAL;
	}
    break;

  case 545:
/* Line 1792 of yacc.c  */
#line 2517 "pars.y"
    {
	    (yyval.pset) = ABSOLUTE;
	}
    break;

  case 546:
/* Line 1792 of yacc.c  */
#line 2520 "pars.y"
    {
	    (yyval.pset) = NEGATE;
	}
    break;

  case 547:
/* Line 1792 of yacc.c  */
#line 2526 "pars.y"
    {
	    (yyval.pset) = UP;
	}
    break;

  case 548:
/* Line 1792 of yacc.c  */
#line 2529 "pars.y"
    {
	    (yyval.pset) = DOWN;
	}
    break;

  case 549:
/* Line 1792 of yacc.c  */
#line 2532 "pars.y"
    {
	    (yyval.pset) = RIGHT;
	}
    break;

  case 550:
/* Line 1792 of yacc.c  */
#line 2535 "pars.y"
    {
	    (yyval.pset) = LEFT;
	}
    break;

  case 551:
/* Line 1792 of yacc.c  */
#line 2538 "pars.y"
    {
	    (yyval.pset) = IN;
	}
    break;

  case 552:
/* Line 1792 of yacc.c  */
#line 2541 "pars.y"
    {
	    (yyval.pset) = OUT;
	}
    break;

  case 553:
/* Line 1792 of yacc.c  */
#line 2547 "pars.y"
    {
	    (yyval.pset) = WORLD;
	}
    break;

  case 554:
/* Line 1792 of yacc.c  */
#line 2550 "pars.y"
    {
	    (yyval.pset) = VIEW;
	}
    break;

  case 555:
/* Line 1792 of yacc.c  */
#line 2556 "pars.y"
    { (yyval.pset) = X; }
    break;

  case 556:
/* Line 1792 of yacc.c  */
#line 2557 "pars.y"
    { (yyval.pset) = Y; }
    break;

  case 557:
/* Line 1792 of yacc.c  */
#line 2558 "pars.y"
    { (yyval.pset) = X0; }
    break;

  case 558:
/* Line 1792 of yacc.c  */
#line 2559 "pars.y"
    { (yyval.pset) = Y0; }
    break;

  case 559:
/* Line 1792 of yacc.c  */
#line 2560 "pars.y"
    { (yyval.pset) = Y1; }
    break;

  case 560:
/* Line 1792 of yacc.c  */
#line 2561 "pars.y"
    { (yyval.pset) = Y2; }
    break;

  case 561:
/* Line 1792 of yacc.c  */
#line 2562 "pars.y"
    { (yyval.pset) = Y3; }
    break;

  case 562:
/* Line 1792 of yacc.c  */
#line 2563 "pars.y"
    { (yyval.pset) = Y4; }
    break;

  case 563:
/* Line 1792 of yacc.c  */
#line 2564 "pars.y"
    { (yyval.pset) = Y5; }
    break;

  case 564:
/* Line 1792 of yacc.c  */
#line 2569 "pars.y"
    {
	    int itmp = (int) (yyvsp[(3) - (6)].val) - 1;
	    if (itmp >= ls) {
		yyerror("subscript out of range");
		return 1;
	    } else {
		(yyvsp[(1) - (6)].ptr)[itmp] = (yyvsp[(6) - (6)].val);
		result = (yyvsp[(6) - (6)].val);
	    }
	}
    break;

  case 565:
/* Line 1792 of yacc.c  */
#line 2580 "pars.y"
    {
	    int itmp = (int) (yyvsp[(3) - (6)].val) - 1;
	    double *ptr = getvptr(cg, curset, (yyvsp[(3) - (6)].val));
	    if (ptr != NULL) {
	        ptr[itmp] = (yyvsp[(6) - (6)].val);
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    result = (yyvsp[(6) - (6)].val);
	    updatesetminmax(cg, curset);
	    update_set_status(cg, curset);
	}
    break;

  case 566:
/* Line 1792 of yacc.c  */
#line 2595 "pars.y"
    {
	    int itmp = (int) (yyvsp[(5) - (8)].val) - 1;
	    double *ptr = getvptr(cg, (yyvsp[(1) - (8)].pset), (yyvsp[(3) - (8)].pset));
	    if (ptr != NULL) {
	        ptr[itmp] = (yyvsp[(8) - (8)].val);
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    result = (yyvsp[(8) - (8)].val);
	    updatesetminmax(cg, (yyvsp[(1) - (8)].pset));
	    update_set_status(cg, (yyvsp[(1) - (8)].pset));
	}
    break;

  case 567:
/* Line 1792 of yacc.c  */
#line 2610 "pars.y"
    {
	    set_prop(cg, SET, SETNUM, (yyvsp[(1) - (5)].pset), SYMBOL, TYPE, (int) (yyvsp[(5) - (5)].val), 0);
	    result = 0;
	}
    break;

  case 568:
/* Line 1792 of yacc.c  */
#line 2622 "pars.y"
    {
	    int itmp = (int) (yyvsp[(7) - (10)].val) - 1;
	    double *ptr = getvptr((yyvsp[(1) - (10)].pset), (yyvsp[(3) - (10)].pset), (yyvsp[(5) - (10)].pset));
	    if (ptr != NULL) {
	        ptr[itmp] = (yyvsp[(10) - (10)].val);
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    result = (yyvsp[(10) - (10)].val);
	    updatesetminmax((yyvsp[(1) - (10)].pset), (yyvsp[(3) - (10)].pset));
	    update_set_status((yyvsp[(1) - (10)].pset), (yyvsp[(3) - (10)].pset));
	}
    break;

  case 569:
/* Line 1792 of yacc.c  */
#line 2640 "pars.y"
    {
	    if ((yyvsp[(3) - (5)].pset) == X) {
		*xx = (yyvsp[(5) - (5)].val);
	    } else {
		*yy = (yyvsp[(5) - (5)].val);
	    }
	}
    break;

  case 570:
/* Line 1792 of yacc.c  */
#line 2651 "pars.y"
    {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[(1) - (3)].ptr)[i] = (yyvsp[(3) - (3)].ptr)[i];
	    }
	    result = (yyvsp[(3) - (3)].ptr)[0];
	}
    break;

  case 571:
/* Line 1792 of yacc.c  */
#line 2659 "pars.y"
    {
	    int i;
	    double *ptr;
	    if (!isactive_set(cg, curset)) {
		activateset(cg, curset);
		setlength(cg, curset, lxy);
		setcomment(cg, curset, "Created");
	    }
	    ptr = getvptr(cg, curset, (yyvsp[(1) - (3)].pset));
	    if (ptr != NULL) {
	        for (i = 0; i < lxy; i++) {
		    ptr[i] = (yyvsp[(3) - (3)].ptr)[i];
	        }
	        result = (yyvsp[(3) - (3)].ptr)[0];
	        updatesetminmax(cg, curset);
	        update_set_status(cg, curset);
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 572:
/* Line 1792 of yacc.c  */
#line 2681 "pars.y"
    {
	    int i;
	    double *ptr;
	    if (!isactive_set(cg, (yyvsp[(1) - (5)].pset))) {
		activateset(cg, (yyvsp[(1) - (5)].pset));
		setlength(cg, (yyvsp[(1) - (5)].pset), lxy);
		setcomment(cg, (yyvsp[(1) - (5)].pset), "Created");
	    }
	    ptr = getvptr(cg, (yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].pset));
	    if (ptr != NULL) {
	        for (i = 0; i < lxy; i++) {
		    ptr[i] = (yyvsp[(5) - (5)].ptr)[i];
	        }
	        result = (yyvsp[(5) - (5)].ptr)[0];
	        updatesetminmax(cg, (yyvsp[(1) - (5)].pset));
	        update_set_status(cg, (yyvsp[(1) - (5)].pset));
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 573:
/* Line 1792 of yacc.c  */
#line 2704 "pars.y"
    {
	    int i;
	    double *ptr;
	    if (!isactive_set((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset))) {
		activateset((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset));
		setlength((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), lxy);
		setcomment((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), "Created");
	    }
	    ptr = getvptr((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), (yyvsp[(5) - (7)].pset));
	    if (ptr != NULL) {
	        for (i = 0; i < lxy; i++) {
		    ptr[i] = (yyvsp[(7) - (7)].ptr)[i];
	        }
	        result = (yyvsp[(7) - (7)].ptr)[0];
	        updatesetminmax((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset));
	        update_set_status((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset));
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 574:
/* Line 1792 of yacc.c  */
#line 2727 "pars.y"
    {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[(1) - (3)].ptr)[i] = (yyvsp[(3) - (3)].val);
	    }
	    result = (yyvsp[(3) - (3)].val);
	}
    break;

  case 575:
/* Line 1792 of yacc.c  */
#line 2735 "pars.y"
    {
	    int i;
	    double *ptr;
	    if (!isactive_set(cg, curset)) {
		activateset(cg, curset);
		setlength(cg, curset, lxy);
		setcomment(cg, curset, "Created");
	    }
	    ptr = getvptr(cg, curset, (yyvsp[(1) - (3)].pset));
	    if (ptr != NULL) {
	        for (i = 0; i < lxy; i++) {
		    ptr[i] = (yyvsp[(3) - (3)].val);
	        }
	        result = (yyvsp[(3) - (3)].val);
	        updatesetminmax(cg, curset);
	        update_set_status(cg, curset);
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 576:
/* Line 1792 of yacc.c  */
#line 2758 "pars.y"
    {
	    int i;
	    double *ptr;
	    if (!isactive_set(cg, (yyvsp[(1) - (5)].pset))) {
		activateset(cg, (yyvsp[(1) - (5)].pset));
		setlength(cg, (yyvsp[(1) - (5)].pset), lxy);
		setcomment(cg, (yyvsp[(1) - (5)].pset), "Created");
	    }
	    ptr = getvptr(cg, (yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].pset));
	    if (ptr != NULL) {
	        for (i = 0; i < lxy; i++) {
		    ptr[i] = (yyvsp[(5) - (5)].val);
	        }
	        result = (yyvsp[(5) - (5)].val);
	        updatesetminmax(cg, (yyvsp[(1) - (5)].pset));
	        update_set_status(cg, (yyvsp[(1) - (5)].pset));
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 577:
/* Line 1792 of yacc.c  */
#line 2781 "pars.y"
    {
	    int i;
	    double *ptr;
	    if (!isactive_set((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset))) {
		activateset((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset));
		setlength((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), lxy);
		setcomment((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), "Created");
	    }
	    ptr = getvptr((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), (yyvsp[(5) - (7)].pset));
	    if (ptr != NULL) {
	        for (i = 0; i < lxy; i++) {
		    ptr[i] = (yyvsp[(7) - (7)].val);
	        }
	        result = (yyvsp[(7) - (7)].val);
	        updatesetminmax((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset));
	        update_set_status((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset));
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 578:
/* Line 1792 of yacc.c  */
#line 2807 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (1)].ptr)[i];
	    }
	}
    break;

  case 579:
/* Line 1792 of yacc.c  */
#line 2816 "pars.y"
    {
	    int i;
	    double *ptr = getvptr(cg, curset, (yyvsp[(1) - (1)].pset));
	    if (ptr == NULL) {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ptr[i];
	    }
	}
    break;

  case 580:
/* Line 1792 of yacc.c  */
#line 2830 "pars.y"
    {
	    int i;
	    double *ptr = getvptr(cg, (yyvsp[(1) - (3)].pset), (yyvsp[(3) - (3)].pset));
	    if (ptr == NULL) {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ptr[i];
	    }
	}
    break;

  case 581:
/* Line 1792 of yacc.c  */
#line 2844 "pars.y"
    {
	    int i;
	    double *ptr = getvptr((yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].pset), (yyvsp[(5) - (5)].pset));
	    if (ptr == NULL) {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ptr[i];
	    }
	}
    break;

  case 582:
/* Line 1792 of yacc.c  */
#line 2858 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (1)].val);
	    }
	}
    break;

  case 583:
/* Line 1792 of yacc.c  */
#line 2867 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 584:
/* Line 1792 of yacc.c  */
#line 2876 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] + (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 585:
/* Line 1792 of yacc.c  */
#line 2885 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 586:
/* Line 1792 of yacc.c  */
#line 2894 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] + (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 587:
/* Line 1792 of yacc.c  */
#line 2903 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 588:
/* Line 1792 of yacc.c  */
#line 2912 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] - (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 589:
/* Line 1792 of yacc.c  */
#line 2921 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 590:
/* Line 1792 of yacc.c  */
#line 2930 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] - (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 591:
/* Line 1792 of yacc.c  */
#line 2939 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 592:
/* Line 1792 of yacc.c  */
#line 2948 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] * (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 593:
/* Line 1792 of yacc.c  */
#line 2957 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 594:
/* Line 1792 of yacc.c  */
#line 2966 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] * (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 595:
/* Line 1792 of yacc.c  */
#line 2975 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[(3) - (3)].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 596:
/* Line 1792 of yacc.c  */
#line 2988 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[(3) - (3)].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] / (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 597:
/* Line 1792 of yacc.c  */
#line 3003 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[(3) - (3)].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 598:
/* Line 1792 of yacc.c  */
#line 3018 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[(3) - (3)].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] / (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 599:
/* Line 1792 of yacc.c  */
#line 3031 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val));
	    }
	}
    break;

  case 600:
/* Line 1792 of yacc.c  */
#line 3040 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].ptr)[i]);
	    }
	}
    break;

  case 601:
/* Line 1792 of yacc.c  */
#line 3049 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].ptr)[i], (yyvsp[(3) - (3)].val));
	    }
	}
    break;

  case 602:
/* Line 1792 of yacc.c  */
#line 3058 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].ptr)[i], (yyvsp[(3) - (3)].ptr)[i]);
	    }
	}
    break;

  case 603:
/* Line 1792 of yacc.c  */
#line 3067 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[(3) - (4)].val));
	    }
	}
    break;

  case 604:
/* Line 1792 of yacc.c  */
#line 3076 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 605:
/* Line 1792 of yacc.c  */
#line 3085 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = acos((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 606:
/* Line 1792 of yacc.c  */
#line 3094 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = asin((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 607:
/* Line 1792 of yacc.c  */
#line 3103 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 608:
/* Line 1792 of yacc.c  */
#line 3112 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan2((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 609:
/* Line 1792 of yacc.c  */
#line 3121 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ceil((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 610:
/* Line 1792 of yacc.c  */
#line 3130 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = cos((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 611:
/* Line 1792 of yacc.c  */
#line 3139 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] *= M_PI / 180.0;
	    }
	}
    break;

  case 612:
/* Line 1792 of yacc.c  */
#line 3148 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erf((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 613:
/* Line 1792 of yacc.c  */
#line 3157 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erfc((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 614:
/* Line 1792 of yacc.c  */
#line 3166 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = exp((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 615:
/* Line 1792 of yacc.c  */
#line 3175 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = floor((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 616:
/* Line 1792 of yacc.c  */
#line 3184 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = my_hypot((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 617:
/* Line 1792 of yacc.c  */
#line 3193 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = my_hypot((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 618:
/* Line 1792 of yacc.c  */
#line 3202 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = my_hypot((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].val));
	    }
	}
    break;

  case 619:
/* Line 1792 of yacc.c  */
#line 3211 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = my_hypot((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	    }
	}
    break;

  case 620:
/* Line 1792 of yacc.c  */
#line 3220 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = i + 1;
	    }
	}
    break;

  case 621:
/* Line 1792 of yacc.c  */
#line 3229 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (1)].func);
	    }
	}
    break;

  case 622:
/* Line 1792 of yacc.c  */
#line 3238 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (int) (yyvsp[(3) - (4)].ptr)[i];
	    }
	}
    break;

  case 623:
/* Line 1792 of yacc.c  */
#line 3247 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lrand48() % (long) ((yyvsp[(3) - (4)].val));
	    }
	}
    break;

  case 624:
/* Line 1792 of yacc.c  */
#line 3256 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lgamma((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 625:
/* Line 1792 of yacc.c  */
#line 3265 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 626:
/* Line 1792 of yacc.c  */
#line 3274 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log10((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 627:
/* Line 1792 of yacc.c  */
#line 3283 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = 1.0 / (1.0 + exp(-((yyvsp[(3) - (8)].ptr)[i] - (yyvsp[(5) - (8)].val))/ (yyvsp[(7) - (8)].val)));
	    }
	}
    break;

  case 628:
/* Line 1792 of yacc.c  */
#line 3292 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(3) - (6)].ptr)[i] >= (yyvsp[(5) - (6)].ptr)[i] ? (yyvsp[(3) - (6)].ptr)[i] : (yyvsp[(5) - (6)].ptr)[i];
	    }
	}
    break;

  case 629:
/* Line 1792 of yacc.c  */
#line 3301 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(3) - (6)].ptr)[i] <= (yyvsp[(5) - (6)].ptr)[i] ? (yyvsp[(3) - (6)].ptr)[i] : (yyvsp[(5) - (6)].ptr)[i];
	    }
	}
    break;

  case 630:
/* Line 1792 of yacc.c  */
#line 3310 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fmod((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 631:
/* Line 1792 of yacc.c  */
#line 3319 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fx((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 632:
/* Line 1792 of yacc.c  */
#line 3328 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI;
	    }
	}
    break;

  case 633:
/* Line 1792 of yacc.c  */
#line 3337 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI / 180.0;
	    }
	}
    break;

  case 634:
/* Line 1792 of yacc.c  */
#line 3346 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (double) drand48();
	    }
	}
    break;

  case 635:
/* Line 1792 of yacc.c  */
#line 3355 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = rnorm((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 636:
/* Line 1792 of yacc.c  */
#line 3364 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = rnorm((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 637:
/* Line 1792 of yacc.c  */
#line 3373 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = rnorm((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].val));
	    }
	}
    break;

  case 638:
/* Line 1792 of yacc.c  */
#line 3382 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = rnorm((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	    }
	}
    break;

  case 639:
/* Line 1792 of yacc.c  */
#line 3391 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sin((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 640:
/* Line 1792 of yacc.c  */
#line 3400 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(3) - (4)].ptr)[i] * (yyvsp[(3) - (4)].ptr)[i];
	    }
	}
    break;

  case 641:
/* Line 1792 of yacc.c  */
#line 3409 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sqrt((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 642:
/* Line 1792 of yacc.c  */
#line 3418 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = tan((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 643:
/* Line 1792 of yacc.c  */
#line 3426 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
	        if ((int) (yyvsp[(1) - (5)].ptr)[i]) {
		    (yyval.ptr)[i] = (yyvsp[(3) - (5)].ptr)[i];
	        } else {
		    (yyval.ptr)[i] = (yyvsp[(5) - (5)].ptr)[i];
	        }
	    }
	}
    break;

  case 644:
/* Line 1792 of yacc.c  */
#line 3439 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] > (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 645:
/* Line 1792 of yacc.c  */
#line 3448 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] < (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 646:
/* Line 1792 of yacc.c  */
#line 3457 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] <= (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 647:
/* Line 1792 of yacc.c  */
#line 3466 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] >= (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 648:
/* Line 1792 of yacc.c  */
#line 3475 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] == (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 649:
/* Line 1792 of yacc.c  */
#line 3484 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] != (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 650:
/* Line 1792 of yacc.c  */
#line 3493 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] && (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 651:
/* Line 1792 of yacc.c  */
#line 3502 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] || (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 652:
/* Line 1792 of yacc.c  */
#line 3511 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = !((yyvsp[(2) - (2)].ptr)[i]);
	    }
	}
    break;

  case 653:
/* Line 1792 of yacc.c  */
#line 3520 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(2) - (3)].ptr)[i];
	    }
	}
    break;

  case 654:
/* Line 1792 of yacc.c  */
#line 3528 "pars.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = -(yyvsp[(2) - (2)].ptr)[i];
	    }
	}
    break;

  case 656:
/* Line 1792 of yacc.c  */
#line 3539 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (4)].ptr)[(int) (yyvsp[(3) - (4)].val)];
	}
    break;

  case 657:
/* Line 1792 of yacc.c  */
#line 3542 "pars.y"
    {
	    (yyval.val) = ((yyvsp[(3) - (3)].pset) == X) ? *xx : *yy;
	}
    break;

  case 658:
/* Line 1792 of yacc.c  */
#line 3545 "pars.y"
    {
	    double *ptr = getvptr(cg, curset, (yyvsp[(1) - (4)].pset));
	    if (ptr != NULL) {
		(yyval.val) = ptr[(int) (yyvsp[(3) - (4)].val) - 1];
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 659:
/* Line 1792 of yacc.c  */
#line 3555 "pars.y"
    {
	    double *ptr = getvptr(cg, (yyvsp[(1) - (6)].pset), (yyvsp[(3) - (6)].pset));
	    if (ptr != NULL) {
		(yyval.val) = ptr[(int) (yyvsp[(5) - (6)].val) - 1];
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 660:
/* Line 1792 of yacc.c  */
#line 3565 "pars.y"
    {
	    double *ptr = getvptr((yyvsp[(1) - (8)].pset), (yyvsp[(3) - (8)].pset), (yyvsp[(5) - (8)].pset));
	    if (ptr != NULL) {
		(yyval.val) = ptr[(int) (yyvsp[(7) - (8)].val) - 1];
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 661:
/* Line 1792 of yacc.c  */
#line 3575 "pars.y"
    {
	    double *ptr = getvptr(cg, (yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].pset));
	    if (ptr == NULL) {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    switch ((yyvsp[(5) - (5)].pset)) {
	    case MINP:
		(yyval.val) = vmin(ptr, g[cg].p[(yyvsp[(1) - (5)].pset)].len);
		break;
	    case MAXP:
		(yyval.val) = vmax(ptr, g[cg].p[(yyvsp[(1) - (5)].pset)].len);
		break;
	    }
	}
    break;

  case 662:
/* Line 1792 of yacc.c  */
#line 3590 "pars.y"
    {
	    double *ptr = getvptr((yyvsp[(1) - (7)].pset), (yyvsp[(3) - (7)].pset), (yyvsp[(5) - (7)].pset));
	    if (ptr == NULL) {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	    switch ((yyvsp[(7) - (7)].pset)) {
	    case MINP:
		(yyval.val) = vmin(ptr, g[(yyvsp[(1) - (7)].pset)].p[(yyvsp[(3) - (7)].pset)].len);
		break;
	    case MAXP:
		(yyval.val) = vmax(ptr, g[(yyvsp[(1) - (7)].pset)].p[(yyvsp[(3) - (7)].pset)].len);
		break;
	    }
	}
    break;

  case 663:
/* Line 1792 of yacc.c  */
#line 3605 "pars.y"
    {
	    (yyval.val) = g[cg].p[(yyvsp[(1) - (3)].pset)].len;
	}
    break;

  case 664:
/* Line 1792 of yacc.c  */
#line 3609 "pars.y"
    {
	    double bar, sd;
	    double *ptr = getvptr(cg, (yyvsp[(1) - (5)].pset), (yyvsp[(3) - (5)].pset));
	    if (ptr != NULL) {
		stasum(ptr, getsetlength(cg, (yyvsp[(1) - (5)].pset)), &bar, &sd, 0);
	        (yyval.val) = bar;
	    }
	    else {
		yyerror("NULL variable, check set type");
		return 1;
	    }
	}
    break;

  case 665:
/* Line 1792 of yacc.c  */
#line 3621 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].val);
	}
    break;

  case 666:
/* Line 1792 of yacc.c  */
#line 3624 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].val);
	}
    break;

  case 667:
/* Line 1792 of yacc.c  */
#line 3627 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].val);
	}
    break;

  case 668:
/* Line 1792 of yacc.c  */
#line 3631 "pars.y"
    {
	    if ((yyvsp[(3) - (3)].val) != 0.0) {
		(yyval.val) = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].val);
	    } else {
		yyerror("Divide by Zero");
		return 1;
	    }
	}
    break;

  case 669:
/* Line 1792 of yacc.c  */
#line 3639 "pars.y"
    {
	    (yyval.val) = fmod((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val));
	}
    break;

  case 670:
/* Line 1792 of yacc.c  */
#line 3642 "pars.y"
    {
	    (yyval.val) = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val));
	}
    break;

  case 671:
/* Line 1792 of yacc.c  */
#line 3645 "pars.y"
    {
	    (yyval.val) = fabs((yyvsp[(3) - (4)].val));
	}
    break;

  case 672:
/* Line 1792 of yacc.c  */
#line 3648 "pars.y"
    {
	    (yyval.val) = acos((yyvsp[(3) - (4)].val));
	}
    break;

  case 673:
/* Line 1792 of yacc.c  */
#line 3651 "pars.y"
    {
	    (yyval.val) = asin((yyvsp[(3) - (4)].val));
	}
    break;

  case 674:
/* Line 1792 of yacc.c  */
#line 3654 "pars.y"
    {
	    (yyval.val) = atan((yyvsp[(3) - (4)].val));
	}
    break;

  case 675:
/* Line 1792 of yacc.c  */
#line 3657 "pars.y"
    {
	    (yyval.val) = atan2((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 676:
/* Line 1792 of yacc.c  */
#line 3660 "pars.y"
    {
	    (yyval.val) = ceil((yyvsp[(3) - (4)].val));
	}
    break;

  case 677:
/* Line 1792 of yacc.c  */
#line 3663 "pars.y"
    {
	    (yyval.val) = cos((yyvsp[(3) - (4)].val));
	}
    break;

  case 678:
/* Line 1792 of yacc.c  */
#line 3666 "pars.y"
    {
	    (yyval.val) = 180.0 / M_PI;
	}
    break;

  case 679:
/* Line 1792 of yacc.c  */
#line 3669 "pars.y"
    {
	    (yyval.val) = erf((yyvsp[(3) - (4)].val));
	}
    break;

  case 680:
/* Line 1792 of yacc.c  */
#line 3672 "pars.y"
    {
	    (yyval.val) = erfc((yyvsp[(3) - (4)].val));
	}
    break;

  case 681:
/* Line 1792 of yacc.c  */
#line 3675 "pars.y"
    {
	    (yyval.val) = exp((yyvsp[(3) - (4)].val));
	}
    break;

  case 682:
/* Line 1792 of yacc.c  */
#line 3678 "pars.y"
    {
	    (yyval.val) = floor((yyvsp[(3) - (4)].val));
	}
    break;

  case 683:
/* Line 1792 of yacc.c  */
#line 3681 "pars.y"
    {
	    (yyval.val) = my_hypot((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 684:
/* Line 1792 of yacc.c  */
#line 3684 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.xv1;
	}
    break;

  case 685:
/* Line 1792 of yacc.c  */
#line 3687 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.xv2;
	}
    break;

  case 686:
/* Line 1792 of yacc.c  */
#line 3690 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.yv1;
	}
    break;

  case 687:
/* Line 1792 of yacc.c  */
#line 3693 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.yv2;
	}
    break;

  case 688:
/* Line 1792 of yacc.c  */
#line 3696 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.xg1;
	}
    break;

  case 689:
/* Line 1792 of yacc.c  */
#line 3699 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.xg2;
	}
    break;

  case 690:
/* Line 1792 of yacc.c  */
#line 3702 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.yg1;
	}
    break;

  case 691:
/* Line 1792 of yacc.c  */
#line 3705 "pars.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.yg2;
	}
    break;

  case 692:
/* Line 1792 of yacc.c  */
#line 3708 "pars.y"
    {
	    (yyval.val) = g[cg].v.xv1;
	}
    break;

  case 693:
/* Line 1792 of yacc.c  */
#line 3711 "pars.y"
    {
	    (yyval.val) = g[cg].v.xv2;
	}
    break;

  case 694:
/* Line 1792 of yacc.c  */
#line 3714 "pars.y"
    {
	    (yyval.val) = g[cg].v.yv1;
	}
    break;

  case 695:
/* Line 1792 of yacc.c  */
#line 3717 "pars.y"
    {
	    (yyval.val) = g[cg].v.yv2;
	}
    break;

  case 696:
/* Line 1792 of yacc.c  */
#line 3720 "pars.y"
    {
	    (yyval.val) = g[cg].w.xg1;
	}
    break;

  case 697:
/* Line 1792 of yacc.c  */
#line 3723 "pars.y"
    {
	    (yyval.val) = g[cg].w.xg2;
	}
    break;

  case 698:
/* Line 1792 of yacc.c  */
#line 3726 "pars.y"
    {
	    (yyval.val) = g[cg].w.yg1;
	}
    break;

  case 699:
/* Line 1792 of yacc.c  */
#line 3729 "pars.y"
    {
	    (yyval.val) = g[cg].w.yg2;
	}
    break;

  case 700:
/* Line 1792 of yacc.c  */
#line 3732 "pars.y"
    {
	    (yyval.val) = g[cg].w.xg2 - g[cg].w.xg1;
	}
    break;

  case 701:
/* Line 1792 of yacc.c  */
#line 3735 "pars.y"
    {
	    (yyval.val) = g[cg].w.yg2 - g[cg].w.yg1;
	}
    break;

  case 702:
/* Line 1792 of yacc.c  */
#line 3738 "pars.y"
    {
	    (yyval.val) = setindex;
	}
    break;

  case 703:
/* Line 1792 of yacc.c  */
#line 3741 "pars.y"
    {
	    (yyval.val) = setsetno;
	}
    break;

  case 704:
/* Line 1792 of yacc.c  */
#line 3744 "pars.y"
    {
	    (yyval.val) = (long) (yyvsp[(3) - (4)].val);
	}
    break;

  case 705:
/* Line 1792 of yacc.c  */
#line 3747 "pars.y"
    {
	    (yyval.val) = lrand48() % (long) ((yyvsp[(3) - (4)].val));
	}
    break;

  case 706:
/* Line 1792 of yacc.c  */
#line 3750 "pars.y"
    {
	    (yyval.val) = lgamma((yyvsp[(3) - (4)].val));
	}
    break;

  case 707:
/* Line 1792 of yacc.c  */
#line 3753 "pars.y"
    {
	    (yyval.val) = log((yyvsp[(3) - (4)].val));
	}
    break;

  case 708:
/* Line 1792 of yacc.c  */
#line 3756 "pars.y"
    {
	    (yyval.val) = log10((yyvsp[(3) - (4)].val));
	}
    break;

  case 709:
/* Line 1792 of yacc.c  */
#line 3760 "pars.y"
    {
	    (yyval.val) = 1.0 / (1.0 + exp(-((yyvsp[(3) - (8)].val) - (yyvsp[(5) - (8)].val))/ (yyvsp[(7) - (8)].val)));
	}
    break;

  case 710:
/* Line 1792 of yacc.c  */
#line 3763 "pars.y"
    {
	    (yyval.val) = (yyvsp[(3) - (6)].val) >= (yyvsp[(5) - (6)].val) ? (yyvsp[(3) - (6)].val) : (yyvsp[(5) - (6)].val);
	}
    break;

  case 711:
/* Line 1792 of yacc.c  */
#line 3766 "pars.y"
    {
	    (yyval.val) = (yyvsp[(3) - (6)].val) <= (yyvsp[(5) - (6)].val) ? (yyvsp[(3) - (6)].val) : (yyvsp[(5) - (6)].val);
	}
    break;

  case 712:
/* Line 1792 of yacc.c  */
#line 3769 "pars.y"
    {
	    (yyval.val) = fmod((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 713:
/* Line 1792 of yacc.c  */
#line 3772 "pars.y"
    {
	    (yyval.val) = fx((yyvsp[(3) - (4)].val));
	}
    break;

  case 714:
/* Line 1792 of yacc.c  */
#line 3775 "pars.y"
    {
	    (yyval.val) = M_PI;
	}
    break;

  case 715:
/* Line 1792 of yacc.c  */
#line 3778 "pars.y"
    {
	    (yyval.val) = M_PI / 180.0;
	}
    break;

  case 716:
/* Line 1792 of yacc.c  */
#line 3781 "pars.y"
    {
	    (yyval.val) = (double) drand48();
	}
    break;

  case 717:
/* Line 1792 of yacc.c  */
#line 3784 "pars.y"
    {
	    (yyval.val) = rnorm((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 718:
/* Line 1792 of yacc.c  */
#line 3787 "pars.y"
    {
	    (yyval.val) = sin((yyvsp[(3) - (4)].val));
	}
    break;

  case 719:
/* Line 1792 of yacc.c  */
#line 3790 "pars.y"
    {
	    (yyval.val) = pow((yyvsp[(3) - (4)].val), 2.0);
	}
    break;

  case 720:
/* Line 1792 of yacc.c  */
#line 3793 "pars.y"
    {
	    (yyval.val) = sqrt((yyvsp[(3) - (4)].val));
	}
    break;

  case 721:
/* Line 1792 of yacc.c  */
#line 3796 "pars.y"
    {
	    (yyval.val) = tan((yyvsp[(3) - (4)].val));
	}
    break;

  case 722:
/* Line 1792 of yacc.c  */
#line 3799 "pars.y"
    {
	    if ((int) (yyvsp[(1) - (5)].val)) {
		(yyval.val) = (yyvsp[(3) - (5)].val);
	    } else {
		(yyval.val) = (yyvsp[(5) - (5)].val);
	    }
	}
    break;

  case 723:
/* Line 1792 of yacc.c  */
#line 3806 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) > (yyvsp[(3) - (3)].val);
	}
    break;

  case 724:
/* Line 1792 of yacc.c  */
#line 3809 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) < (yyvsp[(3) - (3)].val);
	}
    break;

  case 725:
/* Line 1792 of yacc.c  */
#line 3812 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) <= (yyvsp[(3) - (3)].val);
	}
    break;

  case 726:
/* Line 1792 of yacc.c  */
#line 3815 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) >= (yyvsp[(3) - (3)].val);
	}
    break;

  case 727:
/* Line 1792 of yacc.c  */
#line 3818 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) == (yyvsp[(3) - (3)].val);
	}
    break;

  case 728:
/* Line 1792 of yacc.c  */
#line 3821 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) != (yyvsp[(3) - (3)].val);
	}
    break;

  case 729:
/* Line 1792 of yacc.c  */
#line 3824 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) && (yyvsp[(3) - (3)].val);
	}
    break;

  case 730:
/* Line 1792 of yacc.c  */
#line 3827 "pars.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) || (yyvsp[(3) - (3)].val);
	}
    break;

  case 731:
/* Line 1792 of yacc.c  */
#line 3830 "pars.y"
    {
	    (yyval.val) = !((yyvsp[(2) - (2)].val));
	}
    break;

  case 732:
/* Line 1792 of yacc.c  */
#line 3833 "pars.y"
    {
	    (yyval.val) = (yyvsp[(2) - (3)].val);
	}
    break;

  case 733:
/* Line 1792 of yacc.c  */
#line 3836 "pars.y"
    {
	    (yyval.val) = -(yyvsp[(2) - (2)].val);
	}
    break;


/* Line 1792 of yacc.c  */
#line 11508 "y.tab.c"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2055 of yacc.c  */
#line 3840 "pars.y"


void fixupstr(char *val)
{
    int vl = strlen(val);
    lowtoupper(val);
    val[vl + 1] = 0;
    val[vl] = '\n';
}

void scanner(char *s, double *x, double *y, int len, double *a, double *b, double *c, double *d, int lenscr, int i, int setno, int *errpos)
{
    interr = 0;
    whichgraph = cg;
    whichset = setno;
    if (s[0] == '#') {
	return;
    }
    pos = 0;
    aa = a;
    bb = b;
    cc = c;
    dd = d;
    xx = x;
    yy = y;
    lxy = len;
    ls = lenscr;
    setindex = i + 1;
    curset = setsetno = setno;
    strcpy(f_string, s);
    fcnt = 0;
    yyparse();
    *errpos = interr;
    for (i = 0; i < fcnt; i++) {
	free(freelist[i]);
	freelist[i] = NULL;
    }
}

void runbatch(char *bfile)
{
    double x, y, a, b, c, d;
    int i, setno, errpos, lcnt = 1;
    char stext[256];
    FILE *fp;
    if (strcmp("stdin", bfile)) {
	fp = fopen(bfile, "r");
    }
    else {
	fp = stdin;
    }
    if (fp == NULL) {
        fprintf(stderr, "Error opening batch file \"%s\"\n", bfile);
        exit(1);
    }
    while(fgets(stext, 255, fp) != NULL) {
        if (stext[0] == '#') {
            continue;
        }
	if (strlen(stext) == 0) {
	    continue;
	}
        lowtoupper(stext);
/* TODO check on 0, 0 here for index and setno */
        scanner(stext, &x, &y, 1, ax, bx, cx, dx, 1, 0, 0, &errpos);
        stext[0] = 0;
        if (gotparams && paramfile[0]) {
            if (!getparms(cg, paramfile)) {
            }
            gotparams = 0;
        } else if (gotread && readfile[0]) {
            if (getdata(cg, readfile, readsrc, readtype)) {
                drawgraph();
            }
            gotread = 0;
        }
    }
    if (fp != stdin) {
	fclose(fp);
    }
}

symtab_entry key[] = {
	"A", VAR,
	"ABORT", ABORT,
	"ABOVE", ABOVE,
	"ABS", ABS,
	"ABSOLUTE", ABSOLUTE,
	"ACOS", ACOS,
	"ACTIVATE", ACTIVATE,
	"ACTIVE", ACTIVE,
	"ADP", ADP,
	"ALL", ALL,
	"ALT", ALT,
	"ALTERNATE", ALTERNATE,
	"ALTXAXIS", ALTXAXIS,
	"ALTYAXIS", ALTYAXIS,
	"AND", AND,
	"ANGLE", ANGLE,
	"ANNOTATE", ANNOTATE,
	"APPEND", APPEND,
	"AREA", AREA,
	"ARROW", ARROW,
	"ASIN", ASIN,
	"ATAN", ATAN,
	"ATAN2", ATAN2,
	"AUTO", AUTO,
	"AUTOSCALE", AUTOSCALE,
	"AUTOTICKS", AUTOTICKS,
	"AVG", AVG,
	"AXES", AXES,
	"AXIS", AXIS,
	"B", VAR,
	"BACKBUFFER", BACKBUFFER,
	"BACKGROUND", BACKGROUND,
	"BAR", BAR,
	"BATCH", BATCH,
	"BELOW", BELOW,
	"BIN", BIN,
	"BLOCK", BLOCK,
	"BOTH", BOTH,
	"BOTTOM", BOTTOM,
	"BOX", BOX,
	"BOXPLOT", BOXPLOT,
	"BP", BP,
	"C", VAR,
	"CD", CD,
	"CEIL", CEIL,
	"CELLS", CELLS,
	"CENTER", CENTER,
	"CHAR", CHAR,
	"CHRSTR", CHRSTR,
	"CLEAR", CLEAR,
	"CLICK", CLICK,
	"CMAP", CMAP,
	"CO", COLOR,
	"COLOR", COLOR,
	"COMMENT", COMMENT,
	"COPY", COPY,
	"CORIE", CORIE,
	"COS", COS,
	"CTD", CTD,
	"CYCLE", CYCLE,
	"D", VAR,
	"DAYMONTH", DAYMONTH,
	"DAYOFWEEKL", DAYOFWEEKL,
	"DAYOFWEEKS", DAYOFWEEKS,
	"DAYOFYEAR", DAYOFYEAR,
	"DB", DB,
	"DDMMYY", DDMMYY,
	"DDMONTHSYY", DDMONTHSYY,
	"DDMONTHSYYHHMMSS", DDMONTHSYYHHMMSS,
	"DECIMAL", DECIMAL,
	"DEF", DEF,
	"DEFAULT", DEFAULT,
	"DEG", DEG,
	"DEGREESLAT", DEGREESLAT,
	"DEGREESLON", DEGREESLON,
	"DEGREESMMLAT", DEGREESMMLAT,
	"DEGREESMMLON", DEGREESMMLON,
	"DEGREESMMSSLAT", DEGREESMMSSLAT,
	"DEGREESMMSSLON", DEGREESMMSSLON,
	"DELAY", DELAYP,
	"DELETE", DELETE,
	"DEVICE", DEVICE,
	"DFT", DFT,
	"DIFF", DIFFERENCE,
	"DIFFERENCE", DIFFERENCE,
	"DISK", DISK,
	"DOUBLEBUFFER", DOUBLEBUFFER,
	"DOWN", DOWN,
	"DRAW2", DRAW2,
	"DROP", DROP,
	"DX", DX,
	"DXDX", DXDX,
	"DY", DY,
	"DYDY", DYDY,
	"ECHO", ECHO,
	"EDIT", EDIT,
	"ELCIRC", ELCIRC,
	"ELSE", ELSE,
	"END", END,
	"EQ", EQ,
	"ER", ERRORBAR,
	"ERF", ERF,
	"ERFC", ERFC,
	"ERRORBAR", ERRORBAR,
	"EXIT", EXIT,
	"EXP", EXP,
	"EXPONENTIAL", EXPONENTIAL,
	"FALSE", FALSEP,
	"FEGRID", FEGRID,
	"FFT", FFT,
	"FILE", FILEP,
	"FILL", FILL,
	"FIND", FIND,
	"FIXED", FIXED,
	"FIXEDPOINT", FIXEDPOINT,
	"FLOOR", FLOOR,
	"FLUSH", FLUSH,
	"FOCUS", FOCUS,
	"FOLLOWS", FOLLOWS,
	"FONT", FONTP,
	"FOREGROUND", FOREGROUND,
	"FORMAT", FORMAT,
	"FRAME", FRAMEP,
	"FREE", FREE,
	"FRONTBUFFER", FRONTBUFFER,
	"GE", GE,
	"GENERAL", GENERAL,
	"GETP", GETP,
	"GIFL", GIFL,
	"GIFP", GIFP,
	"GRAPH", GRAPH,
	"GRAPHS", GRAPHS,
	"GRAPHTYPE", GRAPHTYPE,
	"GRID", GRID,
	"GT", GT,
	"HARDCOPY", HARDCOPY,
	"HBAR", HBAR,
	"HBOXPLOT", HBOXPLOT,
	"HEAT", HEAT,
	"HGAP", HGAP,
	"HH", HH,
	"HIDDEN", HIDDEN,
	"HISTO", HISTO,
	"HMS", HMS,
	"HORIZONTAL", HORIZONTAL,
	"HYPOT", HYPOT,
	"IF", IF,
	"IGNORE", IGNORE,
	"IHL", IHL,
	"IMAGE", IMAGE,
	"IN", IN,
	"INDEX", INDEX,
	"INIT", INIT,
	"INITGRAPHICS", INITGRAPHICS,
	"INOUT", INOUT,
	"INT", INT,
	"INTEGRATE", INT,
	"INTERP", INTERP,
	"INUM", INUM,
	"INVDFT", INVDFT,
	"INVFFT", INVFFT,
	"INVN", INVN,
	"INVT", INVT,
	"IRAND", IRAND,
	"ISOLINE", ISOLINE,
	"ISOLINES", ISOLINES,
	"JUST", JUST,
	"KILL", KILL,
	"LABEL", LABEL,
	"LANDSCAPE", LANDSCAPE,
	"LAYOUT", LAYOUT,
	"LE", LE,
	"LEAVE", LEAVE,
	"LEAVEGRAPHICS", LEAVEGRAPHICS,
	"LEFT", LEFT,
	"LEGEND", LEGEND,
	"LENGTH", LENGTH,
	"LEVEL", LEVEL,
	"LEVELS", LEVELS,
	"LGAMMA", LGAMMA,
	"LINE", LINE,
	"LINESTYLE", LINESTYLE,
	"LINETO", LINETO,
	"LINEWIDTH", LINEWIDTH,
	"LINK", LINK,
	"LN", LN,
	"LOAD", LOAD,
	"LOCATOR", LOCATOR,
	"LOCATORBAR", LOCATORBAR,
	"LOCTYPE", LOCTYPE,
	"LOG", LOG,
	"LOGISTIC", LOGISTIC,
	"LOGX", LOGX,
	"LOGXY", LOGXY,
	"LOGY", LOGY,
	"LS", LINESTYLE,
	"LT", LT,
	"LW", LINEWIDTH,
	"MAJOR", MAJOR,
	"MAX", MAXP,
	"MIFL", MIFL,
	"MIFP", MIFP,
	"MIN", MINP,
	"MINOR", MINOR,
	"MISSING", MISSINGP,
	"MMDD", MMDD,
	"MMDDHMS", MMDDHMS,
	"MMDDYY", MMDDYY,
	"MMDDYYHMS", MMDDYYHMS,
	"MMSSLAT", MMSSLAT,
	"MMSSLON", MMSSLON,
	"MMYY", MMYY,
	"MOD", MOD,
	"MONITOR", MONITOR,
	"MONTHDAY", MONTHDAY,
	"MONTHL", MONTHL,
	"MONTHS", MONTHS,
	"MOVE", MOVE,
	"MOVE2", MOVE2,
	"MOVETO", MOVETO,
	"NE", NE,
	"NEGATE", NEGATE,
	"NO", NO,
	"NONE", NONE,
	"NORM", NORM,
	"NORMAL", NORMAL,
	"NORMP", NORMP,
	"NOT", NOT,
	"NUMBER", NUMBER,
	"NXY", NXY,
	"OFF", OFF,
	"OFFSETX", OFFSETX,
	"OFFSETY", OFFSETY,
	"ON", ON,
	"OP", OP,
	"OR", OR,
	"ORIENT", ORIENT,
	"OUT", OUT,
	"PAGE", PAGE,
	"PARA", PARA,
	"PARALLEL", PARALLEL,
	"PARAMETERS", PARAMETERS,
	"PARAMS", PARAMS,
	"PATTERN", PATTERN,
	"PERIMETER", PERIMETER,
	"PERP", PERP,
	"PERPENDICULAR", PERPENDICULAR,
	"PI", PI,
	"PIE", PIE,
	"PIPE", PIPE,
	"PLACE", PLACE,
	"POINT", POINT,
	"POLAR", POLAR,
	"POLYI", POLYI,
	"POLYO", POLYO,
	"POP", POP,
	"PORTRAIT", PORTRAIT,
	"POWER", POWER,
	"PREC", PREC,
	"PREPEND", PREPEND,
	"PRINT", PRINT,
	"PS", PS,
	"PSCOLORL", PSCOLORL,
	"PSCOLORP", PSCOLORP,
	"PSMONOL", PSMONOL,
	"PSMONOP", PSMONOP,
	"PUSH", PUSH,
	"PUTP", PUTP,
	"RAD", RAD,
	"RAND", RAND,
	"READ", READ,
	"RECTGRID", RECTGRID,
	"REDRAW", REDRAW,
	"REGRESS", REGRESS,
	"RENDER", RENDER,
	"REVERSE", REVERSE,
	"RIGHT", RIGHT,
	"RISER", RISER,
	"RNORM", RNORM,
	"ROT", ROT,
	"RUNAVG", RUNAVG,
	"RUNMAX", RUNMAX,
	"RUNMED", RUNMED,
	"RUNMIN", RUNMIN,
	"RUNSTD", RUNSTD,
	"SAMPLE", SAMPLE,
	"SAVEALL", SAVEALL,
	"SCALAR", SCALAR,
	"SCALE", SCALE,
	"SCIENTIFIC", SCIENTIFIC,
	"SET", SET,
	"SETNO", SETNO,
	"SETS", SETS,
	"SIGN", SIGN,
	"SIN", SIN,
	"SIZE", SIZE,
	"SKIP", SKIP,
	"SLEEP", SLEEP,
	"SLICE", SLICE,
	"SOURCE", SOURCE,
	"SPEC", SPEC,
	"SPECIFIED", SPECIFIED,
	"SPECTRUM", SPECTRUM,
	"SPLINE", SPLINE,
	"SQR", SQR,
	"SQRT", SQRT,
	"STACK", STACK,
	"STACKEDBAR", STACKEDBAR,
	"STACKEDHBAR", STACKEDHBAR,
	"STACKEDLINE", STACKEDLINE,
	"STAGGER", STAGGER,
	"START", START,
	"STARTTYPE", STARTTYPE,
	"STATUS", STATUS,
	"STATUSBAR", STATUSBAR,
	"STEP", STEP,
	"STOP", STOP,
	"STRING", STRING,
	"SUBTITLE", SUBTITLE,
	"SWAPBUFFER", SWAPBUFFER,
	"SYMBOL", SYMBOL,
	"TAN", TAN,
	"TICK", TICKP,
	"TICKLABEL", TICKLABEL,
	"TICKMARKS", TICKMARKS,
	"TITLE", TITLE,
	"TO", TO,
	"TOOLBAR", TOOLBAR,
	"TOP", TOP,
	"TPC", TPC,
	"TRUE", TRUEP,
	"TYPE", TYPE,
	"UP", UP,
	"VAR", VAR,
	"VECTOR", VECTOR,
	"VELOCITY", VELOCITY,
	"VERTICAL", VERTICAL,
	"VGAP", VGAP,
	"VIEW", VIEW,
	"VX1", VX1,
	"VX2", VX2,
	"VY1", VY1,
	"VY2", VY2,
	"WITH", WITH,
	"WORLD", WORLD,
	"WRITE", WRITE,
	"WX1", WX1,
	"WX2", WX2,
	"WY1", WY1,
	"WY2", WY2,
	"X", X,
	"X0", X0,
	"X1", X1,
	"XAXES", XAXES,
	"XAXIS", XAXIS,
	"XCOR", XCOR,
	"XMAX", XMAX,
	"XMIN", XMIN,
	"XY", XY,
	"XYARC", XYARC,
	"XYBOX", XYBOX,
	"XYBOXPLOT", XYBOXPLOT,
	"XYDX", XYDX,
	"XYDXDX", XYDXDX,
	"XYDXDY", XYDXDY,
	"XYDY", XYDY,
	"XYDYDY", XYDYDY,
	"XYFIXED", XYFIXED,
	"XYHILO", XYHILO,
	"XYRT", XYRT,
	"XYSEG", XYSEG,
	"XYSTRING", XYSTRING,
	"XYUV", XYUV,
	"XYX2Y2", XYX2Y2,
	"XYXX", XYXX,
	"XYYY", XYYY,
	"XYZ", XYZ,
	"XYZW", XYZW,
	"Y", Y,
	"Y0", Y0,
	"Y1", Y1,
	"Y2", Y2,
	"Y3", Y3,
	"Y4", Y4,
	"Y5", Y5,
	"YAXES", YAXES,
	"YAXIS", YAXIS,
	"YES", YES,
	"YMAX", YMAX,
	"YMIN", YMIN,
	"ZEROXAXIS", ZEROXAXIS,
	"ZEROYAXIS", ZEROYAXIS,
};

int maxfunc = sizeof(key) / sizeof(symtab_entry);

int findf(symtab_entry *key, char *s, int tlen)
{

    int low, high, mid;

    low = 0;
    high = tlen - 1;
    while (low <= high) {
	mid = (low + high) / 2;
	if (strcmp(s, key[mid].s) < 0) {
	    high = mid - 1;
	} else {
	    if (strcmp(s, key[mid].s) > 0) {
		low = mid + 1;
	    } else {
		return (mid);
	    }
	}
    }
    return (-1);
}

int getcharstr(void)
{
    if (pos >= strlen(f_string))
	 return EOF;
    return (f_string[pos++]);
}

void ungetchstr(void)
{
    if (pos > 0)
	pos--;
}

int yylex(void)
{
    int c, i;
    int found;
    static char s[256];
    char sbuf[256];
    char *str;

    while ((c = getcharstr()) == ' ' || c == '\t');
    if (c == EOF) {
	return (0);
    }
    if (c == '"') {
	i = 0;
	while ((c = getcharstr()) != '"' && c != EOF) {
	    if (c == '\\') {
		int ctmp;
		ctmp = getcharstr();
		if (ctmp != '"') {
		    ungetchstr();
		}
		else {
		    c = ctmp;
		}
	    }
	    s[i] = c;
	    i++;
	}
	if (c == EOF) {
	    sprintf(sbuf, "Nonterminating string\n");
	    yyerror(sbuf);
	    return 0;
	}
	s[i] = '\0';
	str = (char *) malloc(strlen(s) + 1);
	strcpy(str, s);
	yylval.str = str;
	return CHRSTR;
    }
    if (c == '.' || isdigit(c)) {
	char stmp[80];
	double d;
	int i, gotdot = 0;

	i = 0;
	while (c == '.' || isdigit(c)) {
	    if (c == '.') {
		if (gotdot) {
		    yyerror("Reading number, too many dots");
	    	    return 0;
		} else {
		    gotdot = 1;
		}
	    }
	    stmp[i++] = c;
	    c = getcharstr();
	}
	if (c == 'E' || c == 'e') {
	    stmp[i++] = c;
	    c = getcharstr();
	    if (c == '+' || c == '-') {
		stmp[i++] = c;
		c = getcharstr();
	    }
	    while (isdigit(c)) {
		stmp[i++] = c;
		c = getcharstr();
	    }
	}
	if (gotdot && i == 1) {
	    ungetchstr();
	    return '.';
	}
	stmp[i] = '\0';
	ungetchstr();
	sscanf(stmp, "%lf", &d);
	yylval.val = d;
	return NUMBER;
    }
/* graphs, sets, regions resp. */
    if (c == 'G' || c == 'S' || c == 'R') {
	char stmp[80];
	double d;
	int i = 0, ctmp = c, gn, sn, rn;
	c = getcharstr();
	while (isdigit(c)) {
	    stmp[i++] = c;
	    c = getcharstr();
	}
	if (i == 0) {
	    c = ctmp;
	    ungetchstr();
	} else {
	    ungetchstr();
	    if (ctmp == 'G') {
	        stmp[i] = '\0';
		gn = atoi(stmp);
		if (gn >= 0 && gn < maxgraph) {
		    yylval.ival = gn;
		    whichgraph = gn;
		    return GRAPHNO;
		}
	    } else if (ctmp == 'S') {
	        stmp[i] = '\0';
		sn = atoi(stmp);
		if (sn >= 0 && sn < g[cg].maxplot) {
		    lxy = getsetlength(cg, sn);
		    yylval.ival = sn;
		    whichset = sn;
		    return SETNUM;
		}
	    } else if (ctmp == 'R') {
	        stmp[i] = '\0';
		rn = atoi(stmp);
		if (rn >= 0 && rn < MAXREGION) {
		    yylval.ival = rn;
		    return REGNUM;
		}
	    }
	}
    }
    if (isalpha(c)) {
	char *p = sbuf;
	int gno = -1, setno = -1, xy = -1, elno = -1;

	do {
	    *p++ = c;
	} while ((c = getcharstr()) != EOF && isalnum(c));
	ungetchstr();
	*p = '\0';
        if (debuglevel == 2) {
	    printf("->%s<-\n", sbuf);
	}
	if ((found = findf(key, sbuf, maxfunc)) >= 0) {
	    if (key[found].type == VAR) {
		switch (sbuf[0]) {
		case 'A':
		    yylval.ptr = aa;
		    return VAR;
		case 'B':
		    yylval.ptr = bb;
		    return VAR;
		case 'C':
		    yylval.ptr = cc;
		    return VAR;
		case 'D':
		    yylval.ptr = dd;
		    return VAR;
		}
	    }
	    else { /* set up special cases */
		switch (key[found].type) {
		case XAXIS:
		    naxis = 0;
		    break;
		case YAXIS:
		    naxis = 1;
		    break;
		case ZEROXAXIS:
		    naxis = 2;
		    break;
		case ZEROYAXIS:
		    naxis = 3;
		    break;
		case ALTXAXIS:
		    naxis = 4;
		    break;
		case ALTYAXIS:
		    naxis = 5;
		    break;
		case AXES:
		    naxis = 6;
		    break;
		case XAXES:
		    naxis = 7;
		    break;
		case YAXES:
		    naxis = 8;
		    break;
		case GRAPHS:
		    yylval.ival = -1;
		    whichgraph = -1;
		    return GRAPHS;
		    break;
		case SETS:
		    yylval.ival = -1;
		    whichset = -1;
		    return SETS;
		    break;
		default:
		    break;
		}
	    }
	    yylval.func = key[found].type;
	    return key[found].type;
	} else {
	    strcat(sbuf, ": No such function or variable");
	    yyerror(sbuf);
	    return 0;
	}
    }
    switch (c) {
    case '>':
	return follow('=', GE, GT);
    case '<':
	return follow('=', LE, LT);
    case '=':
	return follow('=', EQ, '=');
    case '!':
	return follow('=', NE, NOT);
    case '|':
	return follow('|', OR, '|');
    case '&':
	return follow('&', AND, '&');
    case '\n':
	return '\n';
    default:
	return c;
    }
}

int follow(int expect, int ifyes, int ifno)
{
    int c = getcharstr();

    if (c == expect) {
	return ifyes;
    }
    ungetchstr();
    return ifno;
}

void yyerror(char *s)
{
    int i;
    char buf[256];
    sprintf(buf, "Error: %s: %s", s, f_string);
    i = strlen(buf);
    buf[i - 1] = 0;
    errwin(buf);
    interr = 1;
}

#define C1 0.1978977093962766
#define C2 0.1352915131768107

double rnorm(double mean, double sdev)
{
    double u = drand48();

    return mean + sdev * (pow(u, C2) - pow(1.0 - u, C2)) / C1;
}

double fx(double x)
{
    return 1.0 / sqrt(2.0 * M_PI) * exp(-x * x * 0.5);
}

/*
 * return a pointer to the array given by v
 */
double *getvptr(int gno, int setno, int v)
{
    switch (v) {
    case X:
    case X0:
	return g[gno].p[setno].ex[0];
	break;
    case Y:
    case Y0:
	return g[gno].p[setno].ex[1];
	break;
    case Y1:
	return g[gno].p[setno].ex[2];
	break;
    case Y2:
	return g[gno].p[setno].ex[3];
	break;
    case Y3:
	return g[gno].p[setno].ex[4];
	break;
    case Y4:
	return g[gno].p[setno].ex[5];
	break;
    }
    return NULL;
}
