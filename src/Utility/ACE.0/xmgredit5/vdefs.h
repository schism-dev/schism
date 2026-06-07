/*
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.- */

#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

extern int triangulate0, sorted, plot, debug;

struct Freenode {
    struct Freenode *nextfree;
};
struct Freelist {
    struct Freenode *head;
    int nodesize;
};
char *getfree(struct Freelist * fl);
char *myalloc(unsigned int n);

extern double vxmin, vxmax, vymin, vymax, deltax, deltay;


struct Point {
    double x, y;
};

/* structure used both for sites and for vertices */
struct Site {
    struct Point coord;
    int sitenbr;
    int refcnt;
    int ind;
};


extern struct Site *sites;
extern int nsites;
extern int siteidx;
extern int sqrt_nsites;
extern int nvertices;
extern struct Freelist sfl;
extern struct Site *bottomsite;


struct Edge {
    double a, b, c;
    struct Site *ep[2];
    struct Site *reg[2];
    int edgenbr;
};

#define le 0
#define re 1
extern int nedges;
extern struct Freelist efl;

int has_endpoint(), right_of();
struct Site *intersect();
double dist(struct Site * s, struct Site * t);
struct Point PQ_min(void);
struct Halfedge *PQextractmin(void);
struct Edge *bisect(struct Site * s1, struct Site * s2);

struct Halfedge {
    struct Halfedge *ELleft, *ELright;
    struct Edge *ELedge;
    int ELrefcnt;
    char ELpm;
    struct Site *vertex;
    double ystar;
    struct Halfedge *PQnext;
};

extern struct Freelist hfl;
extern struct Halfedge *ELleftend, *ELrightend;
extern int ELhashsize;
extern struct Halfedge **ELhash;
struct Halfedge *HEcreate(struct Edge * e, int pm), *ELleft(struct Halfedge * he), *ELright(struct Halfedge * he), *ELleftbnd(struct Point * p);
struct Site *leftreg(struct Halfedge * he), *rightreg(struct Halfedge * he);


extern int PQhashsize;
extern struct Halfedge *PQhash;
struct Halfedge *PQfind();
extern int PQcount;
extern int PQmin;
int PQempty(void);
