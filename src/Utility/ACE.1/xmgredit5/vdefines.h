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

int triangulate0, sorted, plot, debug;

struct Freenode {
    struct Freenode *nextfree;
};
struct Freelist {
    struct Freenode *head;
    int nodesize;
};
char *getfree();
char *malloc(size_t);
char *myalloc();

double vxmin, vxmax, vymin, vymax, deltax, deltay;


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


struct Site *sites;
int nsites;
int siteidx;
int sqrt_nsites;
int nvertices;
struct Freelist sfl;
struct Site *bottomsite;


struct Edge {
    double a, b, c;
    struct Site *ep[2];
    struct Site *reg[2];
    int edgenbr;
};

#define le 0
#define re 1
int nedges;
struct Freelist efl;

int has_endpoint(), right_of();
struct Site *intersect();
double dist();
struct Point PQ_min();
struct Halfedge *PQextractmin();
struct Edge *bisect();

struct Halfedge {
    struct Halfedge *ELleft, *ELright;
    struct Edge *ELedge;
    int ELrefcnt;
    char ELpm;
    struct Site *vertex;
    double ystar;
    struct Halfedge *PQnext;
};

struct Freelist hfl;
struct Halfedge *ELleftend, *ELrightend;
int ELhashsize;
struct Halfedge **ELhash;
struct Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
struct Site *leftreg(), *rightreg();


int PQhashsize;
struct Halfedge *PQhash;
struct Halfedge *PQfind();
int PQcount;
int PQmin;
int PQempty();
