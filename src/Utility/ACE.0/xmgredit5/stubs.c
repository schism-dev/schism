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
 */

#include "symdefs.h"

void set_writemode(int mode)
{
    setcolor(mode);
}

void get_world(int x, int y, double *wx, double *wy)
{
    device2world(x, y, wx, wy);
}

void diamond(double x, double y)
{
    drawpolysym(&x, &y, 1, SYM_DIAMOND, 0, 0, 0.8);
}

void drawdot(double x, double y)
{
    drawpolysym(&x, &y, 1, SYM_DOT, 0, 1, 0.8);
}

void solidbox(double x, double y)
{
    extern int dosolidbox;
    if (dosolidbox)
    drawpolysym(&x, &y, 1, SYM_SQUARE, 0, 1, 0.8);
}

void box(double x, double y)
{
    drawpolysym(&x, &y, 1, SYM_SQUARE, 0, 0, 0.8);
}

void my_circle(double x, double y)
{
    drawpolysym(&x, &y, 1, SYM_CIRCLE, 0, 0, 0.8);
}

void my_circlefilled(double x, double y)
{
    drawpolysym(&x, &y, 1, SYM_CIRCLE, 0, 1, 0.8);
}

void writetext(void)
{
}

void get_device(double x, double y, int *sx, int *sy)
{
    world2deviceabs(x, y, sx, sy);
}
