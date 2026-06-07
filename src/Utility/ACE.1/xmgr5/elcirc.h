/*
 * include for elcirc
 */

#ifndef ELCIRC_H
#define ELCIRC_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

typedef struct _Elcirc {
    int m_type; // 1 or 2 also nitems
    char m_file[2048];
    char m_rundes[32];
    char m_runid[24];
    char m_agrid[24]; // 80 bytes
    int m_ndsetse;
    int m_np;
    float m_dtnspoolge;
    int m_nspoolge;
    int m_irtype; // 20 bytes
    int m_nsteps;
    int m_nlevels;
    FILE *m_fp;
    float n;
    float *t;
    int *it;
    float *d;
} Elcirc;

#endif
