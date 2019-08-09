/****************************************************************
 * fppmain.c
 ****************************************************************/

/******
  Copyright (C) 1993 by Klaus Ehrenfried. 
  
  Permission to use, copy, modify, and distribute this software
  is hereby granted, provided that the above copyright notice appears 
  in all copies and that the software is available to all free of charge. 
  The author disclaims all warranties with regard to this software, 
  including all implied warranties of merchant-ability and fitness. 
  The code is simply distributed as it is.
*******/

#include <stdio.h>
#include <stdlib.h>

#ifdef __EMX__
#include <io.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#define MAIN

#include "fpfli.h"

static int fli_file_created;
static char *fliname;

int exitialise(int error_flag)
{
    if (big_buffer != NULL)
	 free(big_buffer);
    if (pixel_chunk_buffer != NULL)
	free(pixel_chunk_buffer);

    if (output != NULL)
	fclose(output);
    if (input != NULL)
	fclose(input);

    if (error_flag != 0) {
	if (fli_file_created == 1) {
	    fprintf(stderr, "Remove fli-file %s\n", fliname);
	    (void) remove(fliname);
	}
	fprintf(stderr, "abnormal termination\n");
    }
    return (1);
}

/****************************************************************
 * main
 ****************************************************************/

int initialize_fli(char *fname, int nf, int resolution)
{
    FILE *fopen();
    char *listfile, *ppa, *map_file;
    struct stat statbuf;
    char abuff[10];
    int error_flag, answer_flag, big_size, max_chunk_size;
    int i, itest;
    extern int nframes;
    nframes = nf;

    listfile = NULL;
    fliname = NULL;
    map_file = NULL;
    input = 0;
    output = 0;

    fli_file_created = 0;

    big_buffer = NULL;
    pixel_chunk_buffer = NULL;

    error_flag = 0;
    fli_speed = 72;

    Xorigin = 0;
    Xorigin_flag = 0;
    Yorigin = 0;
    Yorigin_flag = 0;

    border_color = 0;
    map_color_flag = 0;
    double_buffer = 1;

    fliname = fname;
    switch (resolution) {
    case 0:
	fli_width = 320;
	fli_height = 200;
	break;
    case 1:
	fli_width = 640;
	fli_height = 400;
	break;
    case 2:
	fli_width = 640;
	fli_height = 480;
	break;
    case 3:
	fli_width = 800;
	fli_height = 600;
	break;
    case 4:
	fli_width = 1024;
	fli_height = 768;
	break;
    case 5:
	fli_width = 1280;
	fli_height = 1024;
	break;
    default:
	fprintf(stderr, "Invalid resolution: -R%d\n", resolution);
	return 0;
    }
    fli_size = fli_width * fli_height;

    if (fli_speed < 0) {
	fprintf(stderr, "Invalid speed: -S%d\n", fli_speed);
	return 0;
    }
    if ((output = fopen(fliname, "wb")) == NULL) {
	fprintf(stderr, "Error opening fli-file %s\n", fliname);
	return 0;
    }
    fli_file_created = 1;

    if (border_color > 0x00FF)
	border_color = 0x00FF;
    if (border_color < 0x0000)
	border_color = 0x0000;

    fprintf(stdout, "Resolution: %dx%d\n", fli_width, fli_height);
    fprintf(stdout, "Origin:     %dx%d\n", Yorigin, Xorigin);
    fprintf(stdout, "Speed:      %d\n", fli_speed);

    max_chunk_size = 8 + fli_height * (2 * fli_width + 10);

    pixel_chunk_buffer = malloc(max_chunk_size);
    if (pixel_chunk_buffer == NULL) {
	fprintf(stderr, "Error: cannot allocate %d bytes\n", max_chunk_size);
	return 0;
    }
    big_size = 5 * fli_size;
    big_buffer = malloc(big_size);
    if (big_buffer == NULL) {
	fprintf(stderr, "Error: cannot allocate %d bytes\n", big_size);
	return 0;
    }
    return 1;
}
