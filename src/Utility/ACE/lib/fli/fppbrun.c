/****************************************************************
 * fppbrun.c
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
#include "fpfli.h"

static int work[FLI_MAX_X];
static int val[FLI_MAX_X];

int make_brun_line(unsigned char *image_line,
		    unsigned char *brun_line);

/****************************************************************
 * make_brun_chunk
 ****************************************************************/

int
 make_brun_chunk
 (
   UBYTE * image		/* first image */
) {
    int chunk_count, j, help;
    unsigned char *brun_line;
    float compression;

    chunk_count = 6;		/* 4 bytes for chunk size */
    /* 2 bytes for chunk type */
    for (j = 0; j < fli_height; j++) {
	brun_line = &pixel_chunk_buffer[chunk_count];
	chunk_count += make_brun_line(image, brun_line);
	image += fli_width;
    }

    if ((chunk_count % 2) == 1)
	add_bytes(pixel_chunk_buffer, &chunk_count, 0x0000, IOM_UBYTE);

    help = 0;
    add_bytes(pixel_chunk_buffer, &help, chunk_count, IOM_LONG);
    add_bytes(pixel_chunk_buffer, &help, FLI_BRUN, IOM_UWORD);

    compression = fli_size / ((float) chunk_count);

    fprintf(stdout, " Brun chunk: %d bytes    compression: %f\n",
	    chunk_count, compression);

    return (chunk_count);
}

/****************************************************************
 * make_brun_line
 ****************************************************************/

int make_brun_line
 (
   unsigned char *image_line,
   unsigned char *brun_line
) {
    int i, j, ipos, packets;
    int size_count, help;
    int mark_pos, mark_jump, min_size, inext;
    int print_flag;

    for (i = 0; i < fli_width; i++)
	val[i] = image_line[i];

    work[fli_width - 1] = -1;

    for (i = (fli_width - 2); i >= 0; i--) {
	if (val[i] == val[i + 1]) {
	    if (work[i + 1] > 0) {
		work[i] = work[i + 1] + 1;
		if (work[i] > 127)
		    work[i] = 1;
	    } else {
		work[i + 1] = 1;
		work[i] = 2;
	    }
	} else {
	    if (work[i + 1] < 0) {
		work[i] = work[i + 1] - 1;
		if (work[i] < -127)
		    work[i] = -1;
	    } else {
		work[i] = -1;
	    }
	}
    }

    print_flag = 0;
    for (min_size = 2; min_size <= 256; min_size++) {	/* check number of
							 * packets and */
	/* reduce it if necessary */
	packets = 0;
	for (i = 0; i < fli_width; i++) {
	    if ((work[i] == 1) || (work[i] == -1))
		packets++;
	}

	if (print_flag == 1) {
	    fprintf(stdout, " reduced to %d packets\n", packets);
	}
	if (packets < 256)
	    break;

	print_flag = 1;
	fprintf(stdout, " Too many packets: %d ... ", packets);

	i = 0;
	mark_pos = -1;
	mark_jump = fli_width;

	while (i < fli_width) {	/* search small packets and *//* merge them
				 * together */
	    inext = (work[i] > 0) ? (i + work[i]) : (i - work[i]);
	    if ((work[i] < -min_size) || (work[i] > min_size)) {	/* big packet */
		/* test if more than one small packet was before this */
		if ((mark_pos >= 0) && (mark_jump < i)) {
		    help = -1;
		    for (j = (i - 1); j >= mark_pos; j--) {
			work[j] = help--;	/* merge */
			if (help < -127)
			    help = -1;
		    }
		}
		mark_pos = -1;
		mark_jump = fli_width;
	    } else if (mark_pos < 0) {
		mark_pos = i;	/* here starts a small packet */
		mark_jump = inext;
	    }
	    i = inext;
	}

	if ((mark_pos >= 0) && (mark_jump < fli_width)) {
	    help = -1;
	    for (j = (fli_width - 1); j >= mark_pos; j--) {
		work[j] = help--;	/* merge */
		if (help < -127)
		    help = -1;
	    }
	}
    }

    ipos = 1;

    packets = 0;
    i = 0;
    while (i < fli_width) {
	size_count = work[i];
	/* fprintf(stdout," %d  %d\n",i,size_count); */

	add_bytes(brun_line, &ipos, size_count, IOM_UBYTE);
	if (size_count > 0) {
	    add_bytes(brun_line, &ipos, val[i], IOM_UBYTE);
	    i += size_count;
	} else {
	    help = i - size_count;
	    while (i < help) {
		add_bytes(brun_line, &ipos, val[i++], IOM_UBYTE);
	    }
	}
	packets++;
    }
    /* fprintf(stdout," packets: %d    ipos: %d\n\n",packets,ipos); */

    help = 0;
    add_bytes(brun_line, &help, packets, IOM_UBYTE);

    return (ipos);
}
