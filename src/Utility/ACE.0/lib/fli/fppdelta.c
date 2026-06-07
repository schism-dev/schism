/****************************************************************
 * fppdelta.c
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
#include "fpfli.h"

static int change[FLI_MAX_X];
static int val[FLI_MAX_X];
static int work[FLI_MAX_X];

int
 make_delta_line
 (
   unsigned char *preprevious_line,
   unsigned char *previous_line,
   unsigned char *current_line,
   unsigned char *delta_line,
   int line_skip_count
);

static int new_pixels;

/****************************************************************
 * make_delta_chunk
 ****************************************************************/

int
 make_delta_chunk
 (
   unsigned char *preprevious,	/* pre previous image */
   unsigned char *previous,	/* previous image */
   unsigned char *current	/* current image */
) {
    int skip_count;
    unsigned char *delta_line;
    int chunk_count, line_count;
    int j, line_size, help;
    float change_factor;

    new_pixels = 0;
    chunk_count = 8;		/* 4 bytes for total size of chunk (header) */
    /* 2 bytes for chunk type (header) */
    /* 2 bytes for number of compressed lines */

    skip_count = 0;
    line_count = 0;

    for (j = 0; j < fli_height; j++) {
	delta_line = &pixel_chunk_buffer[chunk_count];

	/*
	 * printf(" Line: %d   chunk_count: %d   skip: %d\n",
	 * j,chunk_count,skip_count); 
	 */

	line_size = make_delta_line(preprevious, previous, current,
				    delta_line, skip_count);

	preprevious += fli_width;
	previous += fli_width;
	current += fli_width;

	/* printf(" Line: %d  size: %d\n",j,line_size); */

	if (line_size > 0) {	/* yes, we got a new line */
	    chunk_count += line_size;
	    skip_count = 0;
	    line_count++;
	} else {
	    skip_count--;
	}
    }

    if (line_count == 0)	/* no lines no chunk */
	return (0);

    if ((chunk_count % 2) == 1)
	add_bytes(pixel_chunk_buffer, &chunk_count, 0x0000, IOM_UBYTE);

    help = 0;
    add_bytes(pixel_chunk_buffer, &help, chunk_count, IOM_LONG);
    add_bytes(pixel_chunk_buffer, &help, FLI_DELTA, IOM_UWORD);
    add_bytes(pixel_chunk_buffer, &help, line_count, IOM_UWORD);

    change_factor = 100.0 * ((float) new_pixels) / (fli_width * fli_height);

    fprintf(stdout, " Delta chunk: %d bytes    new pixels: %.2f %%\n",
	    chunk_count, change_factor);

    return (chunk_count);
}

/****************************************************************
 * make_delta_line
 ****************************************************************/

int
 make_delta_line
 (
   unsigned char *preprevious_line,
   unsigned char *previous_line,
   unsigned char *current_line,
   unsigned char *delta_line,
   int line_skip_count
) {
    int skip_count, size_count, packets;
    int i, m, ipos, extra, help;
    int w0, w1, w2, wpos0, wpos1, wpos2, wsum, improve_count;

    for (i = 0; i < fli_width; i++) {
	val[i] = current_line[i];
	if ((previous_line[i] != current_line[i]) ||
	((double_buffer == 1) && (preprevious_line[i] != current_line[i]))) {
	    change[i] = 1;	/* yes, update */
	} else {
	    change[i] = 0;	/* no update */
	}
    }

    if ((change[fli_width - 2] != 0) || (change[fli_width - 1] != 0)) {
	work[fli_width - 2] = 1;
	new_pixels += 2;
    } else {
	work[fli_width - 2] = 0;
    }
    work[fli_width - 1] = 0;

    for (i = (fli_width - 4); i >= 0; i -= 2) {
	if ((change[i] != 0) || (change[i + 1] != 0)) {
	    new_pixels += 2;
	    if ((val[i] == val[i + 2]) && (val[i + 1] == val[i + 3])) {
		if (work[i + 2] < 0) {
		    work[i] = work[i + 2] - 1;
		    if (work[i] < -127)
			work[i] = -1;
		} else if (work[i + 2] == 0) {
		    work[i] = 1;
		} else {
		    work[i + 2] = -1;
		    work[i] = -2;
		}
	    } else {		/* count nonequal words */
		if (work[i + 2] > 0) {
		    work[i] = work[i + 2] + 1;
		    if (work[i] > 127)
			work[i] = 1;
		} else {
		    work[i] = 1;
		}
	    }
	} else {
	    work[i] = 0;
	}
	work[i + 1] = 0;
    }

    improve_count = 1;

    /* printf(" ----------------------------\n"; */

    while (improve_count > 0) {
	/* printf(" xxxxxxxxxxx > %d\n",improve_count); */

	improve_count = 0;
	w1 = 0;
	wpos1 = 0;
	w2 = 0;
	wpos2 = 0;
	i = 0;

	while (i < fli_width) {
	    if (work[i] != 0) {
		w0 = w1;
		wpos0 = wpos1;
		w1 = w2;
		wpos1 = wpos2;
		w2 = work[i];
		wpos2 = i;

		if (w2 < 0)
		    i -= 2 * w2;
		else
		    i += 2 * w2;

		/* printf(" %d   %d   %d  --> %d\n",w0,w1,w2,i); */

		if (((w2 == -2) || (w2 == -3)) && (w1 > 0)) {
		    wsum = w1 - w2;
		    if (wsum > 127)
			continue;
		    work[wpos1] = wsum;
		    wpos2 = wpos1;
		    w2 = wsum;
		    w1 = 0;
		    improve_count++;
		} else if ((w2 > 0) && ((w1 == -2) || (w1 == -3))) {
		    wsum = w2 - w1;
		    if (wsum > 127)
			continue;
		    work[wpos1] = wsum;
		    wpos2 = wpos1;
		    w2 = wsum;
		    w1 = 0;
		    improve_count++;
		}
		/****
		else if ((w0 > 0) && (w1 == -3) && (w2 > 0))
		{
		    wsum=w0+3+w2;
		    if (wsum > 127) continue;
		    work[wpos0]=wsum;
		    wpos2=wpos0;
		    w2=wsum;
		    w1=0;
		    improve_count++;
		}
		*****/
		else if ((w2 > 0) && (w1 > 0)) {
		    wsum = w2 + w1;
		    if (wsum > 127)
			continue;
		    work[wpos1] = wsum;
		    wpos2 = wpos1;
		    w2 = wsum;
		    w1 = 0;
		    improve_count++;
		} else if (((w2 == -2) || (w2 == -3)) &&
			   ((w1 == -2) || (w1 == -3))) {
		    wsum = -w1 - w2;
		    work[wpos1] = wsum;
		    wpos2 = wpos1;
		    w2 = wsum;
		    w1 = 0;
		    improve_count++;
		}
	    } else {
		w1 = 0;
		w2 = 0;
		i += 2;
	    }
	}
    }

    packets = 0;
    skip_count = 0;
    extra = 0;
    i = 0;
    if (line_skip_count != 0)
	ipos = 4;
    else
	ipos = 2;

    while (i < fli_width) {	/* assemble output */
	if (work[i] != 0) {	/* add data packet */
	    packets++;
	    extra = 0;
	    size_count = work[i];
	    add_bytes(delta_line, &ipos, skip_count, IOM_UBYTE);
	    add_bytes(delta_line, &ipos, size_count, IOM_SBYTE);
	    if (size_count < 0) {
		m = i;
		add_bytes(delta_line, &ipos, val[i], IOM_UBYTE);
		add_bytes(delta_line, &ipos, val[i + 1], IOM_UBYTE);
		i -= 2 * size_count;
	    } else {
		for (m = 0; m < size_count; m++) {
		    add_bytes(delta_line, &ipos, val[i++], IOM_UBYTE);
		    add_bytes(delta_line, &ipos, val[i++], IOM_UBYTE);
		}
	    }
	    skip_count = 0;
	} else {
	    /* insert extra skip packet to keep */
	    /* skip_count below 256 */
	    skip_count++;
	    i++;
	    if ((skip_count > 254) && (i < fli_width)) {
		if (work[i] == 0) {
		    packets++;
		    extra++;
		    i--;
		    skip_count--;
		    size_count = 1;
		    add_bytes(delta_line, &ipos, skip_count, IOM_UBYTE);
		    add_bytes(delta_line, &ipos, size_count, IOM_SBYTE);
		    add_bytes(delta_line, &ipos, val[i++], IOM_UBYTE);
		    add_bytes(delta_line, &ipos, val[i++], IOM_UBYTE);
		    skip_count = 0;
		}
	    }
	}
    }

    /* remove extra skips if no packet followed */
    while (extra-- > 0) {
	packets--;
	ipos -= 4;
    }


    if (packets == 0)
	return (0);		/* no data packets --> no line */

    help = 0;
    if (line_skip_count != 0) {
	add_bytes(delta_line, &help, line_skip_count, IOM_SWORD);
	add_bytes(delta_line, &help, packets, IOM_SWORD);
    } else {
	add_bytes(delta_line, &help, packets, IOM_SWORD);
    }

    return (ipos);		/* return number of bytes */
}
