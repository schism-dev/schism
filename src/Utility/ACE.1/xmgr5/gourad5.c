/*
   Shaded Triangle for GD

	Implementation based on the 3DICA3 tutorial
	http://tfpsly.free.fr/Docs/3dIca/3dica3.htm
	License New BSD (See gd CVS for a copy of the license)
	(c) 2007 Pierre-Alain Joye
*/


#include <gd.h>
#include <stdio.h>
#include <string.h>

/* That will end in fast ops anyway (at least vc and gcc :) */
#define SWAP(a, b) do {\
  char c[sizeof(a)]; \
  memcpy((void *)&c, (void *)&a, sizeof(c)); \
  memcpy((void *)&a, (void *)&b, sizeof(a)); \
  memcpy((void *)&b, (void *)&c, sizeof(b)); \
} while (0)

#define GD_SHADED_FP 16

void _gdImageShadedHLine_FP(gdImagePtr im, int x1, int x2, int y,
	int c1_r, int c1_g, int c1_b,
	int c2_r, int c2_g, int c2_b)
{
	long d_r, d_g, d_b;
	long c_r, c_g, c_b;
#define RETAIN_SPATIAL_PRECISION 1
#if RETAIN_SPATIAL_PRECISION
    long dval;
#endif
    long x1_cast, x2_cast;
	int i;

	if(x1==x2) {
		 return;
	}

	if(x1>x2) {
		SWAP(x1, x2);
		SWAP(c1_r, c2_r);
		SWAP(c1_g, c2_g);
		SWAP(c1_b, c2_b);
	}

#if RETAIN_SPATIAL_PRECISION
    /* Compute color at rounded values */
    x1_cast = ((x1 + (1<<(GD_SHADED_FP-1))) >> GD_SHADED_FP) << GD_SHADED_FP;
    x2_cast = ((x2 + (1<<(GD_SHADED_FP-1))) >> GD_SHADED_FP) << GD_SHADED_FP;
    dval = (c2_r - c1_r);
    c1_r = dval * (x1_cast - x1) / (x2 - x1) + c1_r;
    c2_r = dval * (x2_cast - x2) / (x2 - x1) + c2_r;
    dval = (c2_g - c1_g);
    c1_g = dval * (x1_cast - x1) / (x2 - x1) + c1_g;
    c2_g = dval * (x2_cast - x2) / (x2 - x1) + c2_g;
    dval = (c2_b - c1_b);
    c1_b = dval * (x1_cast - x1) / (x2 - x1) + c1_b;
    c2_b = dval * (x2_cast - x2) / (x2 - x1) + c2_b;
#endif

	d_r = (c2_r - c1_r) / ((x2-x1) >> GD_SHADED_FP);
	d_g = (c2_g - c1_g) / ((x2-x1) >> GD_SHADED_FP);
	d_b = (c2_b - c1_b) / ((x2-x1) >> GD_SHADED_FP);

	c_r = c1_r;
	c_g = c1_g;
	c_b = c1_b;

    x1 = (x1 + (1<<(GD_SHADED_FP-1))) >> GD_SHADED_FP;
    x2 = (x2 + (1<<(GD_SHADED_FP-1))) >> GD_SHADED_FP;

	i=0;
	while((i++) < (x2 - x1)) {
#if RETAIN_SPATIAL_PRECISION
		gdImageSetPixel(im, x1 + i, y,
		gdTrueColorAlpha((c_r<0?0:c_r >> GD_SHADED_FP)&0xFF,
		(c_g<0?0:c_g >> GD_SHADED_FP)&0xFF,
		(c_b<0?0:c_b >> GD_SHADED_FP)&0xFF,0));
#else
		gdImageSetPixel(im, x1 + i, y,
		gdTrueColorAlpha((c_r >> GD_SHADED_FP)&0xFF,
		(c_g >> GD_SHADED_FP)&0xFF,
		(c_b >> GD_SHADED_FP)&0xFF,0));
#endif
		c_r += d_r;
		c_g += d_g;
		c_b += d_b;
	}
}

void gdImageShadedTriangle(gdImagePtr im, gdPoint v[3], int c[3])
{
	int x1, y1, x2, y2, x3, y3;
	int c1, c2, c3;

	/*
	 sl == scanline dc == delta color dx == delta x
	 */
	long dx12,dx13,dx23,sl_x1,sl_x2;
	long dc12_r, dc13_r, dc23_r, sl_c1_r, sl_c2_r;
	long dc12_g, dc13_g, dc23_g, sl_c1_g, sl_c2_g;
	long dc12_b, dc13_b, dc23_b, sl_c1_b, sl_c2_b;
	int i;

	x1 = v[0].x; y1 = v[0].y;
	x2 = v[1].x; y2 = v[1].y;
	x3 = v[2].x; y3 = v[2].y;
	c1 = c[0]; c2 = c[1]; c3 = c[2];

	/* Sort the vertices */
	if(y1 > y2) {
		SWAP(y1, y2);
		SWAP(x1, x2);
		SWAP(c1, c2);
	}

	if(y2 > y3) {
		SWAP(y2, y3);
		SWAP(x2, x3);
		SWAP(c2, c3);
	}

	if(y1 > y2) {
		SWAP(y2, y1);
		SWAP(x2, x1);
		SWAP(c2, c1);
	}

	if((i= y2 - y1) != 0) {
		dx12 = ((x2 - x1) << GD_SHADED_FP) / i;
		dc12_r = ((gdImageRed(im, c2) - gdImageRed(im, c1)) << GD_SHADED_FP) / i;
		dc12_g = ((gdImageGreen(im, c2) - gdImageGreen(im, c1)) << GD_SHADED_FP) / i;
		dc12_b = ((gdImageBlue(im, c2) - gdImageBlue(im, c1)) << GD_SHADED_FP) / i;
	} else {
		dx12 = 0;
		dc12_r = 0;
		dc12_g = 0;
		dc12_b = 0;
	}

	if((i = y3 - y1) != 0) {
		dx13 = ((x3-x1) << GD_SHADED_FP) / i;
		dc13_r = ((gdImageRed(im, c3) - gdImageRed(im, c1)) << GD_SHADED_FP) / i;
		dc13_g = ((gdImageGreen(im, c3) - gdImageGreen(im, c1)) << GD_SHADED_FP) / i;
		dc13_b = ((gdImageBlue(im, c3) - gdImageBlue(im, c1)) << GD_SHADED_FP) / i;
	} else {
		dx13 = 0;
		dc13_r = 0;
		dc13_g = 0;
		dc13_b = 0;
	}

	if((i = y3 - y2) != 0) {
		dx23 = ((x3-x2) << GD_SHADED_FP) / i;
		dc23_r = ((gdImageRed(im, c3) - gdImageRed(im, c2)) << GD_SHADED_FP) / i;
		dc23_g = ((gdImageGreen(im, c3) - gdImageGreen(im, c2)) << GD_SHADED_FP) / i;
		dc23_b = ((gdImageBlue(im, c3) - gdImageBlue(im, c2)) << GD_SHADED_FP) / i;
	} else {
		dx23 = 0;
		dc23_r = 0;
		dc23_g = 0;
		dc23_b = 0;
	}

	sl_x1 = sl_x2 = x1 << GD_SHADED_FP;
	sl_c1_r = sl_c2_r = gdImageRed(im, c1) << GD_SHADED_FP;
	sl_c1_g = sl_c2_g = gdImageGreen(im, c1) << GD_SHADED_FP;
	sl_c1_b = sl_c2_b = gdImageBlue(im, c1) << GD_SHADED_FP;

	for(i = y1; i < y2; i++) {
		_gdImageShadedHLine_FP(im,
			sl_x1, sl_x2, i,
			sl_c1_r, sl_c1_g, sl_c1_b,
			sl_c2_r, sl_c2_g, sl_c2_b
		);

		sl_x1 += dx13;
		sl_x2 += dx12;

		sl_c1_r += dc13_r;
		sl_c1_g += dc13_g;
		sl_c1_b += dc13_b;
		sl_c2_r += dc12_r;
		sl_c2_g += dc12_g;
		sl_c2_b += dc12_b;
	}

	sl_x2 = x2 << GD_SHADED_FP;
	sl_c2_r = gdImageRed(im, c2) << GD_SHADED_FP;
	sl_c2_g = gdImageGreen(im, c2) << GD_SHADED_FP;
	sl_c2_b = gdImageBlue(im, c2) << GD_SHADED_FP;

	for(i = y2; i < y3; i++) {
		_gdImageShadedHLine_FP(im,
		sl_x1, sl_x2, i,
		sl_c1_r, sl_c1_g, sl_c1_b,
		sl_c2_r, sl_c2_g, sl_c2_b
		);

		sl_x1+=dx13;
		sl_x2+=dx23;
		sl_c1_r += dc13_r;
		sl_c1_g += dc13_g;
		sl_c1_b += dc13_b;
		sl_c2_r += dc23_r;
		sl_c2_g += dc23_g;
		sl_c2_b += dc23_b;
	}
}

/*
int main()
{
	gdImagePtr im;
	int error = 0;
	FILE *fp;

	gdPoint v[3];
	int c[3];

93 413
97 413
97 404


	v[0].x = 10;
	v[0].y = 10;
	c[0] = gdTrueColorAlpha(255, 0, 0, 0);

	v[1].x = 310;
	v[1].y = 10;
	c[1] = gdTrueColorAlpha(0, 255, 0, 0);

	v[2].x = 160;
	v[2].y = 270;
	c[2] = gdTrueColorAlpha(0, 0, 255, 0);

	im = gdImageCreateTrueColor(320, 280);

	gdImageShadedTriangle(im, v, c);


	fp = fopen("gouraud5.png", "wb");
 	gdImagePng(im,fp);
	fclose(fp);

 	gdImageDestroy(im);
	return error;
}
*/
