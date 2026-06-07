
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !												!
 * !	Inverse bilinear mapping for quadrangles						!
 * !	Convexity of the quad must have been checked, and (x,y) must not be outside the quad.	!
 * !												!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

#include <stdio.h>
#include <math.h>

int ibilinear(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double x, double y, double *xi, double *eta, double *shapef)
{
    double axi[2], aet[2], bxy[2], root_xi[2], root_et[2];
    double x0, y0, dxi, deta, dd, beta, delta, dgamma;
    int i, j, icount, icaseno;

    static double SMALL = 1.0e-5;

//!...  Consts.
    x0 = (x1 + x2 + x3 + x4) / 4.0;
    y0 = (y1 + y2 + y3 + y4) / 4.0;
    axi[0] = x2 - x1 + x3 - x4;
    axi[1] = y2 - y1 + y3 - y4;
    aet[0] = x3 + x4 - x1 - x2;
    aet[1] = y3 + y4 - y1 - y2;
    bxy[0] = x1 - x2 + x3 - x4;
    bxy[1] = y1 - y2 + y3 - y4;
    dxi = 2 * ((x3 - x4) * (y1 - y2) - (y3 - y4) * (x1 - x2));
    deta = 2 * ((x4 - x1) * (y3 - y2) - (y4 - y1) * (x3 - x2));

//!...  Inverse mapping
    if ((fabs(bxy[0]) < SMALL && fabs(bxy[1]) < SMALL) || (fabs(dxi) < SMALL && fabs(deta) < SMALL)) {
	icaseno = 1;
//    if (dxi == 0.0 && deta == 0.0) {
	dd = axi[0] * aet[1] - axi[1] * aet[0];
	if (dd == 0.0) {
	    fprintf(stderr, "Case 1 error: %lf\n", dd);
	    return 1;
	}
	*xi = 4 * (aet[1] * (x - x0) - aet[0] * (y - y0)) / dd;
	*eta = 4 * (axi[0] * (y - y0) - axi[1] * (x - x0)) / dd;

    } else if (fabs(dxi) < SMALL && fabs(deta) >= SMALL) {
//    } else if (dxi == 0 && deta != 0) {
	icaseno = 2;
	*eta = 4 * (bxy[1] * (x - x0) - bxy[1] * (y - y0)) / deta;
	dd = (axi[0] + *eta * bxy[0]) * (axi[0] + *eta * bxy[0]) + (axi[1] + *eta * bxy[1]) * (axi[1] + *eta * bxy[1]);
	if (dd == 0) {
	    fprintf(stderr, "Case 2 error: &lf\n", dd);
	    return 1;
	}
	*xi = ((4 * (x - x0) - *eta * aet[0]) * (axi[0] + *eta * bxy[0]) + (4 * (y - y0) - *eta * aet[1]) * (axi[1] + *eta * bxy[1])) / dd;
	;
    } else if (fabs(dxi) >= SMALL && fabs(deta) < SMALL) {
	icaseno = 3;
//    } else if (dxi != 0 && deta == 0) {
	*xi = 4 * (bxy[1] * (x - x0) - bxy[0] * (y - y0)) / dxi;
	dd = (aet[0] + *xi * bxy[0]) * (aet[0] + *xi * bxy[0]) + (aet[1] + *xi * bxy[1]) * (aet[1] + *xi * bxy[1]);
	if (dd == 0) {
	    fprintf(stderr, "Case 3 error: &lf\n", dd);
	    return 1;
	}
	*eta = ((4 * (x - x0) - *xi * axi[0]) * (aet[0] + *xi * bxy[0]) + (4 * (y - y0) - *xi * axi[1]) * (aet[1] + *xi * bxy[1])) / dd;
	;
    } else {
	icaseno = 4;
	beta = aet[1] * axi[0] - aet[0] * axi[1] - 4 * (bxy[1] * (x - x0) - bxy[0] * (y - y0));
	dgamma = 4 * (aet[0] * (y - y0) - aet[1] * (x - x0));
	delta = beta * beta - 4 * dgamma * dxi;
	if (delta == 0) {
	    *xi = (-beta / 2.0 / dxi);
	    *eta = (4 * (bxy[1] * (x - x0) - bxy[0] * (y - y0)) - *xi * dxi) / deta;
	} else if (delta > 0.0) {
	    root_xi[0] = (-beta + sqrt(delta)) / 2 / dxi;
	    root_xi[1] = (-beta - sqrt(delta)) / 2 / dxi;
	    icount = 0;
	    for (i = 0; i < 2; i++) {
		root_et[i] = (4 * (bxy[1] * (x - x0) - bxy[0] * (y - y0)) - root_xi[i] * dxi) / deta;
		if (fabs(root_xi[i]) <= 1.1 && fabs(root_et[i]) <= 1.1) {
		    *xi = root_xi[i];
		    *eta = root_et[i];
		    icount = icount + 1;
		}
	    }
	    if (icount == 2 && fabs(root_xi[0] - root_xi[1]) < SMALL) {
	    } else if (icount != 1) {
		fprintf(stderr, "Abnormal instances %lf %lf %d\n", root_xi[0], root_et[1], icount);
		return 1;
	    }
	} else {
	    fprintf(stderr, "No roots %lf\n", delta);
	    return 1;
	}
    }

    if (fabs(*xi) > 1.1 || fabs(*eta) > 1.1) {
	fprintf(stderr, "Out of bound in ibilinear %lf %lf\n", *xi, *eta);
	return 1;
    }

    *xi = (*xi > -1.0) ? *xi : -1.0;
    *eta = (*eta > -1.0) ? *eta : -1.0;
    *xi = (*xi < 1.0) ? *xi : 1.0;
    *eta = (*eta < 1.0) ? *eta : 1.0;
    //xi=dmin1(1.0,dmax1(xi,-1.0));
    //eta=dmin1(1.0,dmax1(eta,-1.0));

    shapef[0] = (1 - *xi) * (1 - *eta) / 4;
    shapef[1] = (1 + *xi) * (1 - *eta) / 4;
    shapef[2] = (1 + *xi) * (1 + *eta) / 4;
    shapef[3] = (1 - *xi) * (1 + *eta) / 4;
    return 0;
}
