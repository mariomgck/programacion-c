/*
 *Copyright 2002
 *Center for Biological and Computational Learning at MIT and MIT
 *All rights reserved.
 *Permission to copy this software, and its documentation only for internal
 *research use in your organization is hereby granted, provided that this
 *notice is retained thereon and on all copies. A patent protects the
 *underlying algorithm. This software should not be distributed to anyone
 *outside of your organization without explicit written authorization by the
 *author(s) and MIT. It should not be used for commercial purposes without
 *specific permission from the authors and MIT. MIT also requires written
 *authorization by the author(s) to publish results obtained with the data or
 *software and possibly citation of relevant CBCL reference papers.
 *We make no representation as to the suitability and operability of this
 *data or software for any purpose. It is provided "as is" without express or
 *implied warranty.
 */
#include "vec.h"

using namespace std;

const vec vec::ones(1);
const vec vec::zero(0);

const matrix matrix::eye(1,true);
const matrix matrix::zero(0,false);

// adopted from NRiC, pg 45
matrix matrix::inv(bool &worked) const {
	matrix temp;
	int ix[DIM];
	
	if (!LUdecomp(temp,ix)) {
		worked = false;
		return matrix(0,false);
	}
	worked = true;

	matrix ret;
	vec col;
	for(int j=0;j<DIM;j++) {
		col=0;
		col[j] = 1;
		temp.LUbacksub(ix,col);
		for(int i=0;i<DIM;i++) ret.x[i][j] = col[i];
	}
	return ret;
}

// adopted from NRiC, pg 43
int matrix::LUdecomp(matrix &LU, int *ix) const {

	int d=1;
	LU = *this;
	vec vv;

	for(int i=0;i<DIM;i++) {
		vv[i] = x[i].absmax();
		if (vv[i]==(FP)0.0) return 0;
		vv[i] = 1/vv[i];
	}
	FP sum,big,dum;
	int imax;
	for(int j=0;j<DIM;j++) {
		for(int i=0;i<j;i++) {
			sum = LU.x[i][j];
			for(int k=0;k<i;k++) sum -= LU.x[i][k]*LU.x[k][j];
			LU.x[i][j] = sum;
		}
		big = 0;
		for(int i=j;i<DIM;i++) {
			sum = LU.x[i][j];
			for(int k=0;k<j;k++) sum -= LU.x[i][k]*LU.x[k][j];
			LU.x[i][j] = sum;
			if ((dum=vv[i]*fabs(sum))>=big) {
				big = dum;
				imax = i;
			}
		}
		if (j!=imax) {
			vec dv;
			dv = LU.x[imax];
			LU.x[imax] = LU.x[j];
			LU.x[j] = dv;
			d = -d;
			vv[imax] = vv[j];
		}
		ix[j] = imax;
		if (LU.x[j][j] == 0) LU.x[j][j] = (FP)1.0e-20;
		if (j!=DIM-1) {
			dum = 1/LU.x[j][j];
			for(int i=j+1;i<DIM;i++) LU.x[i][j] *= dum;
		}
	}
	return d;
}

void matrix::LUbacksub(int *ix, vec &b) const {
	int ip,ii=-1;
	FP sum;

	for(int i=0;i<DIM;i++) {
		ip = ix[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii!=-1)
			for(int j=ii;j<=i-1;j++) sum -= x[i](j)*b[j];
		else if (sum!=0) ii=i;
		b[i] = sum;
	}
	for(int i=DIM-1;i>=0;i--) {
		sum = b[i];
		for(int j=i+1;j<=DIM-1;j++) sum -= x[i](j)*b[j];
		b[i] = sum/x[i](i);
	}
}

vec matrix::solve(const vec &b, bool &worked) const {
	vec ret(b);
	int ix[DIM];
	matrix a;
	
	if (!LUdecomp(a,ix)) {
		worked = false;
		return vec(0);
	}
	worked=true;
	a.LUbacksub(ix,ret);
	return ret;
}

