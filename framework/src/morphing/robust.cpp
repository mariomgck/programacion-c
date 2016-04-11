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
#include <morphing/robust.h>

void solve(smatrix &A, svec &b, svec &x, svec &w,
		svec &slack, double lambda, double tol, int n) {


	if (lambda < 0) {
		A.solve(b,x,w,tol);
		return;
	}
	bool done;
	svec t(b.s);
	int count;
	double newslack;
	do {
		int i;
		for(i=0;i<b.s;i++)
			t.x[i] = b.x[i] + slack.x[i];
		A.solve(t,x,w,tol);
		A.rightmult(x,t);
		t -= b;
		done = true;
		count = 0;
		for(i=0;i<t.s;i++) {
			if (i>=n) {
				slack.x[i] = 0;
				continue;
			}
			if (t.x[i]<=lambda*w.x[i]) {
				if (t.x[i]>=-lambda*w.x[i]) {
					newslack = 0;
					slack.x[i] = 0;
				} else {
					newslack = t.x[i] + lambda*w.x[i];
					count++;
				}
			} else {
				newslack = t.x[i] - lambda*w.x[i];
				count++;
			}
			if (fabs(newslack - slack.x[i])>0.01) done = false;
			slack.x[i] = newslack;
		}
		cout << "slack variables used: " << count << endl;
	} while (!done);
	cout << "--" << endl;
}
