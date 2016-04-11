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
#ifndef _TRIRENDERER_H_
#define _TRIRENDERER_H_
#include <geometry/trirender.h>

void trirenderer(trirender &renderer, bool split,
	const ilist<double> &x, const ilist<double> &y,
	const ilist<double> &z,
	const ilist<double> &r, const ilist<double> &g,
	const ilist<double> &b,
	const ilist<int> &f, float ambient, float direct);

#endif
