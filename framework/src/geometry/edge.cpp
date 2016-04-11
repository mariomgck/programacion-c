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
#include "edge.h"
#include "face.h"
#include "vertex.h"

FP edge::creaseangle(face *ff, face *rf, vertex *rv1, vertex *rv2) const {

	if (faces.length()<2) return -1;

	vec t,ov1,ov2;
	face *f1,*f2;
	vertex *vv;
	
	if (faces.length()<2) return -1;
	t = (vertices[0]->p - vertices[1]->p).norm();
	f1 = faces(0);
	f2 = faces(1);
	if (f1==ff) f1=rf;
	if (f2==ff) f2=rf;
	int i;
	for(i=0;i<3;i++) {
		vv = f1->vertices[i];
		if (vv!=vertices[0] && vv!=vertices[1] && vv!=rv1 && vv!=rv2) {
			ov1 = vv->p - vertices[1]->p;
			break;
		}
	}
	if (i==3) return -1;
	for(i=0;i<3;i++) {
		vv = f2->vertices[i];
		if (vv!=vertices[0] && vv!=vertices[1] && vv!=rv1 && vv!=rv2) {
			ov2 = vv->p - vertices[1]->p;
			break;
		}
	}
	if (i==3) return -1;
	FP d = -((ov1-(ov1*t)*t).norm() * (ov2-(ov2*t)*t).norm());
	return d;
}

FP edge::area(const vec &p) const {

	FP l1,l2,l3,s;

	l1 = len();
	l2 = (vertices[0]->p-p).len();
	l3 = (vertices[1]->p-p).len();
	s = (l1+l2+l3)/2;
	if (l1>=s || l2>=s || l3>=s) return 0;
	return sqrt(s*(s-l1)*(s-l2)*(s-l3));
}

int edge::closestpoint(const vec &p, vec &ret) const {

	vec tp,pp;
	FP len,l2;

	tp = vertices[1]->p - vertices[0]->p;
	len = tp.len();
	tp = tp/len;
	pp = p-vertices[0]->p;
	l2 = tp*pp;
	if (l2 < len && l2 > 0) {
		pp = tp*l2;
		ret = vertices[0]->p + pp;
		return 0;
	}
	if (l2 < len/2) {
		ret = vertices[0]->p;
		return 1;
	}
	ret = vertices[1]->p;
	return 2;
}
