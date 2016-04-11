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
#include "face.h"
#include "edge.h"
#include "vertex.h"
#include "globaldef.h"
#include <iostream>

using namespace std;
 
FP face::area() const {
	return edges[1]->area(vertices[0]->p);
}

vec face::basepoint() { 
	return vertices[0]->p;
}

vec face::midpoint() const {
	vec ret = 0;

	for(int i=0;i<3;i++)
		ret += vertices[i]->p;
	ret /= 3.0;
	return ret;
}

static inline void swap(FP &a, FP &b) {
	FP c;
	c=a;a=b;b=c;
}

vec face::randompoint() {

	FP ax,ay,bx,by,cx,cy;
	FP chx, chy, maxy,miny;
	vec t1,t2;
	vec tdir;

	if (!tangents(t1,t2)) { // a triangle degenerate to a line or point
		if (t1==0) return vertices[0]->p;
		t1 = vertices[0]->p - vertices[1]->p;
		ax = t1.len2();
		if (ax>0) {
			bx = drand48();
			return vertices[1]->p + bx*t1;
		} 
		t1 = vertices[0]->p - vertices[2]->p;
		ax = t1.len2();
		bx = drand48();
		return vertices[2]->p + bx*t1;
	}
	ax = 0;
	ay = 0;
	tdir = vertices[1]->p - vertices[0]->p;
	//tdir /= tdir.len();
	bx = t1*tdir;
	by = t2*tdir;
	tdir = vertices[2]->p - vertices[0]->p;
	//tdir /= tdir.len();
	cx = t1*tdir;
	cy = t2*tdir;
	// order ax then bx then cx along that axis
	if (cx<ax) {
		swap(ax,cx);
		swap(ay,cy);
	}
	if (bx<ax) {
		swap(ax,bx);
		swap(ay,by);
	}
	if (cx<bx) {
		swap(cx,bx);
		swap(cy,by);
	}

	bool axisflip;
	if (cx-ax < cy-ay) {
		axisflip = true;
		swap(ax,ay);
		swap(bx,by);
		swap(cx,cy);
	} else axisflip = false;
	if (cx-ax == 0.0) return vertices[0]->p;

	// now chose an x (but not uniformly)
	if (drand48() < (bx-ax)/(cx-ax))
		chx = ax + (bx-ax)*sqrt(drand48());
	else chx = bx + (cx-bx)*((FP)1.0-sqrt(drand48()));
	if (chx < bx) // this should be the same as the above test
				// but floating point math is never exact
		miny = ay + (chx - ax)*(by-ay)/(bx-ax);
	else miny = by + (chx - bx)*(cy-by)/(cx-bx);
	maxy = ay + (chx - ax)*(cy-ay)/(cx-ax);
	if (miny>maxy) { swap(miny,maxy); }
	chy = drand48()*(maxy-miny) + miny;

	if (axisflip) {
		swap(chx,chy);
	}
	return vertices[0]->p + chx*t1 + chy*t2;
}

vec face::normal() const {

	vec ret=0;
	FP x1,y1,z1,x2,y2,z2,rx,ry,rz,d;

	x1 = vertices[1]->p[0] - vertices[0]->p[0];
	y1 = vertices[1]->p[1] - vertices[0]->p[1];
	z1 = vertices[1]->p[2] - vertices[0]->p[2];
	x2 = vertices[2]->p[0] - vertices[0]->p[0];
	y2 = vertices[2]->p[1] - vertices[0]->p[1];
	z2 = vertices[2]->p[2] - vertices[0]->p[2];
	rx = y1*z2 - y2*z1;
	ry = z1*x2 - z2*x1;
	rz = x1*y2 - x2*y1;
	d = sqrt(rx*rx+ry*ry+rz*rz);
	if (d==0) {
		ret[0] = 1;
		return ret;
	}
	ret[0] = rx/d;
	ret[1] = ry/d;
	ret[2] = rz/d;
	return ret;
}

bool face::realtangents(vec &t1, vec &t2) {

	FP d;

	tvalid = true;
	vec p1=vertices[1]->p - vertices[0]->p;
	vec p2=vertices[2]->p - vertices[0]->p;
	t1 = p1/p1.len();
	t2 = p2/p2.len();
	if (t1.isvalid() && !t2.isvalid()) {
		t2 = t1;
		t1cache = t1;
		t2cache = t2;
		tretcache = false;
		v0cache = 0;
		u0cache = p1.len()/2;
		r0cache = u0cache*u0cache;
		return false;
	}
	if (!t1.isvalid() && t2.isvalid()) {
		t1 = t2;
		t1cache = t1;
		t2cache = t2;
		tretcache = false;
		v0cache = 0;
		u0cache = p2.len()/2;
		r0cache = u0cache*u0cache;
		return false;
	}
	if (!t1.isvalid() && !t2.isvalid()) {
		t1 = 0;
		t2 = 0;
		t1cache = 0;
		t2cache = 0;
		tretcache = false;
		u0cache = 0;
		v0cache = 0;
		r0cache = 0;
		return false;
		/*
		cout << "invalid tangents for face" << endl;
		exit(1); // this *should* never happen
		*/
	}
	d = t1*t2;
	t2 -= t1*d;
	t2 = t2.norm();
	if (!t2.isvalid()) {
		t2 = t1;
		t1cache = t1;
		t2cache = t2;
		tretcache = false;
		FP mn=0,mx=0,d;
		d = t1*p1; if (d<mn) mn=d; if (d>mx) mx=d;
		d = t1*p2; if (d<mn) mn=d; if (d>mx) mx=d;
		v0cache = 0;
		u0cache = (mx+mn)/2;
		r0cache = u0cache*u0cache;
		return false;
	}
	FP u1,u2,v1,v2,den;
	u1 = p1*t1; u2 = p2*t1; v1 = p1*t2; v2 = p2*t2;
	den = 2*(v2*u1-v1*u2);
	FP s1,s2;
	s1 = v1*v1+u1*u1;
	s2 = v2*v2+u2*u2;
	u0cache = v2*s1-v1*s2/den;
	v0cache = u1*s2-u2*s1/den;
	r0cache = u0cache*u0cache + v0cache*v0cache;
	t1cache = t1;
	t2cache = t2;
	tretcache = true;
	return true;
}

bool face::ptinface(const vec &p) const {

	int s,t,u;
	vec side,su;
	FP l;

	for(s=0;s<3;s++) {
		if (s==2) t=0;
		else t=s+1;
		if (t==2) u=0;
		else u = t+1;
		side = vertices[s]->p - vertices[t]->p;
		su = vertices[u]->p - vertices[t]->p;
		l = side.len();
		if (l==0) return false;
		side /=l;
		su -= (su*side)*side;
		if ((p-vertices[t]->p)*su < 0) return false;
	}
	return true;
}

bool face::ptinface3d(const vec &p) const {

	int s,t,u;
	vec side,su;
	FP l;

	for(s=0;s<3;s++) {
		if (s==2) t=0;
		else t=s+1;
		if (t==2) u=0;
		else u=t+1;
		side = vertices[s]->p - vertices[t]->p;
		su = vertices[u]->p - vertices[t]->p;
		side[3] = side[4] = side[5] = 0;
		su[3] = su[4] = su[5] = 0;
		l = side.len();
		if (l==0) return false;
		side /= l;
		su -= (su*side)*side;
		if ((p-vertices[t]->p)*su < 0) return false;
	}
	return true;
}

bool face::rayintersection(const vec &x0, const vec v, const FP &ceil,
					 FP &dist, vec &ret) {

	vec t1,t2;
	if (!tangents(t1,t2)) return false;
	vec s = x0-vertices[0]->p;
	FP d1 = v*t1;
	FP d2 = v*t2;
	FP dv = v*v;
	FP a = dv-d1*d1-d2*d2;
	if (a<=0) return false; // v is parallel to the plane
	a = ((s*t1)*d1 + (s*t2)*d2 - s*v)/a;
	if (a<0 || a>=ceil) return false;
	vec pt = x0 + a*v;
	if (ptinface(pt)) {
		ret = pt;
		dist = a;
		return true;
	}
	return false;
}

bool face::rayintersection3d(const vec &x0, const vec v, const FP &ceil,
							 FP &dist, vec &ret) const {

	vec n = normal();
	FP a = n*v;

	if (a==0) return false; // ray parallel to plane

	a = (n*(vertices[0]->p - x0))/a;
	if (a<0 || a>=ceil) return false;
	vec pt = x0+a*v;
	if (ptinface3d(pt)) {
		ret = pt;
		dist = a;
		return true;
	}
	return false;
}

vec face::addcolor(const vec &x) {

	vec t1,t2;

	if (!tangents(t1,t2)) return x;
	vec snt1,snt2;
	FP den,dt1,dt2,dt12;

	int i;
	for(i=0;i<3;i++) {
		snt1[i]=t1[i];
		snt2[i]=t2[i];
	}
	for(;i<DIM;i++)
		snt1[i]=snt2[i]=0;

	dt1 = snt1*snt1;
	dt2 = snt2*snt2;
	dt12 = snt1*snt2;
	den = dt1*dt2-dt12*dt12;
	if (den==0) return x;

	vec p;
	FP u,v;
	p = x-vertices[0]->p;
	u = (dt2*(p*snt1) - dt12*(p*snt2))/den;
	v = (dt1*(p*snt2) - dt12*(p*snt1))/den;
	return vertices[0]->p + t1*u + t2*v;
}

int face::closestpoint(const vec &p, vec& ret, const FP &ceil2, FP &dist2, 
				 vertex*& v1, vertex*& v2, const vec &t1, const vec &t2, const FP &anglefact) {

	vec myt1,myt2;
	FP cd; //current best distance
	vec pp; // projected point

	if (tangents(myt1,myt2)) { // otherwise, the triangle has no area
		FP pa;
		if (anglefact>0) {
				pa = anglefact*planeangle(t1,t2,myt1,myt2);
				if (pa>=ceil2) return 0;
		} else pa = 0;
	
		vec dp; // difference vector
		
		pp = p - vertices[0]->p;
		pp = vertices[0]->p + (pp*myt1)*myt1 + (pp*myt2)*myt2;
		dp = pp-p;
		cd = dp.len2() + pa;
		if (cd>=ceil2) return 0;
		FP v,u;
		u = myt1*(p-vertices[0]->p);
		u -= u0cache;
		v = myt2*(p-vertices[0]->p);
		v -= v0cache;
		u = u*u+v*v - r0cache;
		if (cd+u>=ceil2) return 0;
		if (ptinface(pp)) {
			dist2 = cd;
			v1 = v2 = NULL;
			ret = pp;
			return 1;
		}
	} else return 0; // save a few computations -- I think this is okay.
	cd = ceil2;
	
	FP d; // temp
	int j;

	for(int i=0;i<3;i++) {
		j = edges[i]->closestpoint(p,pp);
		d = (p-pp).len2();
		if (d<cd || (d==cd && j)) { // if better match, or an equal match on a corner
			ret = pp;
			cd = d;
			if (j==0) {
				v1 = edges[i]->vertices[0];
				v2 = edges[i]->vertices[1];
			} else if (j==1) {
				v1 = edges[i]->vertices[0];
				v2 = NULL;
			} else {
				v1 = edges[i]->vertices[1];
				v2 = NULL;
			}
		}
	}
	if (cd<ceil2) {
		dist2 = cd;
		return 1;
	}
	return 0;
}

void face::project(const vec &p, FP *alpha) const {
	FP a;
	int i;

	a = area();
	if (a<=0.0) { // there isn't a right answer in this case, but the following
		// will do okay I think (the calling function may be able to determine
		//  the right thing to do but for the moment, we aren't flagging this)
		alpha[0] = alpha[1] = alpha[2] = (FP)1/(FP)3;
		return;
	}
	
	for(int j=0;j<3;j++) {
		if (j==2) i=0;
		else i = j+1;
		alpha[j] = edges[i]->area(p)/a;
	}
}

vec face::unproject(const FP *alpha) const {
	vec ret(0);

	for(int i=0;i<3;i++)
		ret += alpha[i]*vertices[i]->p;
	return ret;

/*
	vec p1,p2,p3,p4;
	int i,j;

	p4 = edges[0]->direction(vertices[0]);
	p3 = vertices[2]->p - vertices[0]->p;
	p3 *= alpha[2];
	p3 += vertices[0]->p;
	for(j=0;j<3;j++) {
		if (j==2) i=0;
		else i=j+1;
		p2 = edges[i]->direction(vertices[i]);
		p1 = vertices[j]->p - vertices[i]->p;
		p1 *= alpha[j];
		p1 += vertices[i]->p;
		ret += lineintersection(p1,p2,p3,p4);
		p3 = p1;
		p4 = p2;
	}
	return ret/3;
*/
}

vec face::unproject(const FP *alpha, const vec &v1, const vec &v2,
			const vec &v3) {

	return alpha[0]*v1 + alpha[1]*v2 + alpha[2]*v3;

/*
	vec ret(0);

	vec p1,p2,p3,p4;
	
	p4 = (v2-v1).norm();
	p3 = v3-v1;
	p3 *= alpha[2];
	p3 += v1;

	p2 = (v3-v2).norm();
	p1 = v1-v2;
	p1 *= alpha[0];
	p1 += v2;
	ret += lineintersection(p1,p2,p3,p4);
	p3 = p1; p4 = p2; // maybe this will be optimized out?

	p2 = (v1-v3).norm();
	p1 = v2-v3;
	p1 *= alpha[1];
	p1 += v3;
	ret += lineintersection(p1,p2,p3,p4);
	p3 = p1; p4 = p2; // ditto?

	p2 = (v2-v1).norm();
	p1 = v3-v1;
	p1 *= alpha[2];
	p1 += v1;
	ret += lineintersection(p1,p2,p3,p4);

	return ret/3;
*/
}
		
void face::updatebounds() {
	for(int i=0;i<DIM;i++)
		mymin[i] = tripmin(vertices[0]->p[i],vertices[1]->p[i],vertices[2]->p[i]);
	for(int i=0;i<DIM;i++)
		mymax[i] = tripmax(vertices[0]->p[i],vertices[1]->p[i],vertices[2]->p[i]);
}

int face::vertexindex(int i) const {
	return vertices[i]->index;
}
