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
/* this is a poorly coded C (not C++) implementation of a triangle
 * depth sorter and splitter -- it passes the triangles off to a renderer
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <containers/ilist.h>
#include <geometry/trirenderer.h>


#define P 1e-10

#define E(a,b) ((b)-(a) < P && (b)-(a) > -P)

class renface {
public:
	double ox1,oy1,oz1,ox2,oy2,oz2,ox3,oy3,oz3;
	double r1,r2,r3,g1,g2,g3,b1,b2,b3;
	double nx,ny,nz;
	double d;
	int p1,p2,p3; // indices into list of renpt's
	int e1,e2,e3;
	int mark;
	ilist<renface *> splitby;
};

class renpt {
public:
	double x, y, z;
};

int hasbeensplitby(renface *f1, renface *f2) {
	for(int i=0;i<f1->splitby.length();i++) {
		if (f1->splitby[i]->nx == f2->nx &&
		    f1->splitby[i]->ny == f2->ny &&
		    f1->splitby[i]->nz == f2->nz &&
		    f1->splitby[i]->d == f2->d) return 1;
	}
	return 0;
}

int sameplane(renface *f1, renface *f2) {
	return (f1->nx == f2->nx && f1->ny == f2->ny && f1->nz == f2->nz &&
		f1->d == f2->d);
}

renface *createnewface(ilist<renpt> &ptlist, int i1, int i2, int i3,
	double r1, double g1, double b1, double r2, double g2, double b2,
	double r3, double g3, double b3) {

	renface *f;

	f = new renface;
	f->ox1 = ptlist[i1].x;
	f->oy1 = ptlist[i1].y;
	f->oz1 = ptlist[i1].z;
	f->ox2 = ptlist[i2].x;
	f->oy2 = ptlist[i2].y;
	f->oz2 = ptlist[i2].z;
	f->ox3 = ptlist[i3].x;
	f->oy3 = ptlist[i3].y;
	f->oz3 = ptlist[i3].z;
/*
	if ((f->ox1==f->ox2 && f->oy1==f->oy2) ||
	    (f->ox1==f->ox3 && f->oy1==f->oy3) ||
	    (f->ox2==f->ox3 && f->oy2==f->oy3)) {
		delete f;
		return NULL;
	}
*/
	f->r1 = r1; f->r2 = r2; f->r3 = r3;
	f->g1 = g1; f->g2 = g2; f->g3 = g3;
	f->b1 = b1; f->b2 = b2; f->b3 = b3;
	{
		double x1,y1,z1,x2,y2,z2,rx,ry,rz,d;
		x1 = f->ox2-f->ox1;
		y1 = f->oy2-f->oy1;
		z1 = f->oz2-f->oz1;
		x2 = f->ox3-f->ox1;
		y2 = f->oy3-f->oy1;
		z2 = f->oz3-f->oz1;
		rx = y1*z2 - y2*z1;
		ry = z1*x2 - z2*x1;
		rz = x1*y2 - x2*y1;
		d = sqrt(rx*rx+ry*ry+rz*rz);
		if (E(d,0)) {
			delete f;
			return NULL;
		}
		if (rz<0) {
			f->nx = -rx/d;
			f->ny = -ry/d;
			f->nz = -rz/d;
		} else {
			f->nx = rx/d;
			f->ny = ry/d;
			f->nz = rz/d;
		}
		f->d = f->nx*f->ox1 + f->ny*f->oy1 + f->nz*f->oz1;
	}
	f->p1 = i1;
	f->p2 = i2;
	f->p3 = i3;
	f->mark = 0;
	f->e1 = f->e2 = f->e3 = 1;
	return f;
}
	
void maketriangles(const ilist<double> &x, const ilist<double> &y,
	const ilist<double> &z, const ilist<double> &r,
	const ilist<double> &g, const ilist<double> &b,
	const ilist<int> &fa, ilist<renpt> &ptlist, ilist<renface *> &flist,
	double &minx, double &maxx, double &miny, double &maxy,
	bool split) {

	int i;
	renpt p;
	renface *f;

	maxx=maxy = -HUGE_VAL; minx=miny = HUGE_VAL;
	for(i=0;i<x.length();i++) {
		p.x = x(i);
		p.y = y(i);
		p.z = z(i);
		if (p.x < minx) minx = p.x;
		if (p.y < miny) miny = p.y;
		if (p.x > maxx) maxx = p.x;
		if (p.y > maxy) maxy = p.y;
		ptlist += p;
	}
	for(i=0;i<fa.length();i+=3) {
		f = createnewface(ptlist,fa(i),fa(i+1),fa(i+2),
			r(fa(i)),g(fa(i)),b(fa(i)),
			r(fa(i+1)),g(fa(i+1)),b(fa(i+1)),
			r(fa(i+2)),g(fa(i+2)),b(fa(i+2)));
		if (f) flist+=f;
	}
}

int linesintersect(double x1, double y1, double x2, double y2,
	double xx1, double yy1, double xx2, double yy2) {

	double den,v;

	x2 -= x1; y2 -= y1;
	xx2 -= xx1; yy2 -= yy1;
	
	den = y2*xx2 - x2*yy2;
	if (E(den,0)) return 0; // parallel -- probably not intersecting ;)
	v = ((y1-yy1)*x2 - (x1-xx1)*y2) / -den;
	if (E(v,0) || E(v,1) || v<0 || v>1) return 0;
	v = ((yy1-y1)*xx2 - (xx1-x1)*yy2) / den;
	return (!E(v,0) && !E(v,1) && v>=0 && v<=1);
}

int ptin(double x1, double y1, double x2, double y2, double x3, 
	double y3, double xx, double yy) {

	double cr1,cr2;
	double dx1,dy1,dx2,dy2,dx3,dy3;

	dx1 = x1-x2; dy1 = y1-y2; dx2 = x3-x2; dy2 = y3-y2;
	dx3 = xx-x2; dy3 = yy-y2;
	cr1 = dx1*dy2 - dx2*dy1;
	cr2 = dx1*dy3 - dx3*dy1;
	if ((E(cr1,0) || cr1>0) && (E(cr2,0) || cr2<0)) return 0;
	if ((E(cr1,0) || cr1<0) && (E(cr2,0) || cr2>0)) return 0;

	dx1 = x2-x3; dy1 = y2-y3; dx2 = x1-x3; dy2 = y1-y3;
	dx3 = xx-x3; dy3 = yy-y3;
	cr1 = dx1*dy2 - dx2*dy1;
	cr2 = dx1*dy3 - dx3*dy1;
	if ((E(cr1,0) || cr1>0) && (E(cr2,0) || cr2<0)) return 0;
	if ((E(cr1,0) || cr1<0) && (E(cr2,0) || cr2>0)) return 0;
	
	dx1 = x3-x1; dy1 = y3-y1; dx2 = x2-x1; dy2 = y2-y1;
	dx3 = xx-x1; dy3 = yy-y1;
	cr1 = dx1*dy2 - dx2*dy1;
	cr2 = dx1*dy3 - dx3*dy1;
	if ((E(cr1,0) || cr1>0) && (E(cr2,0) || cr2<0)) return 0;
	if ((E(cr1,0) || cr1<0) && (E(cr2,0) || cr2>0)) return 0;

	return 1;
}

int realoverlap(renface *f1, renface *f2,ilist<renpt> &ptlist) {

	double x11,y11,x12,y12,x13,y13;
	double x21,y21,x22,y22,x23,y23;

	x11 = ptlist[f1->p1].x;
	y11 = ptlist[f1->p1].y;
	x12 = ptlist[f1->p2].x;
	y12 = ptlist[f1->p2].y;
	x13 = ptlist[f1->p3].x;
	y13 = ptlist[f1->p3].y;
	x21 = ptlist[f2->p1].x;
	y21 = ptlist[f2->p1].y;
	x22 = ptlist[f2->p2].x;
	y22 = ptlist[f2->p2].y;
	x23 = ptlist[f2->p3].x;
	y23 = ptlist[f2->p3].y;
	return (ptin(x11,y11,x12,y12,x13,y13,x21,y21) ||
	        ptin(x11,y11,x12,y12,x13,y13,x22,y22) ||
	        ptin(x11,y11,x12,y12,x13,y13,x23,y23) ||
	        ptin(x21,y21,x22,y22,x23,y23,x11,y11) ||
	        ptin(x21,y21,x22,y22,x23,y23,x13,y12) ||
	        ptin(x21,y21,x22,y22,x23,y23,x13,y13) ||
	        linesintersect(x11,y11,x12,y12,x21,y21,x22,y22) ||
	        linesintersect(x11,y11,x12,y12,x22,y22,x23,y23) ||
	        linesintersect(x11,y11,x12,y12,x21,y21,x23,y23) ||
	        linesintersect(x12,y12,x13,y13,x21,y21,x22,y22) ||
	        linesintersect(x12,y12,x13,y13,x22,y22,x23,y23) ||
	        linesintersect(x12,y12,x13,y13,x21,y21,x23,y23) ||
	        linesintersect(x11,y11,x13,y13,x21,y21,x22,y22) ||
	        linesintersect(x11,y11,x13,y13,x22,y22,x23,y23) ||
	        linesintersect(x11,y11,x13,y13,x21,y21,x23,y23)); 
}
	
int bboverlap(renface *f1, renface *f2,ilist<renpt> &ptlist) {

	double maxx,minx,maxy,miny,maxz,minz;
	double a,b,c;

	maxx = ptlist[f2->p1].x;
	minx = ptlist[f2->p1].x;
	if (ptlist[f2->p2].x < maxx) 
		minx = ptlist[f2->p2].x;
	else maxx = ptlist[f2->p2].x;
	if (ptlist[f2->p3].x > maxx)
		maxx = ptlist[f2->p3].x;
	if (ptlist[f2->p3].x < minx)
		minx = ptlist[f2->p3].x;
	maxy = ptlist[f2->p1].y;
	miny = ptlist[f2->p1].y;
	if (ptlist[f2->p2].y < maxy) 
		miny = ptlist[f2->p2].y;
	else maxy = ptlist[f2->p2].y;
	if (ptlist[f2->p3].y > maxy)
		maxy = ptlist[f2->p3].y;
	if (ptlist[f2->p3].y < miny)
		miny = ptlist[f2->p3].y;
	maxz = ptlist[f2->p1].z;
	minz = ptlist[f2->p1].z;
	if (ptlist[f2->p2].z < maxz) 
		minz = ptlist[f2->p2].z;
	else maxz = ptlist[f2->p2].z;
	if (ptlist[f2->p3].z > maxz)
		maxz = ptlist[f2->p3].z;
	if (ptlist[f2->p3].z < minz)
		minz = ptlist[f2->p3].z;
	a = ptlist[f1->p1].x;
	b = ptlist[f1->p2].x;
	c = ptlist[f1->p3].x;
	if (a <= minx && b <= minx && c <= minx) return 0;
	if (a >= maxx && b >= maxx && c >= maxx) return 0;
	a = ptlist[f1->p1].y;
	b = ptlist[f1->p2].y;
	c = ptlist[f1->p3].y;
	if (a <= miny && b <= miny && c <= miny) return 0;
	if (a >= maxy && b >= maxy && c >= maxy) return 0;
/*
	a = ptlist[f1->p1].z;
	b = ptlist[f1->p2].z;
	c = ptlist[f1->p3].z;
	if (a <= minz && b <= minz && c <= minz)
		return 0;
	if (a >= maxz && b >= maxz && c >= maxz)
		return 0;
*/
	return 1;
}

int breakedge(double x1, double y1, double z1, double x2, double y2,
	double z2, double A, double B, double C, double d, double &x,
	double &y, double &z) {

	double t;
	double dx,dy,dz;

	dx = x2-x1; dy = y2-y1; dz = z2-z1;
	t = A*dx + B*dy + C*dz;
	if (E(t,0)) return 0;
	t = (d - A*x1 - B*y1 - C*z1) / t;
	if (E(t,0)) {
		x = x1; y = y1; z = z1; return 2;
	}
	if (E(t,1)) {
		x = x2; y = y2; z = z2; return 2;
	}
	if (t<0 || t>1) return 0;
	x = x1+t*dx;
	y = y1+t*dy;
	z = z1+t*dz;
	return 1;
}

double planedot(renface *f1, double x, double y, double z) {
	return x*f1->nx + y*f1->ny + z*f1->nz - f1->d;
}

double facezval(renface *f1,ilist<renpt> &p) {

	double min;

	if (p[f1->p1].z<p[f1->p2].z)
		min = p[f1->p1].z;
	else min = p[f1->p2].z;
	if (min<p[f1->p3].z) return min;
	return p[f1->p3].z;
}

ilist<renpt> *ptl;

int fcomp(const void *f1, const void *f2) {

	double d;

	d = facezval(*(renface **)f1,*ptl)-facezval(*(renface **)f2,*ptl);
	if (d<0) return -1;
	if (d>0) return 1;
	return 0;
}

void sorttriangles(ilist<renface *> &flist, ilist<renpt> &ptlist) {

/*
	renface *temp;
	int i,b,j;

	for(i=0;i<flist.length()-1;i++) {
		b = i;
		for(j=i+1;j<flist.length();j++) {
			if (facezval(flist[j],ptlist)<facezval(flist[b],ptlist))
				b = j;
		}
		if (b!=i) {
			temp = flist[i];
			flist[i] = flist[b];
			flist[b] = temp;
		}
	}
*/
	ptl = &ptlist;
	flist.sort(fcomp);
}

void insertflist(ilist<renface *> &flist, ilist<renpt> &ptlist, renface *f, int
	start) {

	int i;
	double v;
	renface *t,*s;

	v = facezval(f,ptlist);
	for(i=start;i<flist.length();i++)
		if (facezval(flist[i],ptlist)>v) break;
	if (i==flist.length()) {
		flist += f;
		return;
	}
	t = flist[i];
	flist[i] = f;
	for(i++;i<flist.length();i++) {
		s = flist[i];
		flist[i] = t;
		t = s;
	}
	flist += t;
}
		
int splitby(ilist<renface *> &flist, int a, int b, ilist<renpt> &ptlist,int s) {

	renface *f1,*f2;
	renpt pt1,pt2;
	double x,y,z,test,origz;
	int c,free,i,j,t,tt;
	renface *nf1, *nf2;

	f1 = flist[a]; f2 = flist[b];
	origz = facezval(f2,ptlist);
	t = 0;
	if (f1->nx==f2->nx && f1->ny==f2->ny && f1->nz==f1->nz) return 0;
	if (tt=breakedge(ptlist[f2->p1].x,ptlist[f2->p1].y,ptlist[f2->p1].z,
	              ptlist[f2->p2].x,ptlist[f2->p2].y,ptlist[f2->p2].z,
	              f1->nx,f1->ny,f1->nz,f1->d,x,y,z)) {
		c=1;
		t += tt;
		pt1.x = x; pt1.y = y; pt1.z = z;
	} else { c=0; free=1; }
	if (tt=breakedge(ptlist[f2->p2].x,ptlist[f2->p2].y,ptlist[f2->p2].z,
	              ptlist[f2->p3].x,ptlist[f2->p3].y,ptlist[f2->p3].z,
	              f1->nx,f1->ny,f1->nz,f1->d,x,y,z)) {
		c++;
		t += tt;
		if (c==1) { pt1.x = x; pt1.y = y; pt1.z = z; }
		else { pt2.x = x; pt2.y = y; pt2.z = z; }
	} else { free=2; if (c==0) return 0; }
	if (tt=breakedge(ptlist[f2->p1].x,ptlist[f2->p1].y,ptlist[f2->p1].z,
	              ptlist[f2->p3].x,ptlist[f2->p3].y,ptlist[f2->p3].z,
	              f1->nx,f1->ny,f1->nz,f1->d,x,y,z)) {
		c++;
		t += tt;
		if (c==1) { pt1.x = x; pt1.y = y; pt1.z = z; }
		else { pt2.x = x; pt2.y = y; pt2.z = z; }
	} else { free=3; if (c<2) return 0; }
	if (t<2 || t>3) return 0;
	if (E(pt1.x,pt2.x) && E(pt1.y,pt2.y) && E(pt1.z,pt2.z)) return 0;
	if (s) {
		if (s==1 && free==2) return 0;
		if (s==2 && free==3) return 0;
		if (s==3 && free==1) return 0;
	}
	nf1 = new renface;
	nf2 = new renface;
	for(i=0;i<f2->splitby.length();i++) {
		nf1->splitby += f2->splitby[i];
		nf2->splitby += f2->splitby[i];
	}
	f2->splitby += f1;
	nf1->splitby += f1;
	nf2->splitby += f1;
	nf1->ox1 = f2->ox1; nf1->oy1 = f2->oy1; nf1->oz1 = f2->oz1;
	nf1->ox2 = f2->ox2; nf1->oy2 = f2->oy2; nf1->oz2 = f2->oz2;
	nf1->ox3 = f2->ox3; nf1->oy3 = f2->oy3; nf1->oz3 = f2->oz3;
	nf2->ox1 = f2->ox1; nf2->oy1 = f2->oy1; nf2->oz1 = f2->oz1;
	nf2->ox2 = f2->ox2; nf2->oy2 = f2->oy2; nf2->oz2 = f2->oz2;
	nf2->ox3 = f2->ox3; nf2->oy3 = f2->oy3; nf2->oz3 = f2->oz3;
	nf1->r1 = f2->r1; nf1->r2 = f2->r2; nf1->r3 = f2->r3;
	nf1->g1 = f2->g1; nf1->g2 = f2->g2; nf1->g3 = f2->g3;
	nf1->b1 = f2->b1; nf1->b2 = f2->b2; nf1->b3 = f2->b3;
	nf2->r1 = f2->r1; nf2->r2 = f2->r2; nf2->r3 = f2->r3;
	nf2->g1 = f2->g1; nf2->g2 = f2->g2; nf2->g3 = f2->g3;
	nf2->b1 = f2->b1; nf2->b2 = f2->b2; nf2->b3 = f2->b3;
	nf1->e1 = nf2->e1 = f2->e1;
	nf1->e2 = nf2->e2 = f2->e2;
	nf1->e3 = nf2->e3 = f2->e3;
	nf1->nx = nf2->nx = f2->nx;
	nf1->ny = nf2->ny = f2->ny;
	nf1->nz = nf2->nz = f2->nz;
	nf1->d = nf2->d = f2->d;
	nf1->mark = nf2->mark = 0;
	i = ptlist.length();
	j = i+1;
	ptlist += pt1;
	ptlist += pt2;
	if (free==1) {
		nf1->p1 = f2->p1;
		nf1->p2 = i;
		nf1->p3 = j;
		nf2->p1 = f2->p1;
		nf2->p2 = f2->p2;
		nf2->p3 = i;
		f2->p1 = j;
		f2->p2 = i;
		f2->e1 = 0;
		nf1->e1 = 0;
		nf1->e2 = 0;
		nf2->e3 = 0;
	} else if (free==2) {
		nf1->p1 = j;
		nf1->p2 = i;
		nf1->p3 = f2->p3;
		nf2->p1 = i;
		nf2->p2 = f2->p2;
		nf2->p3 = f2->p3;
		f2->p2 = i;
		f2->p3 = j;
		f2->e2 = 0;
		nf1->e1 = 0;
		nf1->e2 = 0;
		nf2->e3 = 0;
	} else {
		nf1->p1 = i;
		nf1->p2 = j;
		nf1->p3 = f2->p3;
		nf2->p1 = f2->p1;
		nf2->p2 = i;
		nf2->p3 = f2->p3;
		f2->p1 = i;
		f2->p3 = j;
		f2->e3 = 0;
		nf1->e1 = 0;
		nf1->e3 = 0;
		nf2->e2 = 0;
	}
	if (facezval(f2,ptlist)==origz) {
		insertflist(flist,ptlist,nf1,a);
		insertflist(flist,ptlist,nf2,a);
		return 1;
	} else {
		if (facezval(nf1,ptlist)==origz) {
			flist[b] = nf1;
			insertflist(flist,ptlist,f2,a);
			insertflist(flist,ptlist,nf2,a);
		} else {
			flist[b] = nf2;
			insertflist(flist,ptlist,f2,a);
			insertflist(flist,ptlist,nf1,a);
		}
	}
	return 2;
}

int share2pts(renface *f1, renface *f2, int &s) {
	int c;
	
	c=0; s=0;
	if (f1->p1==f2->p1 || f1->p1==f2->p2 || f1->p1==f2->p3) {
		c++;s=1; }
	if (f1->p2==f2->p1 || f1->p2==f2->p2 || f1->p2==f2->p3) {
		c++;s=2; }
	if (f1->p3==f2->p1 || f1->p3==f2->p2 || f1->p3==f2->p3) {
		c++;s=3; }
	if (c!=1) s=0;
	return c>1;
}

int testordering(renface *f1, renface *f2, ilist<renpt> &ptlist) {

	int i;
	double a,b,c;

	if (sameplane(f1,f2)) return 1;
	a = planedot(f1,ptlist[f2->p1].x,ptlist[f2->p1].y,ptlist[f2->p1].z);
	b = planedot(f1,ptlist[f2->p2].x,ptlist[f2->p2].y,ptlist[f2->p2].z);
	c = planedot(f1,ptlist[f2->p3].x,ptlist[f2->p3].y,ptlist[f2->p3].z);
	if (hasbeensplitby(f2,f1)) {
		if (a>P || b>P || c>P) return 1;
		if (a<P || b<P || c<P) return 2;
	}
	if ((a>=0||E(a,0)) && (b>=0||E(b,0)) && (c>=0||E(c,0)) &&
	    (!E(a,0) || !E(b,0) || !E(c,0))) {
		return 1;
	}
	if ((a<=0||E(a,0)) && (b<=0||E(b,0)) && (c<=0||E(c,0)) &&
	    (!E(a,0) || !E(b,0) || !E(c,0))) {
		return 2;
	}
	a = planedot(f2,ptlist[f1->p1].x,ptlist[f1->p1].y,ptlist[f1->p1].z);
	b = planedot(f2,ptlist[f1->p2].x,ptlist[f1->p2].y,ptlist[f1->p2].z);
	c = planedot(f2,ptlist[f1->p3].x,ptlist[f1->p3].y,ptlist[f1->p3].z);
	if (hasbeensplitby(f1,f2)) {
		if (a>P || b>P || c>P) return 2;
		if (a<P || b<P || c<P) return 1;
	}
	if ((a>=0||E(a,0)) && (b>=0||E(b,0)) && (c>=0||E(c,0)) &&
	    (!E(a,0) || !E(b,0) || !E(c,0))) {
		return 2;
	}
	if ((a<=0||E(a,0)) && (b<=0||E(b,0)) && (c<=0||E(c,0)) &&
	    (!E(a,0) || !E(b,0) || !E(c,0))) {
		return 1;
	}
	return 0;
}

int isinfront(renface *f1, renface *f2,ilist<renpt> &p) {

	double max,abs;
	int ret;

	abs = planedot(f2,p[f1->p1].x,p[f1->p1].y,p[f1->p1].z);
	if (abs<0) {
		max = -abs;
		ret = 0;
	} else {
		max = abs;
		ret = 1;
	}
	abs = planedot(f2,p[f1->p2].x,p[f1->p2].y,p[f1->p2].z);
	if (abs<0) {
		if (-abs > max) {
			max = -abs;
			ret = 0;
		}
	} else {
		if (abs > max) {
			max = abs;
			ret = 1;
		}
	}
	abs = planedot(f2,p[f1->p3].x,p[f1->p3].y,p[f1->p3].z);
	if (abs<0) {
		if (-abs > max) {
			max = -abs;
			ret = 0;
		}
	} else {
		if (abs > max) {
			max = abs;
			ret = 1;
		}
	}
	return ret;
}

#define RESTRICT(v) if (v<0) v=0; else if (v>1) v = 1;

void trirenderer(trirender &rend, bool split,
	const ilist<double> &xs, const ilist<double> &ys,
	const ilist<double> &zs,
	const ilist<double> &rs, const ilist<double> &gs,
	const ilist<double> &bs,
	const ilist<int> &fs, float ambient, float direct) {


	ilist<renpt> ptlist;
	ilist<renface *> flist;
	renface *f;
	double m,currz,nm;
	int h,i,j,k,c,tot,s;
	double maxx,minx,maxy,miny;
	ilist<renface *> unmark;

	double x[3],y[3],r[3],g[3],b[3],cx[3],cy[3];
	bool e[3];

	maketriangles(xs,ys,zs,rs,gs,bs,fs,ptlist,flist,
		minx,maxx,miny,maxy,split);
	rend.start(minx,maxx,miny,maxy);
	sorttriangles(flist,ptlist);
	for(k=0;k<flist.length();flist[k++]->mark=0);
	j = 1;
	for(i=0;i<flist.length();i++) {
		f = flist[i];
		h = 1;
		if (split) {
			if (ptlist[f->p1].z<ptlist[f->p2].z)
				m = ptlist[f->p2].z;
			else m= ptlist[f->p1].z;
			if (m<ptlist[f->p3].z) m = ptlist[f->p3].z;
			for(;j<flist.length()&&
			     facezval(flist[j],ptlist)<m;j++) {
				if (!bboverlap(f,flist[j],ptlist)) continue;
				if (!realoverlap(f,flist[j],ptlist)) continue;
				h = testordering(f,flist[j],ptlist);
				if (h!=1) break;
			}
		}
		if (h==1) {
			for(k=0;k<unmark.length();unmark[k++]->mark=0);
			unmark.setlength(0);
			m = f->nz;
			cx[0] = f->ox1; cx[1] = f->ox2; cx[2] = f->ox3;
			cy[0] = f->oy1; cy[1] = f->oy2; cy[2] = f->oy3;
			r[0] = f->r1*(m*direct+ambient);
			RESTRICT(r[0]);
			r[1] = f->r2*(m*direct+ambient);
			RESTRICT(r[1]);
			r[2] = f->r3*(m*direct+ambient);
			RESTRICT(r[2]);
			g[0] = f->g1*(m*direct+ambient);
			RESTRICT(g[0]);
			g[1] = f->g2*(m*direct+ambient);
			RESTRICT(g[1]);
			g[2] = f->g3*(m*direct+ambient);
			RESTRICT(g[2]);
			b[0] = f->b1*(m*direct+ambient);
			RESTRICT(b[0]);
			b[1] = f->b2*(m*direct+ambient);
			RESTRICT(b[1]);
			b[2] = f->b3*(m*direct+ambient);
			RESTRICT(b[2]);
			x[0] = ptlist[f->p1].x; y[0] = ptlist[f->p1].y;
			x[1] = ptlist[f->p2].x; y[1] = ptlist[f->p2].y;
			x[2] = ptlist[f->p3].x; y[2] = ptlist[f->p3].y;
			e[0] = f->e1; e[1] = f->e2; e[2] = f->e3;
			rend.rendertriangle(x,y,e,cx,cy,r,g,b);
			j = i+2;
		} else {
			if (h==0 || flist[j]->mark) {
				h = 0;
/*
				for(k=0;k<flist[j]->splitby.length();k++)
					if (flist[j]->splitby[k]->nx==f->nx &&
					    flist[j]->splitby[k]->ny==f->ny &&
					    flist[j]->splitby[k]->nz==f->nz &&
					    flist[j]->splitby[k]->d==f->d)
						break;
				s = 0;
				if (flist[j]->splitby.length()!=k ||
				    share2pts(f,flist[j],s)) {
					if (flist[j]->mark ||
					    isinfront(flist[j],f,ptlist)) {
						j++;
						i--;
					} else h=2;
				} else {
*/
					if (!splitby(flist,i,j,ptlist,s)) {
						j++;
					} else {
						for(k=0;k<unmark.length();
						          unmark[k++]->mark=0);
						unmark.setlength(0);
					}
					i--;
/*
				}
*/
			}
			if (h==2) {
				flist[i] = flist[j];
				flist[j] = f;
				unmark += f;
				f->mark = 1;
				i--;
				j = i+2;
			}
		}
	}

	for(i=0;i<flist.length();i++)
		delete flist[i];

	rend.end();
}
