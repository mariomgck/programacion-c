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
#include <warping/warpmesh2.h>
#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>

warpmesh2::warpmesh2(const cyberreader &cr, int subsamplesize, const warpparms &params) :
			awarpmesh(cr,subsamplesize,params) {
	setdata();
}

warpmesh2::warpmesh2(istream &s, const warpparms &params, FP cs, const matrix &A, const vec &b) :
			awarpmesh(s,params,cs,A,b) {
	setdata();
}

warpmesh2::~warpmesh2() {
	delete A;
	int i;
	for(i=0;i<DIM;i++) delete b[i];
}

warpmesh2::warpmesh2(const trimesh &from, const warpparms &params, FP cs,
		const matrix &A, const vec &b) :
			awarpmesh(from,params,cs,A,b) {
	setdata();
}

void warpmesh2::addsprings(ilist<int> *springs) {

	int ttlsp = 0;
	for(int i=0;i<vertices.length();i++) ttlsp += springs[i].length();
	ttlsp /= 2; // we've counted each one twice
	cout << "total number of springs: " << ttlsp << endl;
	A = new smatrix(n+n+ttlsp,vertices.length(),3);
	cout << "matrix created" << endl;
	for(int i=0;i<DIM;i++)
		b[i] = new svec(n+n+ttlsp);
	cout << "vectors created" << endl;
	int in[2];
	FP val[2];
	int j,d,l1,l2;
	int c = n+n;
	vec diff;
	for(int i=0;i<vertices.length();i++) {
		l1 = springs[i].length();
		for(j=0;j<l1;j++,c++) {
			if (springs[i](j)<i) {
				c--;
				continue;
			}
			in[0] = i;
			in[1] = springs[i](j);
			l2 = springs[in[1]].length();
			//cout << l1 << ',' << l2 << endl;
			//val[0] = (l1+l2)/
			// sqrt((origpos[i]-origpos[in[1]]).len())/(l1*l2);
			val[0] = 1/sqrt((origpos[i]-origpos[in[1]]).len());
			val[1] = -val[0];
			A->setrow(c,2,in,val);
			diff = origpos[i]-origpos[in[1]];
			for(d=0;d<DIM;d++) {
				b[d]->x[c] = val[0]*diff[d];
			}
		}
	}
	if (c!=n+n+ttlsp)
		cout << "problem setting up springs" << endl;
}

void warpmesh2::adduserspring(const vec &from, const vec &to) {
	FP alphas[3];
	face *f;
	vec t1,t2,p2;
	vertex *v1,*v2;
	FP d2;
	int in[3];

	p2 = closestpoint(from,d2,f,v1,v2,t1,t2,0);
	f->project(p2,alphas);
	in[0] = f->vertices[0]->index;
	in[1] = f->vertices[1]->index;
	in[2] = f->vertices[2]->index;
	nusp++;
	A->setrow(n+n-nusp,3,in,alphas);
	int i;
	for(i=0;i<DIM;i++)
		b[i]->x[n+n-nusp] = to(i);
}

void warpmesh2::setsystem(trimesh *target) {

	int i,d;
	vec p,p2,t1,t2;
	vertex *v1,*v2;
	face *f,*of;
	FP alphas[3],d2;
	int in[3];

	for(i=0;i<n;i++) {
		p = samplepoint(f,n/4);
		if (rho>0)
			f->tangents(t1,t2);
		p2 = target->closestpoint(p,d2,of,v1,v2,t1,t2,rho);

		f->project(p,alphas);
		for(d=0;d<3;d++)
			in[d] = f->vertices[d]->index;
		A->setrow(i,3,in,alphas);
		for(d=0;d<DIM;d++)
			b[d]->x[i] = p2[d];
	}
	for(i=0;i<n-nusp;i++) {
		p = target->samplepoint(f,n/4);
		if (rho>0)
			f->tangents(t1,t2);
		p2 = closestpoint(p,d2,of,v1,v2,t1,t2,rho);

		of->project(p2,alphas);
		for(d=0;d<3;d++)
			in[d] = of->vertices[d]->index;
		A->setrow(i+n,3,in,alphas);
		for(d=0;d<DIM;d++)
			b[d]->x[i+n] = p[d];
	}
}

void warpmesh2::iterate(trimesh *target) {

	setsystem(target);

	svec *x[DIM];
	svec w(A->getm());
	int d,i;

	for(i=0;i<n+n-nusp;i++)
		w.x[i] = (FP)1/(FP)faces.length();
	for(;i<n+n;i++)
		w.x[i] = zeta;
	for(;i<A->getm();i++)
		w.x[i] = alpha;
	
	for(d=0;d<DIM;d++) {
		x[d] = new svec(vertices.length());
		for(i=0;i<vertices.length();i++) {
			x[d]->x[i] = vertices[i]->p[d];
		}
		A->solve(*(b[d]),*(x[d]),w,
			0.1);
		//	1e-6/(FP)(faces.length()*faces.length()));
	}
	vec pos;
	for(i=0;i<vertices.length();i++) {
		for(d=0;d<DIM;d++)
			pos[d] = x[d]->x[i];
		movevertex(i,pos);
	}
	for(d=0;d<DIM;d++)
		delete x[d];

	alpha *= eta;
	zeta *= xi;
}
