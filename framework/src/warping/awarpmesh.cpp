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
#include <warping/awarpmesh.h>
#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>
#include <iostream>

using namespace std;

awarpmesh::awarpmesh(const cyberreader &cr, int subsamplesize, const warpparms &params) :
			trimesh(cr,subsamplesize) {

	alpha = params.alpha0;
	eta = params.eta;
	xi = params.xi;
	momentum = params.momentum;
	maxmove2 = params.maxmove2;
	rho = params.rho;
	epsilon = params.epsilon;
	n = params.n*numvertices();
	maxsp = params.maxsp;
	zeta = params.zeta;
	//setdata();
}

awarpmesh::awarpmesh(istream &s, const warpparms &params, FP cs, const matrix &A, const vec &b) :
			trimesh(s,cs,true,A,b) {
	alpha = params.alpha0;
	eta = params.eta;
	xi = params.xi;
	momentum = params.momentum;
	maxmove2 = params.maxmove2;
	rho = params.rho;
	epsilon = params.epsilon;
	n = params.n*numvertices();
	maxsp = params.maxsp;
	zeta = params.zeta;
	//setdata();
}

awarpmesh::awarpmesh(const trimesh &from, const warpparms &params, FP cs,
		const matrix &A, const vec &b) :
	trimesh(from,cs,true,A,b) {
	
	alpha = params.alpha0;
	eta = params.eta;
	xi = params.xi;
	momentum = params.momentum;
	maxmove2 = params.maxmove2;
	rho = params.rho;
	epsilon = params.epsilon;
	n = params.n*numvertices();
	maxsp = params.maxsp;
	zeta = params.zeta;
	//setdata();
}

awarpmesh::~awarpmesh() {
	delete []origpos;
}

class findspr : public geoiterate<vertex *> {
public:
	FP lim2;
	int maxnum;
	ilist<FP> *dists;
	ilist<int> *pts;
	vec p;

	virtual FP process(vertex *v) {
		FP d2;

		if (v->position() == p)
			return dists->length()==maxnum ? dists->nth(maxnum-1) : lim2;
		d2 = (v->position() - p).len2();
		if (d2>lim2) 
			return dists->length()==maxnum ? dists->nth(maxnum-1) : lim2;
		if (pts->length()==maxnum && dists->nth(maxnum-1) <= d2)
			return dists->nth(maxnum-1);
		if (pts->length()<maxnum) {
			(*dists) += d2;
			(*pts) += v->getindex();
		}
		int i;
		for(i=pts->length()-2;i>=0;i--) {
			if (dists->nth(i) > d2) {
				dists->setnth(dists->nth(i),i+1);
				pts->setnth(pts->nth(i),i+1);
			} else break;
		}
		dists->setnth(d2,i+1);
		pts->setnth(v->getindex(),i+1);
		return dists->length()==maxnum ? dists->nth(maxnum-1) : lim2;
	}

	virtual FP initaldist() {
		return dists->length()==maxnum ? dists->nth(maxnum-1) : lim2;
	}
};

void awarpmesh::setdata() {
	FP lim2,t;
	vertex *v;
	ilist<FP> td;
	findspr findclose;
	ilist<int> *springs = new ilist<int>[vertices.length()];

	origpos = new vec[vertices.length()];
	for(int i=0;i<vertices.length();i++) {
		v = vertices[i];
		origpos[i] = v->p;
		if (maxsp>0) {
			lim2 = 0;
			for(int j=0;j<v->edges.length();j++)
				if ((t=v->edges[j]->len2())>lim2)
					lim2 = t;
			findclose.maxnum = maxsp;
			findclose.lim2 = lim2*4;
			findclose.pts = springs+i;
			findclose.dists = &td;
			findclose.p = v->p;
			vertexgrid.search(findclose.p,findclose);
			td.setlength(0);
		}
	}
	// add any edges that were omitted due to the maxnum
	//  constraint
	int vi;
	for(int i=0;i<vertices.length();i++) {
		v = vertices[i];
		for(int j=0;j<v->edges.length();j++) {
			if (v->edges[j]->vertices[0] == v)
				vi = v->edges[j]->vertices[1]->index;
			else vi = v->edges[j]->vertices[0]->index;
			if (springs[i].find(vi)==-1)
				springs[i] += vi;
		}
	}
	// make relation reflexive
	for(int i=0;i<vertices.length();i++)
		for(int j=0;j<springs[i].length();j++)
			if (springs[springs[i](j)].find(i) == -1)
				springs[springs[i](j)] += i;
	addsprings(springs);
	delete []springs;
}

vec *awarpmesh::D(bool colorcorrect) const {
	vec *ret;

	ret = new vec[vertices.length()];
	for(int i=0;i<vertices.length();i++) {
		ret[i] = vertices(i)->p - origpos[i];
		if (colorcorrect) {
			for(int j=3;j<DIM;j++)
				ret[i][j] /= colorscale;
		}
	}
	return ret;
}

vec *awarpmesh::P(bool colorcorrect) const {
	vec *ret = new vec[vertices.length()];
	for(int i=0;i<vertices.length();i++) {
		ret[i] = vertices(i)->p;
		if (colorcorrect) {
			for(int j=3;j<DIM;j++)
				ret[i][j] /= colorscale;
		}
	}
	return ret;
}
	

void awarpmesh::warp(vec *D, trimesh *m, bool colorscalecorrected) {

	if (m==NULL) {
		if (!colorscalecorrected) {
			for(int i=0;i<vertices.length();i++)
				movevertex(i,vertices[i]->p+D[i]);
		} else {
			vec tv;

			for(int i=0;i<vertices.length();i++) {
				tv = D[i];
				for(int j=3;j<DIM;j++)
					tv[j] *= colorscale;
				movevertex(i,vertices[i]->p+tv);
			}
		}
	} else {
		face *f;
		vec t1,t2;
		vertex *v1,*v2;
		FP alpha[3],d2;
		vec p;

		for(int i=0;i<vertices.length();i++) {
			p = m->closestpoint(vertices[i]->p,d2,f,v1,v2,t1,t2,0);
			f->project(p,alpha);
			if (!colorscalecorrected) {
				movevertex(i, vertices[i]->p - p +
					f->unproject(alpha,
						f->vertices[0]->p+
						   D[f->vertices[0]->index],
						f->vertices[1]->p+
						   D[f->vertices[1]->index],
						f->vertices[2]->p+
						   D[f->vertices[2]->index]));
			} else {
				vec tv1,tv2,tv3;
				vec tmp;
				tv1 = D[f->vertices[0]->index];
				tv2 = D[f->vertices[1]->index];
				tv3 = D[f->vertices[2]->index];
				for(int j=3;j<DIM;j++) {
					tv1[j] *= colorscale;
					tv2[j] *= colorscale;
					tv3[j] *= colorscale;
				}
				movevertex(i, tmp = vertices[i]->p - p +
					f->unproject(alpha,
						f->vertices[0]->p+tv1,
						f->vertices[1]->p+tv2,
						f->vertices[2]->p+tv3));
			}
		}
	}
}

vec awarpmesh::warp(const vec &pt, vec *D, bool colorscalecorrected) {

	face *f;
	vec t1,t2;
	vertex *v1,*v2;
	FP alpha[3],d2;
	vec p,tv1,tv2,tv3;

	p = closestpoint(pt,d2,f,v1,v2,t1,t2,0);
	f->project(p,alpha);
	tv1 = D[f->vertices[0]->index];
	tv2 = D[f->vertices[1]->index];
	tv3 = D[f->vertices[2]->index];
	if (colorscalecorrected) {
		for(int j=3;j<DIM;j++) {
			tv1[j] *= colorscale; tv2[j] *= colorscale; tv3[j] *= colorscale;
		}
	}
	return pt - p + f->unproject(alpha,f->vertices[0]->p+tv1,
			f->vertices[1]->p+tv2,f->vertices[2]->p+tv3);
}
