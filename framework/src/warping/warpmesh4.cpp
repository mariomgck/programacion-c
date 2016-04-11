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
#include <warping/warpmesh4.h>
#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>
#include <morphing/robust.h>

warpmesh4::warpmesh4(const cyberreader &cr, int subsamplesize, const warpparms &params) :
			awarpmesh(cr,subsamplesize,params) {
	setdata();
}

warpmesh4::warpmesh4(istream &s, const warpparms &params, FP cs, const matrix &A, const vec &b) :
			awarpmesh(s,params,cs,A,b) {
	setdata();
}

warpmesh4::~warpmesh4() {
	delete springCC;
	delete springCq;
	delete userCC;
	delete userCq;
	delete smoothCC;
	delete ptCC;
	delete ptCq;
}

warpmesh4::warpmesh4(const trimesh &from, const warpparms &params, FP cs,
		const matrix &A, const vec &b) :
			awarpmesh(from,params,cs,A,b) {
	setdata();
}

void warpmesh4::addconstraint(smatrix2 *CC, svec *Cq,
	const matrix &Ti, int nels, int *ind, FP *ai, const vec &qi) {

	matrix Tt(Ti.t()); // Tt = Ti'
	matrix TT(Tt*Ti);  // TT = Ti' * Ti
	vec Tq(Tt*qi);     // Tq = Ti' * qi
	
	int i,j,k,l;

	for(i=0;i<nels;i++) for(j=0;j<nels;j++)
		for(k=0;k<DIM;k++) for(l=0;l<DIM;l++) {
			CC->add(DIM*ind[i]+k,DIM*ind[j]+l,
				ai[i]*ai[j]*TT[k][l]);
		}
	for(i=0;i<nels;i++)
		for(k=0;k<DIM;k++)
			Cq->x[DIM*ind[i]+k] += ai[i]*Tq[k];
}
	
void warpmesh4::addspring(int i1, int i2) {
	int i[2];
	i[0] = i1; i[1] = i2;

	FP a[2];
	a[0] = sqrt((origpos[i1]-origpos[i2]).len());
	if (a[0]==0.0) {
		a[0] = 1e6;
	} else {
		a[0] = 1.0/a[0];
	}
	a[1] = -a[0];
	vec diff(origpos[i1]-origpos[i2]);
	addconstraint(springCC,springCq,matrix::eye,2,i,a,diff*a[0]);
}

void warpmesh4::addsmoothspring(int i1, int i2) {
	int i[2];
	i[0] = i1; i[1] = i2;
	FP a[2];
	a[0] = (FP)(vertices[i1]->edges.length()+
	                vertices[i2]->edges.length())/
	           (vertices[i1]->edges.length()*
			vertices[i2]->edges.length());
	a[1] = -a[0];
	addconstraint(smoothCC,springCq,matrix::eye,2,i,a,vec::zero);
}
	
void warpmesh4::addsprings(ilist<int> *springs) {

	int ttlsp = 0;
	for(int i=0;i<vertices.length();i++) ttlsp += springs[i].length();
	ttlsp /= 2; // we've counted each one twice
	//cout << "total number of springs: " << ttlsp << endl;
	springCC = new smatrix2(DIM*vertices.length(),DIM*vertices.length());
	springCC->clear();
	userCC = new smatrix2(DIM*vertices.length(),DIM*vertices.length());
	userCC->clear();
	smoothCC = new smatrix2(DIM*vertices.length(),DIM*vertices.length());
	smoothCC->clear();
	ptCC = new smatrix2(DIM*vertices.length(),DIM*vertices.length());
	usersprings = false;
	//cout << "matrix created" << endl;
	springCq = new svec(DIM*vertices.length());
	springCq->clear();
	userCq = new svec(DIM*vertices.length());
	userCq->clear();
	ptCq = new svec(DIM*vertices.length());
	int l1,j;
	for(int i=0;i<vertices.length();i++) {
		l1 = springs[i].length();
		for(j=0;j<l1;j++) {
			if (springs[i](j)<i) continue;
			addspring(i,springs[i](j));
		}
	}
	for(int i=0;i<edges.length();i++) {
		addsmoothspring(edges[i]->vertices[0]->index,
		                edges[i]->vertices[1]->index);
	}
}

void warpmesh4::adduserspring(const vec &from, const vec &to) {
	FP alphas[3];
	face *f;
	vec t1,t2,p2;
	vertex *v1,*v2;
	FP d2;
	int in[3];

	usersprings = true;
	p2 = closestpoint(from,d2,f,v1,v2,t1,t2,0);
	f->project(p2,alphas);
	in[0] = f->vertices[0]->index;
	in[1] = f->vertices[1]->index;
	in[2] = f->vertices[2]->index;
	
	addconstraint(userCC,userCq,matrix::eye,3,in,alphas,to);
}

void warpmesh4::removesprings() {
	userCC->clear();
	userCq->clear();
	usersprings = false;
}

void warpmesh4::setsystem(trimesh *target) {

	int i,d;
	vec p,p2,t1,t2;
	vertex *v1,*v2;
	face *f,*of;
	FP alphas[3],d2;
	int in[3];
	FP factor = (FP)1.0/faces.length();

	ptCC->clear();
	ptCq->clear();
	matrix Ti;
	if (maxmove2<0) {
		for(i=0;i<3;i++) alphas[i] = 0;
		for(i=0;i<vertices.length();i++) {
			p = vertices[i]->p;
			p2 = target->closestpoint(p,d2,of,v1,v2);
			if (momentum<0) Ti = matrix::eye;
			else {
			if (v1==NULL) {
				of->tangents(t1,t2);
				Ti = matrix::eye-matrix(t1,t1)-matrix(t2,t2);
			} else if (v2==NULL) {
				Ti = matrix::eye;
			} else {
				t1 = target->getedge(v1->findedge(v2))->
							direction(v1);
				Ti = matrix::eye-matrix(t1,t1);
			}
			}

			alphas[0] = 1;
			in[0] = i;
			for(d=0;d<3;d++) {
				alphas[d] *= factor;
			}
			addconstraint(ptCC,ptCq,Ti,1,in,alphas,factor*Ti*p2);
		}
	}
	for(i=0;i<n;i++) {
		if (sampleface.length()<=i) {
			p = samplepoint(f,n/4);
			sampleface += f;
			f->project(p,alphas);
			samplealpha1 += alphas[0];
			samplealpha2 += alphas[1];
			samplealpha3 += alphas[2];
		} else {
			f = sampleface[i];
			alphas[0] = samplealpha1[i];
			alphas[1] = samplealpha2[i];
			alphas[2] = samplealpha3[i];
			p = f->unproject(alphas);
		}
		if (rho>0) {
			f->tangents(t1,t2);
		}
		p2 = target->closestpoint(p,d2,of,v1,v2,t1,t2,rho);
		if (momentum<0) Ti = matrix::eye;
		else {
		if (v1==NULL) {
			of->tangents(t1,t2);
			Ti = matrix::eye - matrix(t1,t1) - matrix(t2,t2);
		} else if (v2==NULL) {
			Ti = matrix::eye;
		} else {
			t1 = target->getedge(v1->findedge(v2))->
						direction(v1);
			if (!t1.isvalid()) Ti = matrix::eye;
			else Ti = matrix::eye - matrix(t1,t1);
		}
		}

		for(d=0;d<3;d++) {
			in[d] = f->vertices[d]->index;
			alphas[d] *= factor;
		}
		addconstraint(ptCC,ptCq,Ti,3,in,alphas,factor*Ti*p2);
	}
	bool tokay;
	for(i=0;i<n;i++) {
		if (sampletface.length()<=i) {
			p = target->samplepoint(f,n/4);
			sampletface += f;
			samplevec += p;
		} else {
			p = samplevec[i];
			f = sampletface[i];
		}
		if (rho>0) {
		tokay = f->tangents(t1,t2);
		}
		p2 = closestpoint(p,d2,of,v1,v2,t1,t2,rho);
		if (momentum<0) Ti = matrix::eye;
		else {
		if (tokay)
			Ti = matrix::eye - matrix(t1,t1) - matrix(t2,t2);
		else
			Ti = matrix::eye;
		}
		of->project(p2,alphas);
		for(d=0;d<3;d++) {
			in[d] = of->vertices[d]->index;
			alphas[d] *= factor;
		}
		addconstraint(ptCC,ptCq,Ti,3,in,alphas,factor*Ti*p);
	}
}

void warpmesh4::iterate(trimesh *target) {

	double len2 = -1.0;
	double baselen2 = -1.0;
	sampleface.setlength(0);
	sampletface.setlength(0);
	samplealpha1.setlength(0);
	samplealpha2.setlength(0);
	samplealpha3.setlength(0);
	samplevec.setlength(0);
	do {
		baselen2 = len2;
		// set the ptCC and ptCq to a new set of sampled points
		setsystem(target);

		svec *w = new svec(DIM*vertices.length());
		int d,i,c;
		vec vp;
		for(i=c=0;i<vertices.length();i++) {
			vp = vertices[i]->p;
			for(d=0;d<DIM;d++,c++)
				w->x[c] = vp[d];
		}

		ptCC->addmult(*springCC,alpha);
		ptCq->addmult(*springCq,alpha);
		ptCC->addmult(*smoothCC,epsilon);
		// smoothCq is the implied zero vector
		if (usersprings) {
			ptCC->addmult(*userCC,zeta);
			ptCq->addmult(*userCq,zeta);
		}
		
		ptCC->solve(*ptCq,*w,1e-5);

		len2 = 0.0;
		for(i=c=0;i<vertices.length();i++) {
			vp = vertexpos(i);
			for(d=0;d<DIM;d++,c++) {
				double temp = vp[d]-w->x[c];
				len2 += temp*temp;
				vp[d] = w->x[c];
			}
			movevertex(i,vp);
		}
		delete w;
		//if (baselen2<0.0) baselen2 = len2/100;
		cout << "len2 = " << len2 << endl;
	} while(baselen2<0.0 || len2<baselen2);
	//} while (len2>baselen2);
	sampleface.setlength(0);
	sampletface.setlength(0);
	samplealpha1.setlength(0);
	samplealpha2.setlength(0);
	samplealpha3.setlength(0);
	samplevec.setlength(0);

	alpha *= eta;
	zeta *= xi;
}
