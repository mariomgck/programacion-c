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
#include <warping/warpmesh.h>
#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>

warpmesh::warpmesh(const cyberreader &cr, int subsamplesize, const warpparms &params) :
			awarpmesh(cr,subsamplesize,params) {
	setdata();
}

warpmesh::warpmesh(istream &s, const warpparms &params, FP cs, const matrix &A, const vec &b) :
			awarpmesh(s,params,cs,A,b) {
	setdata();
}

warpmesh::warpmesh(const trimesh &from, const warpparms &params, FP cs,
		const matrix &A, const vec &b) :
	awarpmesh(from,params,cs,A,b) {
	setdata();
}

warpmesh::~warpmesh() {
	delete []vel;
	delete []oldvel;
	delete []springs;
}

void warpmesh::setdata() {

	vel = new vec[vertices.length()];
	oldvel = new vec[vertices.length()];
	for(int i=0;i<vertices.length();i++)  {
		vel[i] = 0;
		oldvel[i] = 0;
	}
	springs = new ilist<int>[vertices.length()];
	awarpmesh::setdata();
}

void warpmesh::adduserspring(const vec &from, const vec &to) {
	int p = extrasprings.length();
	extrasprings.setlength(extrasprings.length()+1);
	extrasprings[p].to = to;
	vec t1,t2,p2;
	vertex *v1,*v2;
	FP d2;
	p2 = closestpoint(from,d2,extrasprings[p].from,v1,v2,t1,t2,0);
	extrasprings[p].from->project(p2,extrasprings[p].alphas);
}

void warpmesh::derivative(face *f, vertex *v1, vertex *v2, const vec &from,
				const vec &to, FP scale) {

	vec p1(to-from);
	
	if (v1==NULL && v2==NULL) { // the general case, the point lies in the triangle
		FP alpha[3];
		f->project(from,alpha);
		for(int j=0;j<3;j++)
			vel[f->vertices[j]->index] += (scale*alpha[j])*p1;
	} else if (v2==NULL) { // just one vertex
		vel[v1->index] += scale*p1;
	} else { // point is on an edge
		vec p2(from-v1->p);
		FP dd = p2.len2();
		p2 = v2->p - v1->p;
		FP d2 = p2.len2();
		if (d2==0) {
			vel[v2->index] += scale*p1;
			return;
		}
		dd = sqrt(dd/d2);
		vel[v2->index] += (dd*scale)*p1;
		vel[v1->index] += ((1-dd)*scale)*p1;
	}
}

void warpmesh::iterate(trimesh *target) {

	FP dot,d,dold;

	for(int i=0;i<vertices.length();i++) {
		oldvel[i] = vel[i];
		vel[i] *= momentum;
	}
	push(target);
	pull(target);
	spring();
	FP mm2=0.0,t;
	dot = 0.0; d = 0.0; dold = 0.0;
	for(int i=0;i<vertices.length();i++) {
		dot += vel[i]*oldvel[i];
		d += vel[i]*vel[i];
		dold += oldvel[i]*oldvel[i];
		if ((t=vel[i].len2())>mm2) mm2 = t;
	}
	cout << "(" << sqrt(d) << ")(" << dot/(sqrt(d)*sqrt(dold)) << ")";
	if (mm2>maxmove2/(epsilon*epsilon)) {
		mm2 = sqrt(maxmove2/(epsilon*epsilon*mm2));
		cout << "{" << mm2 << "}";
		for(int i=0;i<vertices.length();i++)
			vel[i] *= mm2;
	} else cout << "{1}";
	cout << endl;
	for(int i=0;i<vertices.length();i++)
		movevertex(i,vertices[i]->p+(epsilon*vel[i]));
	alpha *= eta;
	zeta *= xi;
}

void warpmesh::push(trimesh *target) {

	face *f,*of;
	vec p,p2;
	vec t1,t2;
	vertex *v1,*v2;
	FP d2;

	FP min,max,ttl;

	min = hugenum;
	max = -hugenum;
	ttl = 0;
	for(int i=0;i<n;i++) {
		p = samplepoint(f,n/4);
		if (rho>0)
			f->tangents(t1,t2);
		p2 = target->closestpoint(p,d2,of,v1,v2,t1,t2,rho);
		derivative(f,NULL,NULL,p,p2,1);
		ttl += d2;
		if (d2<min) min=d2;
		if (d2>max) max=d2;
	}
	cout << "[" << sqrt(ttl) << "]";
}

void warpmesh::pull(trimesh *target) {

	face *f,*of;
	vec p,p2;
	vec t1,t2;
	vertex *v1,*v2;
	FP d2;

	FP min,max,ttl;

	min = hugenum;
	max = -hugenum;
	ttl = 0;
	for(int i=0;i<n;i++) {
		p = target->samplepoint(f,n/4);
		if (rho>0)
			f->tangents(t1,t2);
		p2 = closestpoint(p,d2,of,v1,v2,t1,t2,rho);
		derivative(of,v1,v2,p2,p,1);
		ttl += d2;
		if (d2<min) min=d2;
		if (d2>max) max=d2;
	}
	cout << "[" << sqrt(ttl) << "]";
}

void warpmesh::spring() {

	vec sum,origv,dis;
	
	FP ttl,min,max,l;

	ttl = 0;
	min = hugenum;
	max = -hugenum;
	for(int i=0;i<vertices.length();i++) {
		sum = 0.0;
		for(int j=0;j<springs[i].length();j++) {
			
			origv = origpos[springs[i][j]]-origpos[i];
			dis = vertices[springs[i][j]]->p - vertices[i]->p;

			dis -= origv;
			l = origv.len(); // new -- to mimic old code
			sum += dis/l;

			l = dis.len2()/l;
			l *= alpha*alpha/
			   ((FP)springs[i].length()*(FP)springs[i].length());
			ttl += l;
			
		}
		vel[i] += (alpha/(FP)springs[i].length())*sum;
	}
	cout << "[" << sqrt(ttl) << "] ";
	ttl = 0;
	face *f;
	for(int i=0;i<extrasprings.length();i++) {
		f = extrasprings[i].from;
		dis = extrasprings[i].to - f->unproject(extrasprings[i].alphas);
		ttl += dis.len2()*zeta*zeta;
		for(int j=0;j<3;j++)
			vel[f->vertices[j]->index] += dis*(zeta*extrasprings[i].alphas[j]);
	}
	cout << "{" << extrasprings.length() << "}" ;
	cout << '[' << sqrt(ttl) << "] ";
}

void warpmesh::addsprings(ilist<int> *springs) {
	for(int i=0;i<vertices.length();i++)
		this->springs[i] += springs[i];
}
