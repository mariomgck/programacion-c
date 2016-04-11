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

#include <morphing/quadric.h>

quadric::quadric(trimesh *m, bool exactplacement) :
	A( (matrix *)malloc( sizeof(matrix)*m->vertices.length() ),
	                   m->vertices.length() ,
			   m->vertices.length() ),

	b((vec *)malloc(sizeof(vec)*m->vertices.length()),
	                m->vertices.length(),m->vertices.length()),

	c((FP *)malloc(sizeof(FP)*m->vertices.length()),
	               m->vertices.length(),m->vertices.length()) {

	this->m = m;
	exact = exactplacement;
	
	for(int i=0; i < m->vertices.length(); i++)
	{
		A[i] = 0;
		b[i] = 0;
		c[i] = 0;
	}

	matrix currA;
	vec currb;
	FP currc;
	int index;

	for(int i=0;i<m->faces.length();i++) 
	{
		facetoquadric(m->faces[i], currA, currb, currc);
		for(int j=0;j<3;j++) {
			index = m->faces[i]->vertices[j]->index;
			A[index] += currA;
			b[index] += currb;
			c[index] += currc;
		}
	}

	for(int i=0;i<m->edges.length();i++) {
		if (m->edges[i]->faces.length()<2) {
			edgetoquadric(m->edges[i],currA,currb,currc);
			index = m->edges[i]->vertices[0]->index;
			A[index] += currA;
			b[index] += currb;
			c[index] += currc;
			index = m->edges[i]->vertices[1]->index;
			A[index] += currA;
			b[index] += currb;
			c[index] += currc;
		}
	}		
	
	processlimit = m->edges.length() < 2000 ? m->edges.length()/4 : 500;
	initqueue();
}

quadric::~quadric(){
	delete q;
}

void quadric::initqueue() {
	togo *forqueue = (togo *)malloc(sizeof(togo)*m->edges.length());
	int count = 0;

	if (forqueue==NULL) {
		cout << "could not allocate memory for queue" << endl;
		cout << "needed space: " << sizeof(togo)*m->edges.length()
			<< endl;
		exit(2);
	}
	for(int i=0;i<m->edges.length();i++) {
		forqueue[count] = consider(i);
		if (forqueue[count].e != NULL) count++;
	}
	q = new pqueue<togo>(forqueue,count,m->edges.length());
}

bool quadric::reduce() {
	togo nextedge;
	if (!getnextedge(nextedge)) return false;


	// record pre information
	int prenumf = m->faces.length();
	int prenume = m->edges.length();

	vertex *v = nextedge.e->vertices[1];

	A.del(v->index);
	b.del(v->index);
	c.del(v->index);

	v = nextedge.e->vertices[0];

	A[v->index] = nextedge.newA;
	b[v->index] = nextedge.newb;
	c[v->index] = nextedge.newc;

	// remove edges from the queue that will disappear forever
	
	edge *e;
	int mi = nextedge.e->faces.length();
	for(int i=0;i<mi;i++) {
		for(int l=0;l<3;l++) {
			e = nextedge.e->faces[i]->edges[l];
			if (e==nextedge.e) continue;
			if (e->extra<0) {
				toprocess.del(-e->extra-1);
				if (-e->extra-1<toprocess.length())
					toprocess[-e->extra-1]->extra=e->extra;
			} else if (e->extra>0) q->remove(e->extra-1);
			e->extra = 0;
		}
	}
	
	// remove edge from mesh
	m->removeedge(nextedge.e->index,nextedge.newpos);

	edge *ed;
	vertex *ov;
	
	for(int i=0;i<v->edges.length();i++) {
		ed = v->edges[i];
		if (ed->extra>0) {
			q->remove(ed->extra-1);
			toprocess += ed;
			ed->extra = -toprocess.length();
		} else if (ed->extra==0) {
			toprocess += ed;
			ed->extra = -toprocess.length();
		}
		if (ed->vertices[0]==v) ov = ed->vertices[1];
		else ov = ed->vertices[0];
		for(int j=0;j<ov->edges.length();j++) {
			if (ov->edges[j]!=ed && ov->edges[j]->extra==0) {
				toprocess += ov->edges[j];
				ov->edges[j]->extra = -toprocess.length();
			}
		}
	}

	// phew... finally done
	return true;
}


int quadric::reduce(int n) {
	int i;
	for(i=0;i<n&&reduce();i++);
	return i; 
}

int quadric::complexity() { 
	return m->vertices.length(); 
}

void quadric::facetoquadric(face *f, matrix &retA, vec &retb, FP &retc) {
	vec t1,t2;

	if (!f->tangents(t1,t2)) {
		retA = matrix::zero;
		retb = 0;
		retc = 0;
	} else {
		FP area = f->area();
		vec p = f->basepoint();
		FP d1,d2;
		retA = (matrix::eye - matrix(t1,t1) - matrix(t2,t2))*area;
		d1 = p*t1;
		d2 = p*t2;
		retb = (d1*t1 + d2*t2 - p)*area;
		retc = (p*p - d1*d1 - d2*d2)*area;
	}
}

void quadric::edgetoquadric(edge *e, matrix &retA, vec &retb, FP &retc) {

	vec t1,t2;
	face *f = e->faces[0];

	if (!f->tangents(t1,t2)) {
		retA = matrix::zero;
		retb = 0;
		retc = 0;
	} else {
		FP len = e->len();
		vec p = e->vertices[0]->p;
		vec e1 = e->direction(e->vertices[0]);
		FP d1,d2,d3;
		d1 = p*t1;
		d2 = p*t2;
		d3 = p*e1;
		retA = (matrix(t1,t1) + matrix(t2,t2) - matrix(e1,e1))*len;
		retb = (d3*e1 - d1*t1 - d2*t2)*len;
		retc = (d1*d1 + d2*d2 - d3*d3)*len;
	}
}

bool quadric::getnextedge(togo &next) {

	bool b;

	if (toprocess.length() > processlimit)
		processedges();
	if (q->size()==0) return false;
	do {
		do { next=q->head(); q->remove(0);
		} while(!(b=removeable(m->edges[next.e->index],next.newpos)) && q->size()>0);
		// note that if we remove one (since it was not a valid
		// position for an edge (it has self-intersecting or reduced
		// past a tetrahedron) we do not have to recalculate it and
		// put it back on the queue.  When one of its neighbors is
		// moved, it will automatically be put back on the queue and
		// if it were put back on sooner, we would only have to remove
		// it again because as its surroundings are the same, the
		// vertex position would be the same and the position still
		// invalid. (this may not quite be true for quadrics --
		// although it is for progressive meshes -- however we are
		// going with it anyway)
		if (!b && toprocess.length()>0) processedges();
	} while (!b && q->size()>0);
	return b;
}

void quadric::processedges() {

	togo nextedge;
	FP pre,post;

	pre = q->head().deltaenergy;
	for(int i=0;i<toprocess.length();i++) {
		nextedge = consider(toprocess[i]->index);
		if (nextedge.e != NULL)
			q->add(nextedge);
		else toprocess[i]->extra = 0;
	}
	toprocess.setlength(0);
	post = q->head().deltaenergy;
	if (post==pre) processlimit++;
	else if (pre>0 && post<pre/2 && processlimit>1) processlimit /= 2;
}

quadric::togo quadric::consider(int e) {

	togo ret;

	if (m->edges[e]->faces.length()==0) {
		ret.e = NULL;
		return ret;
	}
	ret.e = m->edges[e];
	int v1,v2;

	v1 = m->edges[e]->vertices[0]->index;
	v2 = m->edges[e]->vertices[1]->index;
	ret.newA = A[v1] + A[v2];
	ret.newb = b[v1] + b[v2];
	ret.newc = c[v1] + c[v2];

	bool subset;
		bool worked;
	if (exact) {
		ret.newpos = -ret.newA.solve(ret.newb,worked);
		//matrix ainv = ret.newA.inv(worked);
		if (!worked) {
			subset=true;
		} else {
			//ret.newpos = -(ainv*ret.newb);
			ret.deltaenergy = ret.newc + ret.newb*ret.newpos;
			subset=false;
		}
	}
	if (!exact || subset) {
		ret.newpos = m->vertices[v1]->p;
		ret.deltaenergy = ret.newpos*(ret.newA*ret.newpos) +
		                    2*(ret.newb*ret.newpos) + ret.newc;
		vec testp = m->vertices[v2]->p;
		FP testd = testp*(ret.newA*testp) + 2*(ret.newb*testp) +
		                         ret.newc;
		if (testd<ret.deltaenergy) {
			ret.deltaenergy = testd;
			ret.newpos = testp;
		}
	}
	return ret;
}
