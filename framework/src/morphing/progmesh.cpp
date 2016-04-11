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
#include <globaldef.h>

#include <geometry/face.h>
#include <geometry/edge.h>
#include <geometry/vertex.h>

#include <containers/ilist.h>

#include <matrix/vec.h>
#include <morphing/progmesh.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

progmesh::progmesh(trimesh *mesh, FP cpen, int extrapts)
{
	creasepenalty = cpen;
	m = mesh;
	origm = mesh->dup();
	numpts = m->vertices.length()+extrapts;
	pts = new pt[numpts];
	ptsback = new finfo[m->faces.length()];
	pfaces = new int[m->vertices.length()];
	int j=0,mi=m->vertices.length();
	for(int i=0;i<mi;i++) {
		if (m->vertices[i]->faces.length()>0) {
			pts[j].p = m->vertices[i]->p;
			pts[j].d = 0;
			pts[j].face = m->vertices[i]->faces[0]->index;
			ptsback[pts[j].face].pts += j; // add to the reverse lookup too
			pfaces[i] = pts[j].face;
			j++;
		} else { // this means we have a vertex not associated
			// with any face (I'm guessing we should try to 
			// move the mesh to capture this vertex, although
			// what to do isn't exactly clear)
			/*
			face *f;
			vertex *v1,*v2;
			vec t1,t2;

			cout << "having to find closest point (why?)" << endl;
			m->closestpoint(pts[i].p,pts[i].d,f,v1,v2,t1,t2,0);
			pts[i].face = f->index;
			*/
			// okay, I'm just going to leave the point out instead.
			numpts--;
		}
	}
	face *f;

	for(int i=0;i<extrapts;i++) {
		pts[j].p = m->samplepoint(f,extrapts);
		pts[j].face = f->index;
		pts[j].d = 0;
		ptsback[f->index].pts += j;
		j++;
	}
	processlimit = m->edges.length() < 2000 ? m->edges.length()/4 : 500;
	initqueue();
}

progmesh::~progmesh()
{
	delete[] pts;
	delete[] ptsback;
}

int progmesh::reduce(int n) {
	int i;
	for(i=0;i<n&&reduce();i++);
	return i; 
}

int progmesh::complexity() { 
	return m->vertices.length(); 
}

// start dret off with a ceiling if you want to use one (otherwise be sure it is set to
//  hugenum)
bool progmesh::closestpoint(const ilist<face *> &flist, const vec &p, face *&fret, FP &dret) {

	FP curr;
	vec r;
	vertex *v1,*v2;
	vec t1,t2;
	bool ret = false;

	for(int n=0;n<flist.length();n++) {
		if (flist(n)->closestpoint(p,r,dret,curr,v1,v2,t1,t2,0) && curr<dret) {
			dret = curr;
			fret = flist(n);
			ret = true;
		}
	}
	return ret;
	/*
	face *f;
	m->closestpoint(p,curr,f,v1,v2,t1,t2,0);
	if (curr<dret) {
		dret = curr;
		fret = f;
		return true;
	} else return false;
	*/
}

bool progmesh::reduce() {
	togo nextedge;
	if (!getnextedge(nextedge)) return false;

	cout << "removing: " << nextedge.e->index << " (" << nextedge.deltaenergy << ") [" <<
		q.size() << ']';

	// record pre information
	int prenumf = m->faces.length();
	int prenume = m->edges.length();

	ilist<int> adjf;
	int mi=nextedge.e->faces.length();
	for(int i=0;i<mi;i++)
		adjf += nextedge.e->faces[i]->index;

	vertex *v = nextedge.e->vertices[0];

	// remove edges from the queue that will disappear forever
	mi = adjf.length();
	int k;
	edge *e;
	for(int i=0;i<mi;i++) {
		for(int l=0;l<3;l++) {
			e = m->faces[adjf[i]]->edges[l];
			if (e==nextedge.e) continue;
			if (e->extra<0) {
				toprocess.del(-e->extra-1);
				if (-e->extra-1 < toprocess.length())
					toprocess[-e->extra-1]->extra = e->extra;
			} else if (e->extra>0) q.remove(e->extra-1);
			e->extra = 0;
		}
	}

	// adjust the push pface cache
	pfaces[v->index] = nextedge.newpface;
	pfaces[nextedge.e->vertices[1]->index] = pfaces[m->vertices.length()-1];
	
	// remove edge from mesh
	m->removeedge(nextedge.e->index,nextedge.newpos);

	ilist<int> addpts; // points that need to be recomputed from scratch

	// some faces were renumbered so we need to update the pts and ptsback lists
	// accordingly
	mi = adjf.length();
	int mj,j,tt;
	for(int i=0;i<mi;i++) {
		k = adjf[i];
		addpts += ptsback[k].pts;
		ptsback[k].pts.setlength(0);
		tt = m->faces.length()+mi-i-1;
		if (tt==k) continue;
		ptsback[k].pts.swap(ptsback[tt].pts);
		mj = ptsback[k].pts.length();
		for(int j=0;j<mj;j++) {
			pts[ptsback[k].pts[j]].face = k;
		}
		for(int j=i+1;j<mi;j++)
			if (adjf[j]==tt) adjf[j] = k;
	}

	// now we need to recompute the closest faces to some of the points
	// in particular, all of the points projecting onto the neighborhood of
	// v need to be recomputed.  As well, the ones sitting in addpts need to
	// be recomputed.  Also, we will recompute points projecting on faces
	// adjacent to the neighborhood of v.  In general, any point could have
	// changed, but following the lead of Hoppe, DeRose, and Duchamp in their
	// paper Mesh Optimization (see bottom of 8 / top of 9) we will assume that
	// a point doesn't move more than one face between changes.

	// we also need to record which edges are next to faces that changed so that 
	// we can pull them from the queue and readd them.

	ilist<face *> &nhood = v->faces; // the neighborhood
	ilist<face *> nnhood; // the "extended" neighborhood minus the neighborhood
	face *ff,*f;
	int l,n;
	edge *ll;

	// grab the "second order" neighborhood
	// (and add all edges adjacent to the "first order" neighborhood to the
	//  changed edge list)
	for(int i=0;i<nhood.length();i++) {
		f = nhood[i];
		for(int j=0;j<3;j++) {
			if (f->vertices[j]==v) continue;
			for(k=0;k<f->vertices[j]->faces.length();k++) {
				ff=f->vertices[j]->faces[k];
				if (nhood.find(ff)==-1 && nnhood.find(ff)==-1)
					nnhood += ff;
			}

			// the following has nothing to do with the above but needs to be done in
			// basically the same loop, so I put it here

			for(k=0;k<f->vertices[j]->edges.length();k++) {
				ll = f->vertices[j]->edges[k];
				if (ll->extra>0) {
					q.remove(ll->extra-1);
					toprocess += ll;
					ll->extra = -toprocess.length();
				} else if (ll->extra==0) {
					toprocess += ll;
					ll->extra = -toprocess.length();
				}
			}
		}
	}
	

	// update the closest faces to the points which projected into the neighborhood
	for(int i=0;i<nhood.length();i++) {
		k = nhood[i]->index;
		mj = ptsback[k].pts.length();
		for(int j=0;j<mj;j++) {
			l = ptsback[k].pts[j];
			pts[l].d = hugenum;
			{	vertex *v1,*v2; vec t1,t2;
				m->closestpoint(pts[l].p,pts[l].d,f,v1,v2,t1,t2,(FP)0.0);
			}
			//closestpoint(nhood,pts[l].p,f,pts[l].d);
			//closestpoint(nnhood,pts[l].p,f,pts[l].d);
			if (pts[l].face != f->index) {
				pts[l].face = f->index;
				ptsback[k].pts.del(j);
				mj--; j--;
				ptsback[f->index].pts += l;
			}
		}
	}

	// do the same for the points which projected onto the faces which were removed
	for(int i=0;i<addpts.length();i++) {
		k = addpts[i];
		pts[k].d = hugenum;
		//closestpoint(nhood,pts[k].d,f,pts[k].d);
		//closestpoint(nnhood,pts[k].d,f,pts[k].d);
		{ vertex *v1,*v2; vec t1,t2;
			m->closestpoint(pts[k].p,pts[k].d,f,v1,v2,t1,t2,(FP)0.0);
		}
		pts[k].face = f->index;
		ptsback[pts[k].face].pts += k;
		if (nhood.find(f)==-1) { // else it has already been added
			for(int j=0;j<3;j++) {
				for(int l=0;l<f->vertices[j]->edges.length();l++) {
					ll = f->vertices[j]->edges[l];
					if (ll->extra>0) {
						q.remove(ll->extra-1);
						toprocess += ll;
						ll->extra = -toprocess.length();
					} else if (ll->extra==0) {
						toprocess += ll;
						ll->extra = -toprocess.length();
					}
				}
			}
		}
	}
	// do the same for the extended neighborhood
	//for(i=0;i<nnhood.length();i++) {
	//	k = nnhood[i]->index;
	//	mj = ptsback[k].pts.length();
	//	for(j=0;j<mj;j++) {
	//		l = ptsback[k].pts[j];

	for(int l=0;l<numpts;l++) {
		k= pts[l].face;
		if (nhood.find(m->faces[k])==-1 && addpts.find(k)==-1) {

			if (closestpoint(nhood,pts[l].p,f,pts[l].d)) {
				ptsback[k].pts.remove(l);
				ptsback[f->index].pts += l;
				pts[l].face = f->index;
				int b,mb;
				for(int n=0;n<3;n++) {
					mb = m->faces[k]->vertices[n]->edges.length();
					for(b=0;b<mb;b++) {
						ll = m->faces[k]->vertices[n]->edges[b];
						if (ll->extra>0){
							q.remove(ll->extra-1);
							toprocess += ll;
							ll->extra = -toprocess.length();
						} else if (ll->extra==0) {
							toprocess += ll;
							ll->extra = -toprocess.length();
						}
					}
				}
			}
		}
	}
	

	cout << '|' << toprocess.length() << '|' << endl;
	
	// phew... finally done
	return true;
}


void progmesh::processedges() {

	togo nextedge;
	FP pre,post;

	pre = q.head().deltaenergy;
	for(int i=0;i<toprocess.length();i++) {
		nextedge = consider(toprocess[i]->index);
		if (nextedge.e != NULL)
			q.add(nextedge);
		else toprocess[i]->extra = 0;
	}
	toprocess.setlength(0);
	post = q.head().deltaenergy;
	if (post==pre) processlimit++;
	else if (pre>0 && post<pre/2 && processlimit>1) processlimit /= 2;
}

void progmesh::initqueue() {
	togo nextedge;
	for(int i=0;i<m->edges.length();i++) {
		cout << "init queue: " << i << '/' << m->edges.length() << endl;
		nextedge = consider(i);
		if (nextedge.e != NULL)
			q.add(nextedge);
		else m->edges[i]->extra = 0;
	}
}

void progmesh::getfacelist(int e, ilist<face *> &faces, ilist<face *> &ofaces) {

	face *ff;
	face *f1 = m->edges[e]->faces[0];
	face *f2 = m->edges[e]->faces.length()>1 ? m->edges[e]->faces[1] : NULL;

	vertex *v = m->edges[e]->vertices[0];
	int mi=v->faces.length();
	for(int i=0;i<mi;i++)
		faces += v->faces[i];

	vertex *v2 = m->edges[e]->vertices[1];
	mi=v2->faces.length();
	for(int i=0;i<mi;i++) {
		ff = v2->faces[i];
		if (ff!=f1 && ff!=f2) {
			faces += ff;
			ofaces += ff;
		}
	}

	face *f;
	int j,k;
	for(int i=0;i<faces.length();i++) {
		f = faces[i];
		for(j=0;j<3;j++) {
			if (f->vertices[j]==v || f->vertices[j]==v2) continue;
			for(k=0;k<f->vertices[j]->faces.length();k++) {
				ff=f->vertices[j]->faces[k];
				if (ofaces.find(ff)==-1)
					ofaces += ff;
			}
		}
	}
}

FP progmesh::springconstant(const ilist<face *> &faces) {
	
//	int n=0;
//	for(int i=0;i<faces.length();i++)
//		n += ptsback[faces(i)->index].pts.length();
//	n /= faces.length();
//	return n<2 ? (FP)0.1 : (n<4 ? (FP)0.01 : (n<8 ? (FP)0.001 : (FP)0.0001));
	return (FP)0.001;	
}

progmesh::togo progmesh::consider(int e) {

	FP energy;
	vec pos;
	togo ret;
	ilist<face *> faces;
	ilist<face *> ofaces;
	vec v1,v2;

	/*
	if (creasetopologychange(m->edges[e],(FP)-0.99)) { // changes the topology of the mesh
		ret.e = NULL;
		return ret;
	}
	*/
	if (m->edges[e]->faces.length()==0) {
		ret.e = NULL;
		return ret;
	}
	
	v1 = m->edges[e]->vertices[0]->p;
	v2 = m->edges[e]->vertices[1]->p;
	getfacelist(e,faces,ofaces);
	FP K = springconstant(faces);

	ret.e = m->edges[e];
	cpface = pfaces[m->edges[e]->vertices[0]->index];
	FP beginenergy = findenergy(m->edges[e]->vertices[0],m->edges[e]->vertices[1],faces,ofaces,K);
	cpface = pfaces[m->edges[e]->vertices[0]->index];
	findpos(e,v1,ret.newpos,ret.deltaenergy,faces,ofaces,K);
	ret.newpface = cpface;
	cpface = pfaces[m->edges[e]->vertices[1]->index];
	findpos(e,v2,pos,energy,faces,ofaces,K);
	if (energy < ret.deltaenergy) {
		ret.newpface = cpface;
		ret.newpos = pos;
		ret.deltaenergy = energy;
	}
	
	cpface = pfaces[m->edges[e]->vertices[0]->index];
	findpos(e,(v1+v2)/2,pos,energy,faces,ofaces,K);
	if (energy < ret.deltaenergy) {
		ret.newpface = cpface;
		ret.newpos = pos;
		ret.deltaenergy = energy;
	}
	
	m->edges[e]->vertices[0]->p = v1;
	m->edges[e]->vertices[1]->p = v2;
	m->edges[e]->vertices[0]->pointmoved();
	m->edges[e]->vertices[1]->pointmoved();
	//m->movevertex(m->edges[e]->vertices[0]->index,v1);
	//m->movevertex(m->edges[e]->vertices[1]->index,v2);
	ret.deltaenergy -= beginenergy;
	/*
	if (creasepenalty>0 && creasetopologychange(m->edges[e],(FP)0.3)) {
		ret.deltaenergy += m->edges[e]->len()*creasepenalty;
	}
	if (creasepenalty>0) {
		ret.deltaenergy += creasepenalty * creaselengthchange(e,ret.newpos,(FP)0.3,faces);
	}
	*/
	
	return ret;
}

const int ITTNUM=10;

void progmesh::findpos(int e, const vec &start, vec &end,
					  FP &energy, const ilist<face *> &faces, const ilist<face *> &ofaces,
					  const FP &K) {

	vec newpos;
	vertex *v1,*v2;
	FP den;

	end = start;
	int i=0;
	v1 = m->edges[e]->vertices[0];
	v2 = m->edges[e]->vertices[1];
	do {
		v1->p = end;
		v2->p = end;
		v1->pointmoved();
		v2->pointmoved();
		//m->movevertex(v1->index,end);
		//m->movevertex(v2->index,end);
		newpos = 0.0;
		den=pull(v1,v2,faces,ofaces,newpos);
		den+=push(v1,newpos);
		den+=spring(v1,v2,newpos,K);
		newpos /= den;
		if ((newpos-end).len2() < 0.00001 || i>10) {
			end = newpos; 
			v1->p = end;
			v2->p = end;
			v1->pointmoved();
			v2->pointmoved();
			//m->movevertex(v1->index,end);
			//m->movevertex(v2->index,end);
			energy=findenergy(v1,v2,faces,ofaces,K);
			return;
		}
		end = newpos;
		i++;
	} while(1);
}

vec progmesh::neighborhoodclosestpoint(const vec &p, FP &best, const ilist<face *> &faces,
									  face *&fret, vertex *&v1, vertex *&v2) const {

	FP curr;
	vec ret,cv;
	vertex *tv1,*tv2;
	vec t1,t2;

	best=hugenum;
	for(int i=0;i<faces.length();i++) {
		if (faces(i)->closestpoint(p,cv,best,curr,tv1,tv2,t1,t2,(FP)0.0) && curr<best) {
			best = curr;
			v1 = tv1;
			v2 = tv2;
			ret = cv;
			fret =faces(i);
		}
	}
	return ret;
}

FP progmesh::pull(vertex *v1, vertex *v2, const ilist<face *> &faces,
				 const ilist<face *> &ofaces, vec &pos) {

	vec p,t1,t2,pp;
	FP d2,td,ret=0;
	face *f;
	vertex *rv1,*rv2;

	int mj,j, k,mi=faces.length();
	for(int i=0;i<mi;i++) {
		mj = ptsback[faces(i)->index].pts.length();
		for(j=0;j<mj;j++) {
			pp = pts[ptsback[faces(i)->index].pts[j]].p;
			//p = m->closestpoint(pp,d2,f,rv1,rv2,t1,t2,0);
			p = neighborhoodclosestpoint(pp,d2,ofaces,f,rv1,rv2);
			if (faces.find(f)==-1) continue;
			if (rv1==NULL) {
				for(k=0;k<3;k++)
					if (f->edges[k]->vertices[0]!=v1
						&& f->edges[k]->vertices[0]!=v2
						&& f->edges[k]->vertices[1]!=v1
						&& f->edges[k]->vertices[1]!=v2)
							break;
				if(k==3) { // we have a point on one
					// of the two triangles ajoining the edge e (which
					// in *theory* can't happen but does in practice often)
					for(k=0;k<3;k++) 
						if(f->vertices[k]!=v1 && f->vertices[k]!=v2)
							break;
					t1 = f->vertices[k]->p;
				} else { // we have a "general" point on an adjacent triangle
					d2 = f->edges[k]->area(p)/f->area();
					if (!_finite(d2) || p==v1->p) { // the triangle has no area
						// or we picked the point exactly (without it being flagged somehow)
						ret += 1;
						pos += pp;
						continue;
					}
					t2 = lineintersection(f->edges[k]->vertices[0]->p,
						f->edges[k]->vertices[1]->p - f->edges[k]->vertices[0]->p,
						v1->p,p-v1->p);
					if (!t2.isvalid())
						t2 = f->edges[k]->vertices[0]->p;
					pos += d2*(pp-t2) + (d2*d2)*t2;
					ret += d2*d2;
					continue;
				}
			} else if (rv2!=NULL && (rv1==v1 || rv1==v2)) {
				t1 = rv2->p;
			} else if (rv2!=NULL && (rv2==v1 || rv2==v2)) {
				t1 = rv1->p;
			} else if (rv2==NULL && (rv1==v1 || rv1==v2)) {
				pos += pp;
				ret += 1;
				continue;
			} else continue;
			t2 = p-t1;
			d2 = t2*t2;
			t2 = v1->p - t1;
			td = t2*t2;
			if (td==0) { 
				pos += pp;
				ret += (FP)1.0;
			} else {
				d2 = sqrt(d2/td);
				pos += d2*(pp-t1) + (d2*d2)*t1;
				ret += d2*d2;
			}
		}
	}
	return ret;
}

void progmesh::pushproject(vec p, vec &retp, FP &d2) {

	face *f = origm->faces[cpface];
	vec t1,t2;
	vertex *v1,*v2;


	/*
	retp = origm->closestpoint(p,d2,f,v1,v2,t1,t2,0);
	cpface = f->index;
	*/
	
	FP best,curr;
	vec rp;
	best = hugenum;
	for(int i=0;i<3;i++)
		for(int j=0;j<f->vertices[i]->faces.length();j++)
			if (f->vertices[i]->faces[j]->closestpoint(p,rp,best,curr,v1,v2,t1,t2,0) &&
				curr<best) {
				best=curr;
				cpface = f->vertices[i]->faces[j]->index;
				retp = rp;
			}
	d2 = best;
	
}				

FP progmesh::push(vertex *v1, vec &pos) {
	return 0;
	/*
	FP d2;

	pushproject(v1->p,pos,d2);
	return 1.0;
	*/
}

FP progmesh::spring(vertex *v1, vertex *v2, vec &pos, const FP &K) {

	FP ret=0;
	int mi=v1->edges.length();
	for(int i=0;i<mi;i++) {
		if (v2==v1->edges[i]->vertices[0] ||
			v2==v1->edges[i]->vertices[1]) continue;
		if (v1==v1->edges[i]->vertices[0]) {
			pos += v1->edges[i]->vertices[1]->p*K;
			ret += K;
		} else {
			pos += v1->edges[i]->vertices[0]->p*K;
			ret += K;
		}
	}
	mi=v2->edges.length();
	for(int i=0;i<mi;i++) {
		if (v1==v2->edges[i]->vertices[0] ||
			v1==v2->edges[i]->vertices[1]) continue;
		if (v2==v2->edges[i]->vertices[0]) {
			pos += v2->edges[i]->vertices[1]->p*K;
			ret += K;
		} else {
			pos += v2->edges[i]->vertices[0]->p*K;
			ret += K;
		}
	}
	return ret;
}

FP progmesh::findenergy(vertex *v1, vertex *v2, const ilist<face *> &faces, 
					   const ilist<face *> &ofaces, const FP &K) {
	return pullenergy(v1,v2,faces,ofaces)+pushenergy(v1)+springenergy(v1,v2,K);
}

FP progmesh::pullenergy(vertex *v1, vertex *v2, const ilist<face *> &faces,
					   const ilist<face *> &ofaces) {

	FP ret = 0,d2;
	vec p,pp,t1,t2;
	vertex *rv1,*rv2;
	face *f;


	int mj,j,mi=faces.length();
	for(int i=0;i<mi;i++) {
		mj = ptsback[faces(i)->index].pts.length();
		for(j=0;j<mj;j++) {
			pp = pts[ptsback[faces(i)->index].pts[j]].p;
			//p = m->closestpoint(pp,d2,f,rv1,rv2,t1,t2,0);
			p = neighborhoodclosestpoint(pp,d2,ofaces,f,rv1,rv2);
			ret += d2;
		}
	}
	return ret;
}

FP progmesh::pushenergy(vertex *v1) {
	return 0;
	/*
	FP d;
	vec t;

	pushproject(v1->p,t,d);
	return d;
	*/
}

FP progmesh::springenergy(vertex *v1, vertex *v2, const FP &K) {

	FP ret=0;
	int mi=v1->edges.length();
	for(int i=0;i<mi;i++) {
		if (v2==v1->edges[i]->vertices[0] ||
			v2==v1->edges[i]->vertices[1]) continue;
		ret += K*v1->edges[i]->len2();
	}
	mi=v2->edges.length();

	for(int i=0;i<mi;i++) {
		if (v1==v2->edges[i]->vertices[0] ||
			v1==v2->edges[i]->vertices[1]) continue;
		ret += K*v2->edges[i]->len2();
	}
	return ret;
}


bool progmesh::getnextedge(togo &next) {

	bool b;

	if (toprocess.length() > processlimit)
		processedges();
	if (q.size()==0) return false;
	do {
		do { next=q.head(); q.remove(0);
		} while(!(b=removeable(m->edges[next.e->index],next.newpos)) && q.size()>0);
		// note that if we remove one (since it was not a valid position
		// for an edge (it has self-intersecting or reduced past a tetrahedron)
		// we do not have to recalculate it and put it back on the queue.
		// When one of its neighbors is moved, it will automatically be put back
		// on the queue and if it were put back on sooner, we would only have
		// to remove it again because as its surroundings are the same, the
		// vertex position would be the same and the position still invalid.
		if (!b && toprocess.length()>0) processedges();
	} while (!b && q.size()>0);
	return b;
}


bool progmesh::creasetopologychange(edge *e, FP cangle) {
	bool sh;
	int snum,tnum;
	int i,j;
	vertex *v1,*v2,*v3;
	face *f1,*f2;
	edge *esl,*etl,*etr,*esr;
	vertex *vl,*vr;

	
	if (e->faces.length()<1 || e->faces.length()>2) return true;
	v1 = e->vertices[0]; v2 = e->vertices[1];
	f1 = e->faces[0];
	if (e->faces.length()<2) {
		if (v1->numledges() != 2 ||
			v2->numledges() != 2) return true;
		for(i=0;i<3;i++)
			if (f1->vertices[i]!=v1 &&
				f1->vertices[i]!=v2) {
				v3 = f1->vertices[i];
				break;
			}
		ilist<vertex *> vs;
		for(i=0;i<v1->edges.length();i++)
			if (v1->edges[i]->vertices[0]==v1)
				vs += v1->edges[i]->vertices[1];
			else vs += v1->edges[i]->vertices[0];
		for(i=0;i<v2->edges.length();i++)
			if (v2->edges[i]->vertices[0]==v2) {
				if (v2->edges[i]->vertices[0]!=v3 &&
					vs.find(v2->edges[i]->vertices[0])!=-1) return true;
			} else {
				if (v2->edges[i]->vertices[1]!=v3 &&
					vs.find(v2->edges[i]->vertices[1])!=-1) return true;
			}
		tnum = v1->numcreases(cangle);
		snum = v2->numcreases(cangle);
		return ((snum>=3 && tnum>=3) ||
				(snum==1 && tnum!=2) ||
				(tnum==1 && snum!=2));
	}
	f2 = e->faces[1];
	for(i=0;i<3;i++) {
		if (f1->vertices[i] != v1 && f1->vertices[i] != v2) {
			vl = f1->vertices[i];
			if (i==2) j=0;
			else j=i+1;
			if (f1->vertices[j]==v1) {
				esl = f1->edges[i];
				if (i==0) j=2;
				else j=i-1;
				etl = f1->edges[j];
			} else {
				etl=f1->edges[i];
				if (i==0) j=2;
				else j=i-1;
				esl = f1->edges[j];
			}
			break;
		}
	}
	for(i=0;i<3;i++) {
		if (f2->vertices[i] != v1 && f2->vertices[i] != v2) {
			vr = f2->vertices[i];
			if (i==2) j=0;
			else j=i+1;
			if (f2->vertices[j]==v1) {
				esr = f2->edges[i];
				if (i==0) j=2;
				else j=i-1;
				etr = f2->edges[j];
			} else {
				etr = f2->edges[i];
				if (i==0) j=2;
				else j=i-1;
				esr = f2->edges[j];
			}
			break;
		}
	}
	tnum = v1->numcreases(cangle);
	snum = v2->numcreases(cangle);
	sh = e->iscrease(cangle);
	return ((esl->iscrease(cangle) && etl->iscrease(cangle)) ||
			(esr->iscrease(cangle) && etr->iscrease(cangle)) ||
			(!sh && snum>=1 && tnum>=1) ||
			(sh && snum>=3 && tnum>=3) ||
			(sh && snum==1 && tnum!=2) ||
			(sh && tnum==1 && snum!=2));
}

FP progmesh::creaselengthchange(int e, const vec &p, const FP &cangle, const ilist<face *> &nhood) {

	//return (FP)0.0;
	
	edge *e1,*e2,*e3,*e4,*ed;
	face *f;
	ilist<bool> creases;
	ilist<edge *> eds;
	bool b1,b2,b3,b4;

	int j;
	for(int i=0;i<nhood.length();i++) {
		f=nhood(i);
		for(j=0;j<3;j++) {
			if (eds.find(f->edges[j])==-1) {
				eds += f->edges[j];
				creases += f->edges[j]->iscrease(cangle);
			}
		}
	}
	ed = m->edges[e];
	f = ed->faces[0];
	if (f->edges[0]==ed) {
		e1 = f->edges[1];
		e2 = f->edges[2];
	} else if (f->edges[1]==ed) {
		e1 = f->edges[2];
		e2 = f->edges[0];
	} else {
		e1 = f->edges[0];
		e2 = f->edges[1];
	}
	b1 = e1->iscrease(cangle);
	b2 = e2->iscrease(cangle);
	if (ed->faces.length()<2) {
		e3 = e4 = NULL;
		b3 = b4 = false;
	} else {
		f = ed->faces[1];
		if (f->edges[0]==ed) {
			e3 = f->edges[1];
			e4 = f->edges[2];
		} else if (f->edges[1]==ed) {
			e3 = f->edges[2];
			e4 = f->edges[0];
		} else {
			e3 = f->edges[0];
			e4 = f->edges[1];
		}
		b3 = e3->iscrease(cangle);
		b4 = e4->iscrease(cangle);
	}
	FP ret = 0.0;
	vec s1,s2;
	s1 = ed->vertices[0]->p;
	s2 = ed->vertices[1]->p;
	ed->vertices[0]->p = p;
	ed->vertices[1]->p = p;
	for(int i=0;i<eds.length();i++) 
		if (creases[i] != eds[i]->iscrease(cangle))
			ret += eds[i]->len();
	f = ed->faces[0];
	bool check;
	if (e1->faces.length()>1 && e2->faces.length()>1) {
		if (e2->faces[0] == f)
			check = e1->creaseangle(f,e2->faces[1],ed->vertices[0],ed->vertices[1])<=cangle;
		else check = e1->creaseangle(f,e2->faces[0],ed->vertices[0],ed->vertices[1])<=cangle;
		if (check != (b1 || b2)) ret += (e1->len() + e2->len())/2;
	}
	if (e3!=NULL && e3->faces.length()>1 && e4->faces.length()>1) {
		if (e4->faces[0] == f)
			check = e3->creaseangle(f,e4->faces[1],ed->vertices[0],ed->vertices[1])<=cangle;
		else check = e3->creaseangle(f,e4->faces[0],ed->vertices[0],ed->vertices[1])<=cangle;
		if (check != (b3 || b4)) ret += (e3->len() + e4->len())/2;
	}
	ed->vertices[0]->p = s1;
	ed->vertices[1]->p = s2;
	return ret;
	
}

