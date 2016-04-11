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

#include <morphing/reducer.h>

#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>

static const FP SELFINT = (FP)(-0.7);

bool reducer::removeable(edge *ed, const vec &to) {

	// first check that we won't be reducing past a tetrahedron
	int j=0;
	vertex *v1 = ed->vertices[0];
	vertex *v2 = ed->vertices[1];
	int mi=v1->faces.length();
	ilist<face *> f1;
	for(int i=0;i<mi;i++)
		if (v2->faces.find(v1->faces[i])==-1)
			f1 += v1->faces[i];
		else j++;
	if (j!=2 && j!=1) {
		//printf ("missed on 1\n");
		return false;
	}
	if (f1.length()==0) {
		//printf ("missed on 2\n");
		return false;
	}
	j=0;
	mi=v2->faces.length();
	ilist<face *> f2;
	for(int i=0;i<mi;i++)
		if (v1->faces.find(v2->faces[i])==-1)
			f2 += v2->faces[i];
		else j++;
	if (j!=2 && j!=1) {
		//printf ("missed on 3\n");
		return false;
	}
	if (f2.length()==0) {
		//printf ("missed on 4\n");
		return false;
	}

	ilist<edge *> e1;
	mi=f1.length();
	for(int i=0;i<mi;i++)
		for(int j=0;j<3;j++)
			if(e1.find(f1[i]->edges[j])==-1)
				e1 += f1[i]->edges[j];
	ilist<edge *> e2;
	mi=f2.length();
	for(int i=0;i<mi;i++)
		for(int j=0;j<3;j++)
			if(e2.find(f2[i]->edges[j])==-1)
				e2 += f2[i]->edges[j];
	// this all comes down to the fact that the faces from
	// the first vertex should share no edges with the faces
	// from the second vertex with the exception of those edges
	// which are part of the triangles attached to both vertices
	for(int i=0;i<e1.length();i++)
		for(int j=0;j<e2.length();j++)
			if (e1[i]==e2[j]) {
				//printf ("missed on 5\n");
				return false;
			}
	

	// Okay... now that that's done, we must determine whether
	// any of the dihedral angles are too small (self-intersection)
	ilist<edge *> e3;
	face *f;
	mi = ed->faces.length();
	for(int i=0;i<mi;i++) {
		f = ed->faces[i];
		for(int j=0;j<3;j++)
			if (e3.find(f->edges[j])==-1)
				e3 += f->edges[j];
	}
	vec store1=v1->p,store2=v2->p;
	v1->p = to; v2->p = to;
	mi=e1.length();
	for(int i=0;i<mi;i++)
		if (e3.find(e1[i])==-1 && e1[i]->faces.length()>1 &&
			e1[i]->iscrease(SELFINT)) {
			v1->p = store1;
			v2->p = store2;
			//printf ("missed on 6\n");
			return false;
		}
	mi=e2.length();
	for(int i=0;i<mi;i++)
		if (e3.find(e2[i])==-1 && e2[i]->faces.length()>1 &&
			e2[i]->iscrease(SELFINT)) {
			v1->p = store1;
			v2->p = store2;
			//printf ("missed on 7\n");
			return false;
		}
	mi=ed->faces.length();
	face *ff2;
	edge *ee1,*ee2;
	// we have to do a complex switch on the edges and such for the
	// two pairs of edges which will merge since we don't want to 
	// actually perform the topology change yet
	for(int i=0;i<mi;i++) {
		f = ed->faces[i];
		if (f->edges[1]==ed) {
			ee2 = f->edges[2];
			ee1 = f->edges[0];
		} else if (f->edges[2]==ed) {
			ee2 = f->edges[0];
			ee1 = f->edges[1];
		} else {
			ee2 = f->edges[1];
			ee1 = f->edges[2];
		}
		if (ee1->faces.length()<2 || ee2->faces.length()<2) continue;
		
		for(int j=0;j<ee2->faces.length();j++)
			if(ee2->faces[j]!=f) { ff2 = ee2->faces[j]; break; }
		if (j==ee2->faces.length()) continue;
		if (ee1->faces.length()>1 &&
			ee1->creaseangle(f,ff2,v1,v2)<SELFINT) {
			v1->p = store1;
			v2->p = store2;
			//printf ("missed on 8\n");
			return false;
		}
	}
	v1->p = store1;
	v2->p = store2;
	//printf ("succeeded\n");
	return true;

}
