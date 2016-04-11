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
#ifndef _QUADRIC_H_
#define _QUADRIC_H_

#include <morphing/reducer.h>

#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>
#include <geometry/trimesh.h>

#include <containers/ilist.h>
#include <containers/pqueue.h>
#include <containers/dlist.h>

// quadric implements Garland's and Heckbert's Quadric Error Metrics
// method for reducing polygonal models.

class quadric : public reducer  
{
public:
	
	quadric(trimesh *m, bool exactplacement=true);

	virtual ~quadric();

	virtual bool reduce();

	virtual int reduce(int n);

	virtual int complexity();

private:

	class togo { // everything you need to know to remove an edge
	public:
		edge *e;
		FP deltaenergy;
		vec newpos;
		//int index; // in pqueue

		matrix newA;
		vec newb;
		FP newc;
		
		// needed just to satisfy the ilist template
		inline bool operator==(const togo &a) { return e==a.e; }

		inline FP rank() { return -deltaenergy; }
		inline void setposition(int i) { e->extra = i+1; }
	};

	void initqueue();
	bool getnextedge(togo &next);
	void processedges();
	togo consider(int e);

	void facetoquadric(face *f, matrix &retA, vec &retb, FP &retc);
	void edgetoquadric(edge *e, matrix &retA, vec &retb, FP &retc);
	trimesh *m;

	bool exact;

	dlist<matrix> A;
	dlist<vec> b;
	dlist<FP> c;

	pqueue<togo> *q;
	dlist<edge *> toprocess;
	int processlimit;
};

#endif
