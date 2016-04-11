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
#ifndef _PROGMESH_H_
#define _PROGMESH_H_

#include <iostream>
#include <geometry/trimesh.h>
#include <geometry/edge.h>
#include <containers/pqueue.h>
#include <morphing/reducer.h>

// implements Hoppe's progressive meshes (see SIGGRAPH '96)
// (well, at the moment, it only reduces, it does not perform the reconstruction
//  which would require a few changes -- note that the reduce function would have
//  to return the lost information instead of a bool)
// If I need to I'll implement that later (let's hope I don't).

class progmesh : public reducer  
{
public:
	// from the time the reducer is created until it is destroyed, it has 
	//  control over the mesh given
	progmesh(trimesh *mesh, FP cpen=0, int extrapts=0);
	// the mesh is *not* destroyed here
	virtual ~progmesh();

	// removes one edge (which implies removing a vertex, two faces, and another two edges)
	// returns false if no edges can be reduced
	virtual bool reduce();

	// removes n edges and returns the number actually reduced
	virtual int reduce(int n);

	// returns a measure of complexity
	virtual int complexity();

private:
	void initqueue();

	FP creasepenalty;

	trimesh *m;
	trimesh *origm;

	ilist<edge *> toprocess;
	int processlimit;

	int *pfaces;
	int cpface;

	void processedges();

	class pt {
	public:
		vec p;
		int face; // face on the mesh onto which it projects
		FP d; // square of the distance from the point to that face
	};
	
	pt *pts; // points sampled from the original mesh
	int numpts;

	class finfo { // face information
	public:
		ilist<int> pts; // points which project onto this face
	};


	finfo *ptsback; // back pointing list of points

	class togo { // everything you need to know to remove an edge
	public:
		edge *e;
		FP deltaenergy;
		vec newpos;
		int index; // in pqueue
		int newpface;

		// needed just to satisfy the ilist template
		inline bool operator==(const togo &a) { return e==a.e; }

		inline FP rank() { return -deltaenergy; }
		inline void setposition(int i) { e->extra = i+1; }
	};

	pqueue<togo> q;
	
	bool getnextedge(togo &next);
	togo consider(int e); // returns the best position for removing edge e
	void findpos(int e, const vec &start, vec &end,
		FP &energy, const ilist<face *> &faces, const ilist<face *> &ofaces, const FP &K);
	FP pull(vertex *v1, vertex *v2, const ilist<face *> &faces, const ilist<face *> &ofaces,
		vec &pos);
	FP push(vertex *v1, vec &pos);
	FP spring(vertex *v1, vertex *v2, vec &pos, const FP &K);
	FP springconstant(const ilist<face *> &faces);

	FP findenergy(vertex *v1, vertex *v2, const ilist<face *> &faces, 
		const ilist<face *> &ofaces, const FP &K);
	FP pullenergy(vertex *v1, vertex *v2, const ilist<face *> &faces, const ilist<face *> &ofaces);
	FP pushenergy(vertex *v1);
	FP springenergy(vertex *v1, vertex *v2, const FP &K);
	void getfacelist(int e, ilist<face *> &faces, ilist<face *> &ofaces);
	
	bool creasetopologychange(edge *e, FP cangle);

	vec neighborhoodclosestpoint(const vec &p, FP &best, const ilist<face *> &faces,
									  face *&fret, vertex *&v1, vertex *&v2) const;

	bool closestpoint(const ilist<face *> &flist, const vec &p, 
		face *&fret, FP &dret);

	FP creaselengthchange(int e, const vec &p, const FP &cangle, const ilist<face *> &nhood);

	void pushproject(vec p, vec &retp, FP &d2);
};

#endif
