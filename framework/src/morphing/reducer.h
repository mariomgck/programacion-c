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
#ifndef _REDUCER_H_
#define _REDUCER_H_

#include <matrix/vec.h>
#include <geometry/edge.h>

// this is the superclass of all polygon reducers (progmesh and quadric)
class reducer  
{
public:

	virtual ~reducer() {};

	// reduce by one edge (later will make this return the information needed
	// to undo the edge collapse
	virtual bool reduce() =0;

	// reduce by n edges
	virtual int reduce(int n) =0;

	// return a measure of the complexity of the surface (usually the number of vertices)
	virtual int complexity() =0;

protected:
	// this method is provided for the subclasses to determine if an edge can be
	// reduced without changing the topology of the mesh (or producing self-intersections)
	// It is not exact (such a check for self-intersections would be very time consuming)
	// but works well
	bool removeable(edge *ed, const vec &to);
};

#endif
