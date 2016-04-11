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
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <matrix/vec.h>
#include "globaldef.h"
#include "face.h"
#include "vertex.h"
#include <iostream>

#include <containers/ilist.h>

/*class edge;
class reducer;
class progmesh;
class quadric;
class trimesh;
*/
 
class vertex  
{
/*
friend class trimesh;
friend class edge;
friend class face;
friend class reducer;
friend class quadric;
friend class progmesh;
friend class awarpmesh;
friend class warpmesh;
friend class warpmesh2;
friend class warpmesh3;
friend class warpmesh4;
friend class invshape;
*/

public:

	// number of edges leaving this vertex with an angle with a cosine of
	// less than angle
	int numcreases(FP angle=0.3) const;
	// number of edges leaving this vertex which are truly a ledge
	// (on the edge) -- they have only one face
	int numledges() const;

	// the position of the vertex
	inline vec position() const { return p; }

	int findedge(const vertex * const v) const;

	//static int count;

	// needed so that it can be used in a geohash
	inline vec max() const { return p; }
	inline vec min() const { return p; }

	inline int getindex() const { return index; }

	// signal that the point has moved
	inline void pointmoved() const {
		for(int i=0;i<faces.length();i++) faces(i)->pointsmoved();
	}

	inline FP area() const { 
		FP ret = 0;
		for(int i=0;i<faces.length();i++)
			ret += faces(i)->area();
		return ret;
	}

	inline FP sqrtarea() const {
		FP ret = 0;
		for(int i=0;i<faces.length();i++)
			ret += faces(i)->area();
		return sqrt(ret);
	}


//private:
	vertex() : edges(6,1), faces(6,1) { /*count++;*/ }
		// cout<<"vertices: "<<count<<endl;} ;
	virtual ~vertex() {/*count--;*/} ;

	vec p;

	ilist<edge *> edges;
	ilist<face *> faces;
	int index;
};

#endif
