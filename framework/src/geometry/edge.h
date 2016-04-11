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
#ifndef _EDGE_H_
#define _EDGE_H_

#include <matrix/vec.h>
#include "vertex.h"
//#include "trimesh.h"
#include "globaldef.h"
#include <iostream>
#include "face.h"

/*class reducer;
class progmesh;
class quadric;
*/
class edge  
{
/*friend class trimesh;
friend class vertex;
friend class face;
friend class reducer;
friend class quadric;
friend class progmesh;
friend class awarpmesh;
friend class warpmesh;
friend class warpmesh2;
friend class warpmesh3;
friend class warpmesh4;
*/

public:
	// returns if the edge is a crease (considering only the first two
	//  connected faces) given critangle as the maximum cos of the angle
	//  for a crease.
	inline bool iscrease(FP critangle = 0.3) const {
		return faces.length()<=1 || creaseangle()<=critangle;
	}

	inline bool is3dcrease(FP critangle = 0.3) const {
		return faces.length()<=1 || faces(0)->normal()*faces(1)->normal()<=critangle;
	}

	inline bool isend() const {
		return faces.length()<=1;
	}

	// find the cos of the crease angle
	// replaces the connected face ff with rf
	// considers rv1 and rv2 to not be part of either adjacent face
	FP creaseangle(face *ff=NULL, face *rf=NULL, vertex *rv1=NULL, vertex *rv2=NULL) const;

	// returns the closest point to the edge in ret.
	// a zero return code indicated this point is on the edge
	// a one return code means that the returned point is the first vertex
	// a two for a return code means that the returned point is the second vertex
	int closestpoint(const vec &p, vec &ret) const;

	// returned the area using the two vertices of the edge and the given
	// point p as corners of a triangle (this is more needed/useful than one
	// might think!)
	FP area(const vec &p) const;

	// the midpoint
	inline vec mid(void) const {
		vec m;
		m = (vertices[0]->p + vertices[1]->p)/2;
		return m;
	}
 
	// the wieghted midpoint
	inline vec mid(FP frac) const {
		vec m;
		m = frac*vertices[1]->p + (1-frac)*vertices[0]->p;
		return m;
	}

	// the square of the length
	inline FP len2() const {
		return (vertices[0]->p-vertices[1]->p).len2();
	}

	// the length
	inline FP len() const {
		return (vertices[0]->p-vertices[1]->p).len();
	}

	// do the two vertices match the endpoints of this edge?
	inline bool match(const vertex * const v1, const vertex * const v2) {
		return (v1==vertices[0]&& v2==vertices[1]) || (v1==vertices[1] && v2==vertices[0]);
	}

	// given one of the vertices, return a unit vector in the direction from
	// that vertex to the other
	inline vec direction(vertex *v) {
		vec dir;

		if (v==vertices[0]) dir = vertices[1]->p - vertices[0]->p;
		else dir = vertices[0]->p - vertices[1]->p;
		return dir.norm();
	}

	inline int vertexindex(int i) {
		return vertices[i]->index;
	}

	// this is an extra integer to keep track of whatever you want (well, provided
	//  no one else is using it :) )
	int extra; // I hate to add this, but I need this for the reducer and this is the best
	// place to store the extra integer (the index in the pqueue).  This means that in use outside
	// of the reducer, I'm wasting these four bytes and I really killing the abstraction
	// barrier.  But, keeping a list of these in the reducer class would be a pain since the
	// index of the edges keeps changing.

//private:

	inline edge() : faces(2,1) { /*count++;*/ } // cout << "edges: " << count << endl; }
	inline ~edge() { /*count--;*/}

	vertex *vertices[2];
	ilist<face *> faces;
	int index;
	

};

#endif
