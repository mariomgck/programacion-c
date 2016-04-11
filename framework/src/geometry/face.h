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
#ifndef _FACE_H_
#define _FACE_H_

#include <matrix/vec.h>
#include "globaldef.h"
#include <iostream>

class edge;
//class trimesh;
//class reducer;
//class progmesh;
//class quadric;
class vertex;

class face {
/*friend class trimesh;
friend class vertex;
friend class edge;
friend class reducer;
friend class quadric;
friend class progmesh;
friend class awarpmesh;
friend class warpmesh;
friend class warpmesh2;
friend class warpmesh3;
friend class warpmesh4;
friend class invshape;*/

public:

	// the middle of the triangle
	vec midpoint() const;
	// a random point sampled uniformly from the triangle
	vec randompoint();

	// finds the closest point on the triangle to the point p.  Returns
	//  whether it found it (this is needed since it uses ceil2 as the
	// largest possible value -- if the answer *must* be greater than
	// small2 it terminates the algorithm in mid-computation and returns
	// 0).  dist2 is the squared distance.
	int closestpoint(const vec &p, vec& ret, const FP &ceil2, FP &dist2, 
		vertex*& v1, vertex*& v2, const vec &t1, const vec &t2,
		const FP &anglefact = 0);

	// finds the point of intersection of a ray and the face
	// (note that ceil2 is measured in units of v -- if v isn't of unit
	// length this can make a difference)
	bool rayintersection(const vec &x0, const vec v, const FP &ceil,
		FP &dist, vec &ret);
	// and now the version ignoring color and pretending everything is
	//  in 3d (the last dim-3 dimensions of each vector are ignored)
	bool rayintersection3d(const vec &x0, const vec v, const FP &ceil,
		FP &dist, vec &ret) const;

	// the area of the triangle
	FP area() const;

	// returns two orthonormal vectors that are tangent to the triangle
	bool tangents(vec &t1, vec &t2) {
		if (tvalid) {
			t1 = t1cache;
			t2 = t2cache;
			return tretcache;
		} else return realtangents(t1,t2);
	}

	// a point on the triangle (to be use with the above to establish
	// a plane)
	vec basepoint();

	// normal returns the R3 normal of the triangle projected into
	// R3 (the last dim-3 components of its vertices' positions are
	// set to zero)
	vec normal() const;
 
	// the corner of the axis aligned bounding box with maximal coordinates
	// (used by the geohash)
	inline vec max() const {
		return mymax;
	}
	// smae as above, but with the mimunum coordinates
	inline vec min() const {
		return mymin;
	}

	// is the given point (which must lie in the plane of the triangle)
	// inside the triangle?
	bool ptinface(const vec &p) const;
	// same as above, but only considering the projection onto the
	// first 3 axes
	bool ptinface3d(const vec &p) const;

	// given a point with only the first three coordinates filled in,
	// recontruct the other three
	vec addcolor(const vec &x);

	// It is necessary that the bounds be cached so that they can be
	// updated as needed by trimesh::movevertex as the faces are moved
	// one-by-one in the geohash structure.
	void updatebounds();

	// takes a vector p and places in alpha (a pointer to an array of
	// 3 FP's) the coordinates of that point with respect to the vertices
	// (the triangle centered coordinate system)
	void project(const vec &p, FP *alpha) const;
		// changes alpha (alpha is the return value)
	// unproject the projection of above
	vec unproject(const FP *alpha) const;

	// I don't like having this version here, but I'm not sure exactly
	// where to put this function -- it is the same as the unproject
	// above, except that we replace the vertices' positions with p1, p2,
	// and p3. (so really it doesn't reply on any data members of the face)
	static vec unproject(const FP *alpha, const vec &v1, const vec &v2,
		const vec &v3); 

	// tells the triangle that the points have moved and it should
	// invalidate its cache
	inline void pointsmoved() { tvalid = false; }

	int vertexindex(int i) const;

//private:
	bool realtangents(vec &t1, vec &t2);
	// only meshes can make new triangles
	inline face() { tvalid = false; gone=false; }
	inline ~face() { }

	inline FP tripmax(FP a, FP b, FP c) const {
		return (a>b)?((a>c)?a:c):((b>c)?b:c);
	}
	inline FP tripmin(FP a, FP b, FP c) const {
		return (a<b)?((a<c)?a:c):((b<c)?b:c);
	}

	// the code is so tailored to triangles, there is no sense in
	// making this general
	edge *edges[3];
	vertex *vertices[3];
	int index;
	bool gone;
	vec mymin,mymax;
	vec t1cache,t2cache;
	FP u0cache,v0cache,r0cache;
	bool tretcache;
	bool tvalid;

};

#endif
