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
#ifndef _TRIMESH_H_
#define _TRIMESH_H_

#include <iostream>
#include "globaldef.h"

#include <geometry/face.h>
#include <geometry/vertex.h>
#include <geometry/trirender.h>
#include <geometry/edge.h>

#include <matrix/vec.h>


#include <loaders/cyberreader.h>
#include <loaders/ppmreader.h>
#include <image/image.h>

#include <containers/ilist.h>
#include <containers/dlist.h>
#include <containers/geohash.h>

/*class reducer;
class progmesh;
class quadric;
*/

class vecconvert {
public:
	virtual void convertpt(vec in, double &x, double &y, double &z,
		double &r, double &g, double &b) = 0;
};

class straightvecconvert : public vecconvert {
public:
	virtual void convertpt(vec in, double &x, double &y, double &z,
                double &r, double &g, double &b) {
		x = in[0%DIM];
		y = in[1%DIM];
		z = in[2%DIM];
		r = in[3%DIM];
		g = in[4%DIM];
		b = in[5%DIM];
	}
};

// a general class for keeping track of meshes of triangles of any topology
class trimesh  
{
/*friend class reducer;
friend class progmesh;
friend class quadric;
friend class invshape;
*/
public:
	// create a blank mesh with the given bounding box (divsize is ignored)
	// hashon sets whether to create a geometric hashing structure for
	// queries (if it is set to be off, then you cannot call closestpoint
	// point or castray)
	trimesh(const vec &minpt, const vec &maxpt, const vec &divsize = 0,
		bool hashon = true);

	// a new trimesh from a cyberware scan (the reader is passed in)
	trimesh(const cyberreader &cr, int subsamplesize=1, bool hashon = true);

	// a new mesh from a ppm image
	trimesh(const ppmreader &pr, int subsamplesize=1, const FP &scale = 1,
		bool hashon = true);

	typedef enum {tmff1, vrml1, oldff} fileformat;
	// save the file in one of the previous file formats
	// (tmff1 = the file format for this class
	//  vrml1 = VRML 1.0 (cannot be read -- only written )
	//  oldff = the old mesh file formet)

	// a new mesh loaded from a stream
	trimesh(istream &s, FP cs = 1.0, bool hashon = true,
		const matrix &A = matrix::eye, const vec &b = vec::zero,
		const fileformat ft=tmff1);

	// copy constructor
	trimesh(const trimesh &from, FP cs = 1.0, bool hashon = true,
		const matrix &A = matrix::eye, const vec &b = vec::zero);

	virtual ~trimesh();
	
	virtual trimesh *dup() const; // this may not work on all compilers
		// but it should -- if you get warnings this may have to be
		// adjusted to return through the parameter list instead
		// see comment on warpmesh.h

	// returns the vector of the closest point on the mesh (not the
	// closest vertex -- although it might end up being a vertex).  This
	// can only be used if the mesh was created with hashon set to true.
	// The distance from p to the returned point, squared, is returned in
	// dist2.  fret is the face on which the returned point lies.  v1 and
	// v2 are set to indicate whether the returned point is a vertex (v1
	// is the vertex, v2 is NULL), on an edge (v1 and v2 are the vertices
	// of the edge) or on a triangle (v1 and v2 are both NULL).
	// t1 and t2 are tangent vectors to a plane if you which to match
	// orientation (with the weight of anglefact -- see Shelton's MS
	// thesis:  anglefact is the same as rho in that thesis)
	// maxdist2 can be set to stop the search at a certain radius (well,
	// the radius is the sqrt of maxdist2)
	vec closestpoint(const vec &p, FP& dist2, face*& fret,
		vertex* &v1, vertex*& v2, const vec &t1=vec::zero,
		const vec &t2=vec::zero, const FP &anglefact=0,
		const FP &maxdist2=hugenum) const;

	// this finds the closest *vertex* to p and returns it along with the
	// squared distance in dist2 (maxdist2 again can be used to set a
	// maximum search radius)
	vertex *closestvertex(const vec &p, FP &dist2,
		FP maxdist2=hugenum) const;

	// the next two functions cast a ray and find the first point of
	// intersection.  The ray is p + av for any positive values of a.  It
	// returns whether anything was hit and ret gives the position of the
	// intersection with fret giving the face of interesection.  maxd
	// gives a maximal a (from the eq. above) for the search.
	bool castray(const vec &p, const vec &v, vec &ret,
		face *&fret, FP &retd, const FP &maxd=hugenum) const;
	// same as above, but *everything* is projected onto the first 3 axes
	bool castray3d(const vec &p, const vec &v, vec &ret, face *&fret,
		FP &retd, const FP &maxd=hugenum, face *hint=NULL) const;
	
	// removes the edge given by edgenum and moves the resulting vertex
	// to newpos
	void removeedge(int edgenum, const vec &newpos);

	// moves a vertex
	void movevertex(int vertexnum, const vec &newpos);
	// returns the position of a vertex
	vec vertexpos(int vertexnum) const;
	inline vertex *vertexnum(int vertexn) const
		{ return vertices(vertexn); }

	// adds a face connecting the given points (and adds any needed
	// edges and vertices)
	face *addface(const vec &p1, const vec &p2, const vec &p3);

	// on a trial basis -- these may get moved back to the private
	// section
	face *addface(int v1, int v2, int v3);
	int addvertex(const vec &pos);

	// removes a face and any edges that are left unattached.  It does
	// *not* remove edges or points which are left "unattached" if the
	// corresponding argument is set to false
	void removeface(int fnum, bool removeedge=true, bool removept=true);
	// This one removes a face connecting the three vertices listed
	//   if it exists
	void removeface(int v1, int v2, int v3,
		bool removeedge=true, bool removept=true);

	// removes the entire connected component.  Otherwise just as
	// removeface.
	void removecomponent(int fnum, bool removeedge=true,
			bool removept=true);
	void removecomponent(face *f, bool removeedge=true,
			bool removept=true);

	// numpts should be set to roughly 1/2 or so of the number of
	//  points you wish to sample.  This allows the algorithm to be
	//  faster if it knows you want a set of samples and are only
	//  considering whether the collection is random and not whether 
	//  their ordering is random.  Put differently, only collections
	//  of about numpts points are random.  The ordering within a set of
	//  numpts points is not random, but if considered as a set (not
	//  a list), they are random.
	// If you don't want to recompute the areas (regardless of whether the
	//  cached values are correct), set checkvalid to false.

	// Since I have a cache, this can't be const -- which may turn out to
	// be a pain and require me to make the cache a pointer off of the
	// class later.
	vec samplepoint(face *& fret, int numpts=1, bool checkvalid=true);
 
	// these operate on pointers since in order to create a mesh you must
	// already know the bounding box (which isn't strictly true anymore,
	// but I've had no need to change this)
	friend ostream &operator<<(ostream &s, const trimesh * const m);
	friend istream &operator>>(istream &s, trimesh *&m);

	void saveit(ostream &s, fileformat ft=tmff1) const;

	// the number of vertices (if n is the value, the valid indices for
	//    vertices are 0 through n-1)
	inline int numvertices() { return vertices.length(); }

	// make a virtual cyberware scan by ray casting --
	//    returns an image with channels red, green, blue, radius.
	image cyberwarescan(FP miny, FP maxy,
			int numlat=256,int numlon=512) const;

	inline void extremes(vec &minp, vec &maxp) {
		minp = vertexgrid.getmin(); maxp = vertexgrid.getmax();
	}
	void realextremes(vec &minp, vec &maxp);

	// this function completes a mesh in which the "color" sides were
	// not included (like one read in from an inventor file)
	void colorcomplete(void);

	inline FP getcolorscale(void) { return colorscale; }

	FP distance(trimesh *m, int npts);
	FP maxdistance(trimesh *m);
	FP mindistance(trimesh *m);
	FP width();

	static straightvecconvert svc;

	void render(trirender &renderer,bool split=true, 
		vecconvert *vc = &(trimesh::svc), FP ambient=0.25,
		bool colorized=true, FP direct=0.75);

	inline const face *facenum(int fnum) const { return faces(fnum); }
	inline int numfaces() const { return faces.length(); }

	inline void printlevelstats() const {
		if (facegrid.numq > 0 && facegrid.numq%1000 == 0) {
			for(int i=0;i<=facegrid.levellimit;i++)
				 cout << i << ": " << facegrid.counts[i] << endl;
			cout << facegrid.numq << '/' << facegrid.nums << endl;
		}
	}

	edge *getedge(int i) const { return edges(i); }
	inline int numedges() const { return edges.length(); }

public:
	dlist<face *> faces;
	dlist<edge *> edges;
	dlist<vertex *> vertices;

	bool geohashing;

	geohash<face *> facegrid;
	geohash<vertex *> vertexgrid;
	FP colorscale;

public:

	void finishloading(istream &s, const matrix &A = matrix::eye,
		const vec &b = vec::zero, const fileformat = tmff1);

	// a bunch of stuff for keeping track of point sampling and
	// caching areas

	bool areasvalid;
	ilist<FP> areas;
	ilist<FP> areas10;
	ilist<FP> areas100;
	int areaface;
	FP arearem,totalarea;
};

#endif
