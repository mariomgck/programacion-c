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
#ifndef _WARPMESH3_H_
#define _WARPMESH3_H_

#include <geometry/trimesh.h>
#include <iostream>
#include <warping/warpparms.h>
#include <matrix/sparse.h>
#include <warping/awarpmesh.h>

using namespace std;

// This is a subclass of awarpmesh designed to add functionality needed to
// move the mesh around according to minimizing the energy function detailed
// in Shelton's MS Thesis.
 
// NOTE:  removeedge must not be called on a warpmesh!  This could be fixed
// later, but I haven't done so yet
class warpmesh3 : public awarpmesh
{
public:
	// Load a new warpmesh from a cyberware scan
	// The subsamplesize refers to the compression per side needed --
	// see trimesh
	warpmesh3(const cyberreader &cr, int subsamplesize, const warpparms &params);

	// Load a new mesh from a stream (this mesh should have been saved
	// with the trimesh save method with type trimesh::tmff1)
	// cs is the color to shape ratio to be applied after it is loaded
	// and the entire mesh is to be tranformed by the affine
	// transformation defined by A and b: x <- Ax + b
	warpmesh3(istream &s, const warpparms &params, FP cs=1.0,
		const matrix &A = matrix::eye, const vec &b = vec::zero);

	warpmesh3(const trimesh &from, const warpparms &params, FP cs=1.0,
		const matrix &A = matrix::eye, const vec &b = vec::zero);

	virtual ~warpmesh3();

	// This adds a spring from the point from on the mesh to another
	// point "to."  "to" is the user defined correspondence for the
	// point "from."
	virtual void adduserspring(const vec &from, const vec &to);
	// removes all springs from "adduserspring."
	virtual void removesprings() { nusp=0; }

	// performs one iteration of matching the warpmesh to the target viai
	// gradient descent on the energy function in Shelton's MS Thesis.
	virtual void iterate(trimesh *target);

public:
	void addsprings(ilist<int> *springs);

private:

	int nusp;

	smatrix *A;
	svec *b[DIM];
	svec *slack[DIM];

	virtual void setsystem(trimesh *target);
};

#endif
