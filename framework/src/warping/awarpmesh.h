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
#ifndef _AWARPMESH_H_
#define _AWARPMESH_H_


#include <iostream>
#include <geometry/trimesh.h>
#include <warping/warpparms.h>

// This is a subclass of trimesh designed to add functionality needed to
// move the mesh around according to minimizing the energy function detailed
// in Shelton's MS Thesis.
 
// NOTE:  removeedge must not be called on a warpmesh!  This could be fixed
// later, but I haven't done so yet
class awarpmesh : public trimesh  
{
public:
	// Load a new warpmesh from a cyberware scan
	// The subsamplesize refers to the compression per side needed --
	// see trimesh
	awarpmesh(const cyberreader &cr, int subsamplesize,
		const warpparms &params);

	// Load a new mesh from a stream (this mesh should have been saved
	//     with the trimesh save method with type trimesh::tmff1)
	// cs is the color to shape ratio to be applied after it is loaded and
	// the entire mesh is to be tranformed by the affine transformation
	// defined by A and b
	// x <- Ax + b
	awarpmesh(istream &s, const warpparms &params, FP cs=1.0,
		const matrix &A = matrix::eye, const vec &b = vec::zero);
	awarpmesh(const trimesh &from, const warpparms &params, FP cs=1.0,
		const matrix &A = matrix::eye, const vec &b = vec::zero);

	virtual ~awarpmesh();

	// This adds a spring from the point from on the mesh to another
	// point "to."  "to"
	// is the user defined correspondence for the point "from."
	virtual void adduserspring(const vec &from, const vec &to) = 0;
	// removes all springs from "adduserspring."
	virtual void removesprings() = 0;

	// returns the displacements of the vertices (the callee now owns
	// the returned memory)
	virtual vec *D(bool colorcorrect=true) const;
	// returns the abslution position of the vertices (otherwise, like D)
	virtual vec *P(bool colorcorrect=true) const;

	// returns the new position of the point pt when the mesh is warped by
	// the displacement field D (colorscalecorrected is true if D needs
	// to have gamma applied first and false if gamma has already been
	// applied -- just leave it alone).
	virtual vec warp(const vec &pt, vec *D, bool colorscalecorrected=true);

	// performs one iteration of matching the warpmesh to the target via
	// gradient descent on the energy function in Shelton's MS Thesis.
	virtual void iterate(trimesh *target) = 0;

	// Warp the whole mesh by the displacement field D.  If m is NULL, D
	// is taken to be a displacement field on the vertices of this (as
	// normal).  If m is not NULL, D is taken to be a displacement field
	// on m.
	virtual void warp(vec *D, trimesh *m = NULL,
		bool colorscalecorrected=true);

	// returns the original poisition of the vertex i.
	virtual vec origvertexpos(int i) const { return origpos[i]; }


public:
	virtual void addsprings(ilist<int> *springs) = 0;
	virtual void setdata();

	FP alpha, eta, xi, momentum, maxmove2, rho, epsilon, zeta; // parameters
	int n; // more parameters
	int maxsp; // still more parameters

	vec *origpos; // recording of the starting position

};

#endif
