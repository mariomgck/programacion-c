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
#ifndef _WARPER_H_
#define _WARPER_H_

#include <containers/ilist.h>
#include <matrix/vec.h>
#include <warping/awarpmesh.h>
#include <warping/warperparams.h>

// the warper class control the matching of two "pyramids" of meshes
class warper  
{
public:
	// create a warper class ready to match the pyramids listed in from
	// and to.  The two lists need not be the same length.  The warp will,
	// in the end, match the first of the from list to the first of the
	// to list.  The lists are listed in the order from most complex to
	// least complex and should be the names of files saved by a trimesh
	// class (with type trimesh::tmff1)
	warper(const ilist<char *> &from, const ilist<char *> &to,
	       const char * const frompts, const char * const topts,
		   const warperparams & params,int vernum=0);
	warper(const ilist<trimesh *> *from, const ilist<char *> &to,
		const ilist<vec> &frompts, const char * const topts,
		const warperparams &params,int vernum=0);
	warper(const ilist<trimesh *> *from, const ilist<trimesh *> *to,
		const ilist<vec> &frompts, const ilist<vec> &topts,
		const warperparams &params,int vernum=0);
	warper(const ilist<char *> &from, const ilist<trimesh *> *to,
		const char * const frompts, const ilist<vec> &topts,
		const warperparams &params,int vernum=0);
	virtual ~warper();

	// take a step in the process
	bool iterate(); // returns true if not yet done and false if done

	// do the whole process
	inline void run() { while(iterate()); }

	// returns the displacement vector for all of the vertices of the
	// current A.  Just leave colorcorrect to true (this makes sure that
	// the vector has been corrected for the gamma setting properly and
	// the output is the same regardless of gamma).  size gets the length
	// of the returned array:  the callee now owns the array's memory.
	inline vec *D(int &size, bool colorcorrect=true) const
		{ size = A->numvertices();
		  return addaffine(A->D(colorcorrect),size); }

	// returns the absolution position of the vertices of the current A
	// (in a similar format to the D of above)
	inline vec *P(int &size, bool colorcorrect=true) const
		{ size = A->numvertices(); return A->P(colorcorrect); }

	// reads in a file (which is just a list of the vectors returned by
	// the above) and warps the current mesh accordingly.
	inline void applywarp(istream &s) {
		vec *D = new vec[A->numvertices()];
		for(int i=0;i<A->numvertices();i++)
			s >> D[i];
		A->warp(D);
		delete []D;
	}

	// snap all the points on A to their closest points on B
	// (often done after the iterations are complete.).
	void snap();

	inline awarpmesh *getA(void) { return A; }
	inline trimesh *getB(void) { return B; }

private:

	bool nextlevel();
	void loadpts(const char * const frompts, const char * const topts);
	void loadpts(const ilist<vec> &frompts, const char * const topts);
	void loadpts(const ilist<vec> &frompts, const ilist<vec> &topts);
	void loadpts(const char * const frompts, const ilist<vec> &topts);
	void completeload();
	void findaffine();

	vec *addaffine(vec *x, int s) const {
		matrix Ainv = affineA.inv();
		for(int i=0;i<s;i++)
			x[i] += A->origvertexpos(i) -
				Ainv*(A->origvertexpos(i)-affineb);
		return x;
	}

	matrix affineA;
	vec affineb;

	awarpmesh *A;
	trimesh *Aold;
	trimesh *B;

	warperparams p;
	int ittnum;
	int level;

	ilist<char *> As;
	ilist<char *> Bs;
	ilist<trimesh *> *Ams;
	ilist<trimesh *> *Bms;

	ilist<vec> tovec;
	ilist<vec> fromvec;
	
	int warpv;
};

#endif
