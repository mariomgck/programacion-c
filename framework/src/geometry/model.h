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
#ifndef _MODEL_H_
#define _MODEL_H_

#include <geometry/trimesh.h>
#include <warping/warperparams.h>
#include <fstream>
#include <matrix/gmatrix.h>

using namespace std;

class model {
public:
	model(const trimesh &basemesh, const warperparams &p,
		FP cutlevel=0.115,
		ilist<vec> *pts = NULL, bool split=false);
	model(const ilist<char *> &basefiles, const warperparams &p,
		const char * pts = NULL, bool split=false);
	model(istream &is, bool split=false);

	~model();

	// ref is whether the meshes should be saved as references
	//     to files or not (if not, they are saved completely inside
	//     of the stream)
	// complete is whether the full set of meshes and the parameters
	//     should be saved (so that it would be possible to add new
	//     objects later) or not (in which case, only the full version
	//     of the base object is saved).
	// notice that depending on how the model was created, some
	// combinations of these parameters may not make sense.  In that case,
	// the function will return false.  Otherwise, the function will
	// return true.
	bool save(ostream &os, bool ref=false, bool complete=false);

	bool iscomplete() { return !incomplete; }

	void changeparams(const warperparams &p);

	void addobject(const ilist<char *> &filelist,
		const char * pts = NULL,
		bool split=false, int vernum=0);
	void addobject(trimesh *m, FP cutlevel=0.115, ilist<vec> *pts = NULL,
		bool split=false, int vernum=0);
		
	void addvector(const ilist<vec> &v, bool split=false);
	void addvector(vec *v, bool split=false);

	int numdim();
	int veclength();

	void setparameter(int i, FP v);
	void setparameters(const ilist<FP> &vs);
	FP getparameter(int i);
	void getparameters(ilist<FP> &vs, bool append=false);

	void getvector(int i, ilist<vec> &v, bool append=false);
	void getvector(int i, vec *v);

	void getorigin(ilist<vec> &v, bool append=false);
	void getorigin(vec *v);

	trimesh *getmesh(); // the returned value will change as the
		// parameters are adjusted

	void matchshape(trimesh *shape, FP pchange=0.01,
		int nitt=20, int npts=200, FP lambda=0,
		void (*cbfn)(void *,int,int) = NULL,void *cbdata = NULL);

	// assumes that the dimensions have been orthongonalized by
	// the second method!
	void matchparms(vec *v);

	void zeroorigin();
	void setcurrenttoorigin(bool changebasis=false);
	void orthogonalize();
	void orthogonalize2(FP **sret);

	void split(int v=-1);
	void removedim(int d);
	void removelast(int d);

private:

	vec splitvec(const vec &v, int sn) {
		int i;
		vec ret(v);
		if (sn==1) for(i=3;i<DIM;i++) ret[i] = 0;
		else for(i=0;i<3;i++) ret[i] = 0;
		return ret;
	}
	void warpbase();

	void setbasevecs(bool split);
	void completelist();
	void setpmatrix(gmatrix &m);
	void pfrommatrix(const gmatrix &m,FP *scales=NULL);
	
	ilist<char *> files;
	warperparams params;
	trimesh *base;
	ilist<trimesh *> m;
	ilist<vec *> p;
	vec *origin;
	bool incomplete;
	ilist<vec> upts;
	char *uptsname;
	ilist<FP> alpha;
	FP basecutlevel;
};

#endif
	

	
