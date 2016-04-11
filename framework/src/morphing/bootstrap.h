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
#ifndef _BOOTSTRAP_H_
#define _BOOTSTRAP_H_

#include <geometry/model.h>
#include <warping/warper.h>

class bootstrap {
public:
	
	bootstrap(const warperparams &wp, FP cutlevel,
		trimesh *baseobj);
	bootstrap(const warperparams &wp, FP cutlevel,
		trimesh *baseobj, istream &basemodel);
	~bootstrap();

	
	inline model *getmodel() const { return m; }

	bool iterate(FP percent=0.9,int parallel=1);

	// these mesh must exist for the life of the bootstrap object
	void addobject(trimesh *tm);

	int numobjects();
	int nummodeldim();

private:
	bootstrap(const bootstrap &) { } // no can do (just yet)

	class matchprms {
	public:
		trimesh *from;
		trimesh *to;
		trimesh *ret;
		warperparams wp;
		FP cutlevel;
	};
	
	static void *match(void *voidmp);

	model *m;
	ilist<trimesh *> examples;
	warperparams p;
	trimesh *base;
	FP cl;
	
};

#endif
