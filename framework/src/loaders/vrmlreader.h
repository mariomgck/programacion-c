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
#ifndef _VRMLREADER_H_
#define _VRMLREADER_H_

#include <globaldef.h>
#include <geometry/trimesh.h>

class vrmlreader {
public:
	vrmlreader(const char *filename, bool hashon=true);
	virtual ~vrmlreader();

	inline trimesh *mesh() { return ret; }

private:
	trimesh *ret;
};

#endif
