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
#ifndef _PPMREADER_H_
#define _PPMREADER_H_

#include <iostream>
#include "globaldef.h"
#include <image/image.h>

// clearly this isn't a good reader for some apps since it stores things
// in floating point -- however, this is what I need now and it can be changed
// later if needed
// (clearly I might also want << and >> operators and the like later!)

// a class that reads in a ppm file.. Just create one and pass in the ppm file
// stream (which *must* be opened with ios::binary as the second argument)
// You can then access the inidividual pixels or the channels.
class ppmreader  
{
public:
	// in this case the ppmreader does not own the channel memory
	ppmreader(istream &s,FP *r, FP *g, FP *b, FP mult=1.0);
	// in this case the ppmreader does own the channel memory and it will be
	//  deleted when the reader is destroyed (unless release is called first)
	ppmreader(istream &s, FP mult=1.0);
	// again, the ppmreader owns the channel memory
	ppmreader(int width, int height, double initr=0, double initg=0, double initb=0);
	// and here it doesn't
	ppmreader(int width, int height, FP *r, FP *g, FP *b);

	// this makes one from a 1 or 3 channel image:
	ppmreader(const image &im);

	~ppmreader();

	inline int height() const { return h; }
	inline int width() const { return w; }

	// get references to the pixel values
	inline FP &red(int x, int y) { return r[y*w+x]; }
	inline FP &green(int x, int y) { return g[y*w+x]; }
	inline FP &blue(int x, int y) { return b[y*w+x]; }
	inline FP &red(int i) { return r[i]; }
	inline FP &green(int i) { return g[i]; }
	inline FP &blue(int i) { return b[i]; }

	// get the pixel values
	inline FP getr(int x,int y) const { return r[y*w+x]; }
	inline FP getg(int x, int y) const { return g[y*w+x]; }
	inline FP getb(int x, int y) const { return b[y*w+x]; }
	inline FP getr(int i) const { return r[i]; }
	inline FP getg(int i) const { return g[i]; }
	inline FP getb(int i) const { return b[i]; }

	// get the channels
	inline FP *redchannel() const { return r; }
	inline FP *greenchannel() const { return g; }
	inline FP *bluechannel() const { return b; }

	// tell the ppmreader is no longer owns the channel memory
	inline void release() { owndata = false; }

	// save the channels to a ppm file (which must be opened with ios::binary
	// permissions if ascii==false
	void save(ostream &s, int range=255, bool ascii=true, bool color=true,
		FP min=0, FP max=1);

private:
	ppmreader(const ppmreader &) { }; // don't want this to happen!
	void loadit(istream &s, FP mult);

	bool owndata;
	FP *r,*g,*b;
	int w,h;
};

#endif
