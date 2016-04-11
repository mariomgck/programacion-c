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
#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <globaldef.h>
#include <containers/ilist.h>
#include <math.h>
#include "float.h"

// a general image class (all image pixels are stored at floats/doubles)
// Note that copying an image (i.e. image im1 = im2) yields two images
// that share the same memory (thus changing im1 will change im2).  However,
// there are reference coutns on each image channel and so you can delete im1
// or im2 first and it won't matter.  Use copy to get a deep copy of the image
class image  
{
public:
	inline image(int w=-1, int h=-1, bool xw=false, bool yw=false)
		{ width = w; height = h; xwrap=xw; ywrap=yw; }
	inline image(int w, int h, const ilist<FP *> &ch, bool xw=false, bool yw=false) {
		width = w; height = h;
		for(int i=0;i<ch.length();i++)
			chs += new chptr(ch(i),false);
		xwrap = xw; ywrap = yw;
	}
	inline image(int w, int h, FP **ch, int numch, bool xw=false, bool yw=false) {
		width = w; height = h; for(int i=0;i<numch;i++) chs += new chptr(ch[i],false);
		xwrap = xw; ywrap = yw;
	}
	inline image(const image &im) {
		width = im.width; height = im.height; xwrap = im.xwrap; ywrap = im.ywrap;
		for(int i=0;i<im.chs.length();i++)
			chs += new chptr(*im.chs(i));
	}
	inline image &operator=(const image &im) {
		for(int i=0;i<chs.length();i++)
			delete chs[i];
		chs.setlength(0);
		width = im.width; height = im.height;
		xwrap = im.xwrap; ywrap = im.ywrap;
		for(int i=0;i<im.chs.length();i++)
			chs += new chptr(*im.chs(i));
		return *this;
	}

	// this returns true only if the image channel point to the same memory (it
	// does not check actual content)
	inline bool operator==(const image &im) {
		if (im.width!=width || im.height!=height || im.xwrap != xwrap ||
			im.ywrap!=ywrap || im.chs.length() != chs.length()) return false;
		for(int i=0;i<chs.length();i++)
			if (chs(i)->cptr!=im.chs(i)->cptr) return false;
		return true;
	}

	// deep copy only explicitly
	image copy() const;

	// remove all channels and reset the size and wraps
	inline void reset(int neww, int newh, bool newxw=false, bool newyw=false) {
		for(int i=0;i<chs.length();i++) delete chs[i];
		width = neww;
		height = newh;
		xwrap = newxw;
		ywrap = newyw;
	}

	inline ~image() { for(int i=0;i<chs.length();i++) delete chs[i]; }

	inline int getwidth() const { return width; }
	inline int getheight() const { return height; }
	// return whether the image wraps in a particular direction
	inline bool wrapx() const { return xwrap; }
	inline bool wrapy() const { return ywrap; }


	inline void addchannel() { chs += new chptr(new FP[width*height],true); }

	inline int numchannels() const { return chs.length(); }

	// returns a pointer (which the image still owns!) to a particular channel
	inline FP *channel(int c) const { return chs(c)->cptr->c; }

	inline FP value(int x, int y, int c=0) const { return channel(c)[y*width+x]; }
	inline FP v(int i, int c=0) const { return channel(c)[i]; }

	inline FP &pixel(int x, int y, int c=0) { return channel(c)[y*width+x]; }
	inline FP &p(int i, int c=0) { return channel(c)[i]; }

	// returns a vector of the channel values at a point (in ret)
	inline void getvec(int x, int y, FP *&ret) const { 
		int index = y*width+x;
		for(int i=0;i<chs.length();i++)
			ret[i] = channel(i)[index];
	}
	inline void getvec(int index, FP *&ret) const {
		for(int i=0;i<chs.length();i++)
			ret[i] = channel(i)[index];
	}

	// set all values for the image (or one channel) to the same number
	inline void set(const FP &v) {
		for(int c=0;c<chs.length();c++)
			set(v,c);
	}
	inline void set(const FP &v, int c) {
		int m = width*height;
		for(int i=0;i<m;i++) channel(c)[i] = v;
	}

	// add a constant
	inline void add(const FP &v) {
		for(int c=0;c<chs.length();c++)
			add(v,c);
	}
	inline void add(const FP &v, int c) {
		int m = width*height;
		for(int i=0;i<m;i++) channel(c)[i] += v;
	}

	// subtract a constant
	inline void sub(const FP &v) {
		for(int c=0;c<chs.length();c++)
			sub(v,c);
	}
	inline void sub(const FP &v, int c) {
		int m = width*height;
		for(int i=0;i<m;i++) channel(c)[i] -= v;
	}

	// multiply by a constant
	inline void mult(const FP &v) {
		for(int c=0;c<chs.length();c++)
			mult(v,c);
	}
	inline void mult(const FP &v, int c) {
		int m = width*height;
		for(int i=0;i<m;i++) channel(c)[i] *= v;
	}

	// for historical reasons: this does the same as channel
	inline FP *getchannel(int c) const { return channel(c); }

	// add a image to this one (on a pixel by pixel basis)
	inline bool add(const image &im) {
		if (im.chs.length()!=chs.length() ||
			im.width != width || im.height != height) return false;
		int size = width*height;
		for(int c=0;c<chs.length();c++)
			for(int i=0;i<size;i++)
				channel(c)[i] += im.channel(c)[i];
		return true;
	}

	// subtract an image from this one
	inline bool sub(const image &im) {
		if (im.chs.length()!=chs.length() ||
			im.width != width || im.height != height) return false;
		int size = width*height;
		for(int c=0;c<chs.length();c++)
			for(int i=0;i<size;i++)
				channel(c)[i] -= im.channel(c)[i];
		return true;
	}

	// multiply an image by this one
	inline bool mult(const image &im) {
		if (im.chs.length()!=chs.length() ||
			im.width != width || im.height != height) return false;
		int size = width*height;
		for(int c=0;c<chs.length();c++)
			for(int i=0;i<size;i++)
				channel(c)[i] *= im.channel(c)[i];
		return true;
	}

	// set this image to be the same as another (shallow copy)
	inline bool set(const image &im) {
		if (im.chs.length()!=chs.length() ||
			im.width != width || im.height != height) return false;
		int size = width*height;
		for(int c=0;c<chs.length();c++)
			for(int i=0;i<size;i++)
				channel(c)[i] = im.channel(c)[i];
		return true;
	}

	// reduce to one channel (being the sum of all of the others)
	inline void flatten() {
		int size=height*width;
		for(int i=0;i<size;i++)
			for(int c=1;c<chs.length();c++)
				channel(c)[0] += channel(c)[i];
		for(int i=1;i<chs.length();i++)
			delete chs[i];
		chs.setlength(1);
	}

	// blur by exanding and reducing using the taps passed in
	inline image blur(int ntaps, FP *taps, int level=1) const {
		ilist<int> w,h;
		image ret=*this;
		for(int i=0;i<level;i++) {
			w += ret.width;
			h += ret.height;
			ret = ret.reduce(ntaps,taps);
		}
		for(int i=level-1;i>=0;i--)
			ret = ret.expand(ntaps,taps,w(i),h(i));
		return ret;
	}

	// find the max and min over all channels
	inline void limit(FP min, FP max) {
		for(int i=0;i<chs.length();i++)
			limit(min,max,i);
	}
	// find the max and min over one channel
	inline void limit(FP min, FP max, int c) {
		int size=width*height;
		for(int i=0;i<size;i++)
			if (channel(c)[i]<min) channel(c)[i]=min;
			else if (channel(c)[i]>max) channel(c)[i]=max;
	}

	// reduce using the taps passed in (the new returned image is 1/2 the width and 1/2 the height)
	image reduce(int ntaps, FP *taps) const;
	// expand using the taps passed in (the new returned image is twice the width and twice the height --
	//  although these can be changed by setting w and h -- they should only be changed by
	//  one pixel one way or the other)
	image expand(int ntaps, FP *taps, int w=-1, int h=-1) const;
	image backwarp(const image &backflow, FP scale=1) const;
	// return an image of the x axis derivative of this image on a per channel basis
	image dx() const;
	// ditto for the y derivative
	image dy() const;
	// on a per channel basis, compute the sum of the sqaures of the elements inside a size by size
	// window centered on a given pixel and place the value in that pixel.  Returns a new image
	// and assumes that size<width and size<height
	image squaresum(int size, bool onlyfinite=true) const;

	// saves the image as a ppm file with channels c1=red c2=green c3=blue.  min and max
	// are used to scale the values into the needed [0,1] range (-hugenum or hugenum implies
	// the value should be set to the global min or max -- global across all channels, not
	// just c1, c2, and c3)
	void save(const char * const name, int c1, int c2, int c3, FP min=-hugenum, FP max=hugenum);
	// load an image from a ppm file (the range is 0 to 1)
	bool load(const char * const name); // returns success

	// return the max (onlyfinite ignores values that aren't finite)
	FP max(bool onlyfinite=true) {
		FP ret = -hugenum,t;
		for(int i=0;i<chs.length();i++)
			if ((t=max(i,onlyfinite))>ret) ret = t;
		return ret;
	}
	// max of a channel
	FP max(int c,bool onlyfinite=true) {
		FP ret = -hugenum;
		int size=height*width;
		for(int i=0;i<size;i++)
			if (channel(c)[i] > ret && (!onlyfinite || _finite(channel(c)[i]))) ret = channel(c)[i];
		return ret;
	}
	// same with mins:
	FP min(bool onlyfinite=true) {
		FP ret = hugenum,t;
		for(int i=0;i<chs.length();i++)
			if ((t=min(i,onlyfinite))<ret) ret = t;
		return ret;
	}
	FP min(int c,bool onlyfinite=true) {
		FP ret = hugenum;
		int size=height*width;
		for(int i=0;i<size;i++)
			if (channel(c)[i] < ret && (!onlyfinite || _finite(channel(c)[i]))) ret = channel(c)[i];
		return ret;
	}

private:
	// the next two private classes are for doing garbage collection and memory management
	// for the image channels (this allows a shallow copy)
	class chtype {
	public:
		inline chtype(FP *ch, bool own) { mine=own; ref=0; c = ch; }
		inline ~chtype() { if (mine) { ref=0; delete c; } }

		FP *c;
		int ref;
		bool mine;
	};
	class chptr {
	public:
		inline chptr() { cptr=NULL; }
		inline chptr(const chptr &cp) { cptr = cp.cptr; if (cptr!=NULL) cptr->ref++; }
		inline chptr(chtype *cp) { cptr=cp; if (cptr!=NULL) cptr->ref++; }

		inline ~chptr() { if (cptr==NULL) return;
					cptr->ref--; if (cptr->ref==0) delete cptr; cptr = NULL; }
		
		inline chptr(FP *c,bool own) { cptr = new chtype(c,own); cptr->ref++; }

		inline chptr &operator=(const chptr& cp) {
			if (cptr!=NULL) {
				cptr->ref--;
				if (cptr->ref==0) delete cptr;
			}
			cptr = cp.cptr;
			if (cptr!=NULL)
				cptr->ref++;
			return *this;
		}
		chtype *cptr;
	};

	// these have to be pointer to chptr's since the ilist object does not use new,
	//  but rather malloc/realloc to allocate memory.
	ilist<chptr *> chs;
	int width,height;
	bool xwrap,ywrap;

	FP *filterrows(FP *ch, int w, int h, int ntaps, FP *taps) const;
	FP *filtercolumns(FP *ch, int w, int h, int ntaps, FP *taps) const;
	FP *halfrows(FP *ch, int w, int h) const;
	FP *halfcolumns(FP *ch, int w, int h) const;
	FP *doublerows(FP *ch, int w, int h, int newh) const;
	FP *doublecolumns(FP *ch, int w, int h, int neww) const;
};

#endif
