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
#ifndef _GEOHASH_H_
#define _GEOHASH_H_

#include <matrix/vec.h>
#include <containers/ilist.h>
#include <globaldef.h>
#include <iostream>
#include <stdlib.h>
#include <search.h>

// this is not a simple class.  It is designed to store a collection of
// objects which can find their own axis-aligned bounding boxes and then
// find the ones which intersect an ever-growing sphere (for finding closest
// points and the like) or a growing ray (for ray tracing or the like)

// this class is passed in to either of the geometric queries -- rather
// a subclass of this class is.  process is called with each stored object
// in turn.  initaldist and process both return the maximal squared distance
// that needs to be searched (this can change each time as more objects are
// searched).  Addendum -- it is only squared distance for the method search,
// for castray and castray3d it is straight distance.
// The best way of understanding this is to look at the closestpoint method
// of trimesh

template <class T>
class geoiterate {
public:
	inline geoiterate() {}
	virtual ~geoiterate() {}

	// This function is called with a series of items of type T
	//  from a geohash.  It returns a new cutoff value for the search
	//  (a cufoff value here is the square of the maximal distance
	//  that needs to be searched).
	virtual FP process(T) = 0;
	virtual FP initaldist() = 0;
};


// This class is hardcoded pretty severly to work with 6D vectors.
// Sorry.
template <class T>
class geohash  
{
public:
	// min and max define the bounding box of the whole hash structure
	// divsize is ignored and l should be 1
	geohash(const vec &min, const vec &max, const vec &divsize, short l=1);
	// for creating a geohash before you know the bounding box (which will
	// lead to troubles if you don't call initialize before doing anything
	// else)
	inline geohash() { }
		// you must call the following function before any other if
		// you use this constructor!
	void initialize(const vec &min, const vec &max, const vec &divsize,
		short l=1);

	virtual ~geohash();

	// the type T must respond to max and min methods (which should
	// return vecs)
	void add(T t);
	void remove(T t);
	// during the move the type T can return any value for
	// max and min queries: the answers will be ignored (or rather the
	// methods never called)
	void move(T t, const vec &oldmin, const vec &oldmax,
		const vec &newmin, const vec &newmax);

	// search a ever growing sphere center center calling itt->process
	// on each found item until the value returned by itt->process is
	// less than the square of the radius of the completely searched
	// sphere (returns the smallest value returned by itt->process)
	FP search(const vec &center, geoiterate<T> &itt) const;

	// search along the ray starting at x0 and moving along v in
	// the positive direction calling itt->process until the value
	// returned is less than the shortest distance from x0 completely
	// searched (returned the smallest value
	// returned by itt->process)
	FP castray(const vec &x0, const vec &v, geoiterate<T> &itt) const;
	// same as above, but ignoring all but the first 3 axes of the vectors
	FP castray3d(const vec &x0, const vec &v, geoiterate<T> &itt) const;

	// return the maximal and minimal corners of the bounding box
	inline void getbb(vec &min, vec &max) const
		{ min = mincorner; max = maxcorner; }
	inline vec getmin() const { return mincorner; }
	inline vec getmax() const { return maxcorner; }
	// div is ignored -- don't worry about it
	inline vec getdiv() const { return (maxcorner-mincorner)/2; }

	// these need to be defined for each geohash typed created in
	// a file somewhere (see globals.cpp)

	static short globallimit;
	static int levellimit;
	static int pushlimit;
	static FP overlap;
	static int *counts;
	static int totalnum;
	static int numq;
	static int nums;

private:

	void add(T t, const vec &min, const vec &max);
	void remove(T t, const vec &min, const vec &max);
	
	class cell {
	public:
		cell() { contents = NULL; subhash = NULL; }
		~cell() { } // do not remove contents, since it may have
			// been copied to another cell (we really don't want
			// to reallocate it if we can avoid it) -- so it is
			// geohash's responsibility to take care of the
			// allocation/deallocation of contents

		ilist<T> *contents;
		geohash *subhash;
	};

	class ucfp {
	public:
		unsigned char index;
		FP value;
	};

	static int spsort(const void *e1, const void *e2) {
		const ucfp *i1 = (const ucfp *)e1;
		const ucfp *i2 = (const ucfp *)e2;
		if (i1->value < i2->value) return -1;
		if (i1->value > i2->value) return 1;
		return 0;
	}

	ilist<T> contents; // for those items which are *not* to be
		// pushed through the cells

	cell storage[64];

	short limit; // on the cell class before we move from a list to
		// a recursive hash
	short level;

	vec mincorner,maxcorner;
	
	inline unsigned char findindex(const vec &p1, const vec &p2) const {
		unsigned char index=0;
		FP divider,d1,d2;
		for(int i=0;i<DIM;i++) {
			divider = (maxcorner(i)+mincorner(i))/2;
			if (p1(i)<divider && p2(i)<divider) continue;
			if (p1(i)>=divider && p2(i)>=divider) {
				index |= (1<<i);
				continue;
			}
			d1 = overlap*(maxcorner(i)-mincorner(i))/2;
			d2 = divider+d1;
			d1 = divider-d1;
			if (p1(i)>=d1 && p2(i)<d2) {
				if (divider-p1(i) < p2(i)-divider)
					index |= (1<<i);
				// just a "random" little test (must depend
				// only on p1 and p2 though so that it
				// evaluates the same way each time (for
				// addition and removal for instance -- true
				// random would be bad) -- this puts about
				// equal numbers on each side
				continue;
			} 
			if (p2(i)<d2) continue;
			if (p1(i)>=d1) {
				index |= (1<<i);
				continue;
			}
			return 64;// can't put it into either one unfortunately
		}
		return index;
	}

	inline vec cellmin(unsigned char index) const {
		vec ret;
		for(int i=0;i<DIM;i++) {
			if (index&(1<<i))
				ret[i]=(mincorner(i)+maxcorner(i))/2-
				       overlap*(maxcorner(i)-mincorner(i))/2;
			else ret[i]=mincorner(i);
		}
		return ret;
	}
	inline vec cellmax(unsigned char index) const {
		vec ret;
		for(int i=0;i<DIM;i++) {
			if (index&(1<<i)) ret[i] = maxcorner(i);
			else ret[i]=(mincorner(i)+maxcorner(i))/2+
			             overlap*(maxcorner(i)-mincorner(i))/2;
		}
		return ret;
	}
	inline FP dist2tobox(const vec &minc, const vec &maxc,
				const vec &p) const {
		FP ret = 0;
		FP a,x,y;
		for(int i=0;i<DIM;i++) {
			a = p(i);
			x = minc(i);
			y = maxc(i);
			if (a<x && a<y) {
				a-=x;
				ret += a*a;
			} else if (a>x && a>y) {
				a-=y;
				ret += a*a;
			}
		}
		return ret;
	}

	inline FP raydisttobox(const vec &minc, const vec &maxc,
				const vec &x0, const vec &v) const {
		vec front,back;
		FP t,s;

		front = (minc-x0)/v;
		back = (maxc-x0)/v;
		for(int i=0;i<DIM;i++) {
			if (v(i)==0) {
				if (v(i)>=minc(i) && v(i)<=maxc(i)) {
					front[i] = -hugenum;
					back[i] = hugenum;
				} else return hugenum;
			} else if (front[i]>back[i]) {
				t = back[i];
				back[i] = front[i];
				front[i] = t;
			}
		}
		t = front.max();
		s = back.min();
		if (t>s) return hugenum; // ray does not intersect;
		if (s<0) return hugenum; // the intersection is in the
					// "back half" of the ray
		return t<0 ? 0 : t;
	}

	inline FP raydisttobox3d(const vec &minc, const vec &maxc,
			const vec &x0, const vec &v) const {
		vec front,back;
		FP t,s;

		front = (minc-x0)/v;
		back = (maxc-x0)/v;
		for(int i=0;i<DIM;i++) {
			if (i>2) {
				front[i] = front[0];
				back[i] = back[0];
			} else if (v(i)==0) {
				if (x0(i)>=minc(i) && x0(i)<=maxc(i)) {
					front[i] = -hugenum;
					back[i] = hugenum;
				} else return hugenum;
			} else if (front[i]>back[i]) {
				t = back[i];
				back[i] = front[i];
				front[i] = t;
			}
		}
		t = front.max();
		s = back.min();
		if (t>s) return hugenum; // ray does not intersect;
		if (s<0) return hugenum; // the intersection is in the
				// "back half" of the ray
		return t<0 ? 0 : t;
	}

	void addtocell(T t, unsigned char index,
		const vec &min, const vec &max);
	void removefromcell(T t, unsigned char index,
		const vec &min,const vec &max);
};

#include <stdlib.h>
inline FP rand48() { return (FP)rand()/(FP)RAND_MAX; }

template <class T>
geohash<T>::geohash(const vec &min, const vec &max,
			const vec &divsize, short l) {
	initialize(min,max,divsize,l);
}

template <class T>
void geohash<T>::initialize(const vec &min, const vec &max,
			const vec &divsize, short l)
{

	if (totalnum==-1) {
		totalnum = 0;
		for(int i=0;i<=levellimit;i++) counts[i] = 0;
	}
	totalnum++;
	counts[l]++;
	//if (totalnum%100==0)
	//	cout << "number of hash classes: " << totalnum << endl;
	// ignore divsize in this version:
	mincorner = min;
	maxcorner = max;
	limit = globallimit;
	level = l;
}

template <class T>
geohash<T>::~geohash()
{
	counts[level]--;
	for(int i=0;i<64;i++) {
		if (storage[i].contents != NULL) {
			delete storage[i].contents;
		}
		if (storage[i].subhash != NULL) delete storage[i].subhash;
	}
	//totalnum--;
}

#define check64(min,max,spot,cmd1,cmd2) \
{ \
	if (pushlimit==0) \
		cmd2; \
	else { \
		unsigned char mnsp=findindex(min,min), \
		              mxsp=findindex(max,max), \
		              msk=1; \
		int i,c=0; \
		for(i=0;i<DIM;i++,msk<<=1) \
			if ((mnsp&msk)!=(mxsp&msk)) c++; \
		if (c>pushlimit) \
			cmd2; \
		else { \
			unsigned char mobile=mnsp^mxsp; \
			for(spot=0;spot<1<<DIM;spot++) { \
				if (((spot^mnsp)&~mobile)==0) \
					cmd1; \
			} \
		} \
	} \
}

template <class T>
void geohash<T>::add(T t) {

	vec min = t->min();
	vec max = t->max();
	unsigned char spot = findindex(min,max);
	if (spot==64) {
		check64(min,max,spot,addtocell(t,spot,min,max),contents+=t);
	} else {
		addtocell(t,spot,min,max);
	}
}

template <class T>
void geohash<T>::add(T t, const vec &min, const vec &max) {

	unsigned char spot = findindex(min,max);
	if (spot==64) {
		check64(min,max,spot,addtocell(t,spot,min,max),contents+=t);
	} else {
		addtocell(t,spot,min,max);
	}
}

template <class T>
void geohash<T>::remove(T t) {
	vec min = t->min();
	vec max = t->max();
	unsigned char spot = findindex(min,max);
	if (spot==64)  {
		check64(min,max,spot,addtocell(t,spot,min,max),contents+=t);
	} else removefromcell(t,spot,min,max);
}

template <class T>
void geohash<T>::remove(T t,const vec &min, const vec &max) {
	unsigned char spot = findindex(min,max);
	if (spot==64) {
		check64(min,max,spot,addtocell(t,spot,min,max),contents+=t);
	} else removefromcell(t,spot,min,max);
}

// this is a complex function since we must make sure that max and min are
// not called on the object t (since it is in a state of flux)... this may
// sound trivial, but if it is sitting in the contents of a cell which then
// gets converted into a subhash, this will happen.  However, if we make sure
// to remove from a cell before adding, we should be okay.  Note that this
// also means calling the public versions of add or remove is out of the
// question.  Fortunately, the addtocell calls the private version when it
// converts a  list to a subhash.
// (I guess now with the new version of the cells it isn't so complex)
template <class T>
void geohash<T>::move(T t, const vec &oldmin, const vec &oldmax,
				   const vec &newmin, const vec &newmax) {

	unsigned char a,r;
	
	r = findindex(oldmin,oldmax);
	a = findindex(newmin,newmax);
	unsigned char amn,amx,rmn,rmx;
	unsigned char msk;
	int i;
	int ac=0,rc=0;
	if (a==64) {
		if (pushlimit==0) ac=1;
		else {
			amn=findindex(newmin,newmin);
			amx=findindex(newmax,newmax);
			for(i=0,msk=1;i<DIM;i++,msk<<=1)
				if ((amn&msk)!=(amx&msk)) ac++;
		}
	}
	if (r==64) {
		if (pushlimit==0) rc=1;
		else {
			rmn=findindex(oldmin,oldmin);
			rmx=findindex(oldmax,oldmax);
			for(i=0,msk=1;i<DIM;i++,msk<<=1)
				if ((rmn&msk)!=(rmx&msk)) rc++;
		}
	}
	if (r==64 && rc>pushlimit) {
		if (a==64 && ac>pushlimit) return;
		contents.remove(t);
		if (a==64) add(t,newmin,newmax);
		else addtocell(t,a,newmin,newmax);
	} else if (a==64 && ac>pushlimit) {
		if (r==64) remove(t,oldmin,oldmax);
		else removefromcell(t,r,oldmin,oldmax);
		contents += t;
	} else if (a==64) {
		if (r==64) {
			unsigned char spot;
			unsigned char amobile = amn^amx;
			unsigned char rmobile = rmn^rmx;
			for(spot=0;spot<1<<DIM;spot++) {
				if (((spot^amn)&~amobile)==0) {
				   if (((spot^rmn)&~rmobile)==0) {
				      if (storage[spot].subhash!=NULL)
				         storage[spot].subhash->
				           move(t,oldmin,oldmax,newmin,newmax);
				   } else {
				      addtocell(t,spot,newmin,newmax);
				   }
				} else {
				   if (((spot^rmn)&~rmobile)==0) {
				      removefromcell(t,spot,oldmin,oldmax);
				   }
				}
			}
		} else {
			unsigned char spot;
			unsigned char amobile = amn^amx;
			for(spot=0;spot<1<<DIM;spot++) {
				if (((spot^amn)&~amobile)==0) {
				   if (spot==r) {
				      if (storage[spot].subhash!=NULL)
				         storage[spot].subhash->
				           move(t,oldmin,oldmax,newmin,newmax);
				   } else {
				      addtocell(t,spot,newmin,newmax);
				   }
				} else {
				   if (spot==r) {
				      removefromcell(t,spot,oldmin,oldmax);
				   }
				}
			}
		}
	} else if (r==64) {
		unsigned char spot;
		unsigned char rmobile = rmn^rmx;
		for(spot=0;spot<1<<DIM;spot++) {
			if (spot==a) {
			   if (((spot^rmn)&~rmobile)==0) {
			      if (storage[spot].subhash!=NULL)
				 storage[spot].subhash->
				   move(t,oldmin,oldmax,newmin,newmax);
			   } else {
			      addtocell(t,spot,newmin,newmax);
			   }
			} else {
			   if (((spot^rmn)&~rmobile)==0) {
			      removefromcell(t,spot,oldmin,oldmax);
			   }
			}
		}
	} else {
		if (a==r) {
			if (storage[a].subhash!=NULL)
				storage[a].subhash->move(t,oldmin,oldmax,
				                             newmin,newmax);
		} else {
			removefromcell(t,r,oldmin,oldmax);
			addtocell(t,a,newmin,newmax);
		}
	}
}


template <class T>
void geohash<T>::addtocell(T t, unsigned char index,
		const vec &min,const vec &max) {

	if (storage[index].contents!=NULL) {
		// if the list just overflowed, make the sub geohash
		// and move everything
		if (storage[index].contents->length()==limit &&
		    level<levellimit) {
			storage[index].subhash =
				new geohash<T>(cellmin(index),
				               cellmax(index),0,level+1);
			storage[index].subhash->add(t,min,max);
			for(int i=0;i<storage[index].contents->length();i++)
				storage[index].subhash->
					add((*storage[index].contents)(i));
			delete storage[index].contents;
			storage[index].contents = NULL;
		} else {
			(*storage[index].contents) += t;
		}
	// if we have a sub geohash, add it to that
	} else if (storage[index].subhash!=NULL) {
		storage[index].subhash->add(t,min,max);
	} else {
		storage[index].contents = new ilist<T>;
		(*storage[index].contents) += t;
	}
}

template <class T>
void geohash<T>::removefromcell(T t, unsigned char index,const vec &min,const vec &max) {

	if (storage[index].contents!=NULL) {
		storage[index].contents->remove(t);
	} else {
		storage[index].subhash->remove(t,min,max);
	}
}

template <class T>
FP geohash<T>::search(const vec &center, geoiterate<T> &itt) const {

	FP currr = dist2tobox(mincorner,maxcorner,center);
	FP bestr = itt.initaldist();
	FP tr;
	
	if (level==1) nums++;
	int mi = contents.length(),i;
	for(i=0;i<mi && bestr>currr;i++) {
		numq++;
		bestr = itt.process(contents(i));
	}
	if (currr>=bestr) {
		//if (level==1) cout << endl;
		return bestr;
	}
	ucfp indices[1<<DIM];
	int c;
	for(i=0,c=0;i<1<<DIM;i++) {
		if (storage[i].contents==NULL && storage[i].subhash==NULL)
			continue;
		indices[c].index = i;
		indices[c].value = dist2tobox(cellmin(i),cellmax(i),center);
		c++;
	}
	qsort(indices,c,sizeof(ucfp),geohash<T>::spsort);
	for(i=0,currr=indices[0].value;    i<c && bestr>currr;
				currr=indices[++i].value) {
		mi = indices[i].index;
		if (storage[mi].contents) {
			int mj=storage[mi].contents->length();
			int bi = 0;
			for(int j=0;j<mj&&bestr>currr;j++) {
				numq++;
				tr =
				   itt.process(storage[mi].contents->nth(j));
				if (tr<bestr) {
					bestr = tr;
					bi = j;
				}
			}
			if (bi!=0) {
				T temp = storage[mi].contents->nth(bi-1);
				storage[mi].contents->
				   setnth(storage[mi].contents->nth(bi),bi-1);
				storage[mi].contents->setnth(temp,bi);
			}
		} else if (storage[mi].subhash) {
			tr=storage[mi].subhash->search(center,itt);
			if (tr<bestr) bestr=tr;
		}
	}
	return bestr;
}

template <class T>
FP geohash<T>::castray(const vec &x0, const vec &v, geoiterate<T> &itt) const {

	FP currt = raydisttobox(mincorner,maxcorner,x0,v);
	if (currt==hugenum) return itt.initaldist();
	FP bestt = itt.initaldist();
	FP tt;

	int mi = contents.length(),i;
	for(i=0;i<mi && bestt>currt;i++)
		bestt = itt.process(contents(i));
	if (currt>=bestt) {
		return bestt;
	}
	ucfp indices[1<<DIM];
	for(i=0;i<1<<DIM;i++) {
		indices[i].index = i;
		indices[i].value = raydisttobox(cellmin(i),cellmax(i),x0,v);
	}
	qsort(indices,1<<DIM,sizeof(ucfp),geohash<T>::spsort);
	for(i=0,currt=indices[0].value;i<1<<DIM && bestt>currt &&
				currt!=hugenum;currt=indices[++i].value) {
		mi = indices[i].index;
		if (storage[mi].contents) {
			int mj=storage[mi].contents->length();
			for(int j=0;j<mj&bestt>currt;j++)
				bestt =
				   itt.process(storage[mi].contents->nth(j));
		} else if (storage[mi].subhash) {
			tt = storage[mi].subhash->castray(x0,v,itt);
			if (tt<bestt) bestt = tt;
		}
	}
	return bestt;
}

template <class T>
FP geohash<T>::castray3d(const vec &x0, const vec &v, geoiterate<T> &itt) const {

	FP currt = raydisttobox3d(mincorner,maxcorner,x0,v);
	if (currt==hugenum) return itt.initaldist();
	FP bestt = itt.initaldist();
	FP tt;

	int mi = contents.length(),i;
	for(i=0;i<mi && bestt>currt;i++)
		bestt = itt.process(contents(i));
	if (currt>=bestt) {
		return bestt;
	}
	ucfp indices[1<<DIM];
	for(i=0;i<1<<DIM;i++) {
		indices[i].index = i;
		indices[i].value = raydisttobox3d(cellmin(i),cellmax(i),x0,v);
	}
	qsort(indices,1<<DIM,sizeof(ucfp),geohash<T>::spsort);
	for(i=0,currt=indices[0].value;i<1<<DIM && bestt>currt &&
				currt!=hugenum;currt=indices[++i].value) {
		mi = indices[i].index;
		if (storage[mi].contents) {
			int mj=storage[mi].contents->length();
			for(int j=0;j<mj&bestt>currt;j++)
				bestt =
				   itt.process(storage[mi].contents->nth(j));
		} else if (storage[mi].subhash) {
			tt = storage[mi].subhash->castray3d(x0,v,itt);
			if (tt<bestt) bestt = tt;
		}
	}
	return bestt;
}
		
#endif
