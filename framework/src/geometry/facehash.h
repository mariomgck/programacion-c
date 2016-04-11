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
#ifndef _FACEHASH_H_
#define _FACEHASH_H_


#include <globaldef.h>
#include <iostream>
#include <stdlib.h>

#include <geometry/face.h>
#include <matrix/vec.h>

#include <containers/ilist.h>
#include <containers/pqueue.h>
#include <containers/longset.h>

class geoiterateface {
public:
	inline geoiterateface() {}
	virtual ~geoiterateface() {}

	virtual FP process(face *) = 0;
	virtual FP initaldist() = 0;
};


class facehash  
{
public:
	// min and max define the bounding box of the whole hash structure
	facehash(vec &min, vec &max, vec &divsize, int l=0);
	facehash(int l=0);
	void initialize(vec &min, vec &max, vec &divsize);

	virtual ~facehash();

	void add(face *t);
	void remove(face *t);

	// search a ever growing sphere center center calling itt->process
	// on each found item until the value returned by itt->process is
	// less than the square of the radius of the completely searched
	// sphere (returns the smallest value returned by itt->process)
	float search(vec &center, geoiterateface &itt);

	// search along the ray starting at x0 and moving along v in
	// the positive direction calling itt->process until the value
	// returned is less than the shortest distance from x0 completely
	// searched (returned the smallest value
	// returned by itt->process)
	float castray(vec &x0, vec &v, geoiterateface &itt) ;
	// same as above, but ignoring all but the first 3 axes of the vectors
	float castray3d(vec &x0, vec &v, geoiterateface &itt) ;

	// return the maximal and minimal corners of the bounding box
	inline void getbb(vec &min, vec &max) const
		{ min = mincorner; max = maxcorner; }
	inline vec getmin() { return mincorner; }
	inline vec getmax() { return maxcorner; }
	// div is ignored -- don't worry about it
	inline vec getdiv() { return div; }

	static int ncells[DIM];
	static int nums,numq;
	static int maxlevel;
	static int maxcontents;

	inline void printstats() {
		if (vsize==0) return;
		if (numq%1000 != 0) return;
		cout << "queries: " << numq << '/' << nums << endl;
		cout << "cell num: " << vsize << '/' << asize << endl;
	}

private:

	vec mincorner, maxcorner, div;
	ilist<face *> toadd;
	longset toaddset;
	longset toremoveset;

	int level;

	void addall();

	class cell {
	public:
		inline cell() { contents = NULL; subhash = NULL; }
		short i[DIM];
		ilist<face *> *contents;
		facehash *subhash;
	};

	class qcell {
	public:
                short i[DIM];
                FP celldist;
                inline FP rank() { return -celldist; }
                inline void setposition(int i) { }
        };

	class qcellptr {
	public:
		inline qcellptr() { p = NULL; }
		inline qcellptr(const qcellptr &qc) { p = qc.p; } 
		inline qcellptr(qcell *ptr) { p = ptr; }
		inline FP rank() { return -p->celldist; }
		inline void setposition(int i) { }

		qcell *p;
	};

	qcell *tocell(vec &c, longset &ls) {
		qcell *r = new qcell;
		FP t;
		r->celldist = 0;
		long uniqueid = 0;
		for(int i=0;i<DIM;i++) {
			r->i[i] = (short)floor((c[i]-mincorner[i])/div[i]);
			if (r->i[i]<0) r->i[i] = 0;
			if (r->i[i]>=ncells[i]) r->i[i] = ncells[i]-1;
			if (c[i] < (t=r->i[i]*div[i]+mincorner[i])) {
				t -= c[i];
				r->celldist += t*t;
			} else if (c[i] > (t+=div[i])) {
				t = c[i]-t;
				r->celldist += t*t;
			}
			uniqueid = uniqueid*ncells[i] + r->i[i];
		}
		ls.add(uniqueid);
		return r;
	}

	qcell *tocell(vec &c, short i[DIM], short d, short add,
		longset &ls) {
		if (add<0 && i[d]==0) return NULL;
		if (add>0 && i[d]==ncells[d]-1) return NULL;
		qcell *r = new qcell;
		r->celldist = 0;
		FP t;
		long uniqueid = 0;
		for(int k=0;k<DIM;k++) {
			r->i[k] = k==d ? i[k]+add : i[k];
			uniqueid = uniqueid*ncells[k] + r->i[k];
			if (c[k] < (t=r->i[k]*div[k]+mincorner[k])) {
				t -= c[k];
				r->celldist += t*t;
			} else if (c[k] > (t+=div[k])) {
				t = c[k]-t;
				r->celldist += t*t;
			}
		}
		if (ls.ismember(uniqueid)) {
			delete r;
			return NULL;
		}
		ls.add(uniqueid);
		return r;
	}
		
	cell *storage;
	int vsize;
	int asize;

	inline bool equal(short i[DIM],short j[DIM]) {
		for(int k=0;k<DIM;k++) if (i[k]!=j[k]) return false;
		return true;
	}
	
	inline unsigned long hash1(short i[DIM]) {
		unsigned long n = 0;
		for(int k=0;k<DIM;k++)
			n = n*ncells[k]+i[k];
		n *= 2654435769L;
		return n%asize;
	}
	inline unsigned long hash2(short i[DIM]) {
		unsigned long n=0;
		for(int k=DIM-1;k>=0;k--)
			n = n*ncells[k]+i[k];
		n *= 2654435769L;
		n %= asize;
		if ((n&1) == 0) n++;
		return n;
	}

	inline unsigned int hashfn(short i[DIM], bool alloc) {
		unsigned long h1=hash1(i),h2=hash2(i);
		unsigned long long p;
		unsigned int index;
		p = h1;
		for(int j=0,p=h1;j<vsize+1;j++,p+=h2) {
			index = p%asize;
			if (!storage[index].contents &&
			    !storage[index].subhash) {
				if (alloc) {
					if (vsize>asize/2) {
						doublesize();
						return hashfn(i,alloc);
					}
					for(int k=0;k<DIM;k++)
						storage[index].i[k] = i[k];
					vsize++;
					return index;
				} else return -1;
			} else if (equal(i,storage[index].i)) return index;
		}
		cout << h1 << ' ' << h2 << endl;
		cout << vsize << ' ' << asize << endl;
		cerr << "hashing function problem!" << endl;
		return -1;
	}

	inline void doublesize() {
		//cout << "doubling (" << level << "): " << vsize << ' ' << asize << endl;
		cell *oldstorage = storage;
		storage = new cell[asize*2];
		int oldsize = asize;
		vsize = 0;
		asize *= 2;
		unsigned int p;
		for(int i=0;i<oldsize;i++) {
			if (oldstorage[i].contents &&
			    oldstorage[i].contents->length() > 0) {
				p = hashfn(oldstorage[i].i,true);
				storage[p].contents = oldstorage[i].contents;
			} else if (oldstorage[i].contents) {
				delete oldstorage[i].contents;
			} else if (oldstorage[i].subhash) {
				p = hashfn(oldstorage[i].i,true);
				storage[p].subhash = oldstorage[i].subhash;
			}
		}
		//cout << "done doubling: " << vsize << ' ' << asize << endl;
		delete []oldstorage;
	}

	void dorange(face *f, void (*fn)(short *,face *));
	void addcell(short i[DIM], face *f);
	void removecell(short i[DIM], face *f);
};

#endif
