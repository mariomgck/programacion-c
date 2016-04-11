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
#ifndef _LONGSET_H_
#define _LONGSET_H_

class longset {
public:
	inline longset() {
		storage = new cell[16];
		asize = 16;
		vsize = 0;
	}
	inline ~longset() {
		delete []storage;
	}

	inline void add(long i) {
		storage[hashfn(i,true)].in = true;
	}

	inline bool remove(long i) {
		int p = hashfn(i,false);
		if (p==-1) return false;
		storage[p].in = false;
		return true;
	}

	inline bool ismember(long i) {
		int p = hashfn(i,false);
		return p!=-1 && storage[p].in;
	}

	inline void clear() {
		delete []storage;
		storage = new cell[16];
		asize = 16;
		vsize = 0;
	}

private:

	class cell {
	public:
		inline cell() { used=false; }
		long i;
		bool used;
		bool in;
	};


	cell *storage;
	int asize;
	int vsize;
	
	inline unsigned long hash1(long i) {
		unsigned long n = i;
		n *= 2654435769L;
		return n%asize;
	}
	inline unsigned long hash2(long i) {
		unsigned long n = i*i;
		n *= 2654435769L;
		n %= asize;
		if ((n&1)==0) n++;
		return n;
	}

	inline int hashfn(long i, bool alloc) {
		unsigned long h1=hash1(i),h2=hash2(i);
		unsigned long long p;
		unsigned int index;
		p = h1;
		for(int j=0,p=h1;j<vsize+1;j++,p+=h2) {
			index = p%asize;
			if (!storage[index].used) {
				if (alloc) {
					if (vsize>asize/2) {
						doublesize();
						return hashfn(i,alloc);
					}
					storage[index].used = true;
					storage[index].i = i;
					vsize++;
					return index;
				} else return -1;
			} else if (i==storage[index].i) return index;
		}
		cout << "error in set hashing" << endl;
		return -1;
	}

	inline void doublesize() {
		cell *oldstorage = storage;
		storage = new cell[asize*2];
		int oldsize = asize;
		vsize = 0;
		asize *= 2;
		int p;
		for(int i=0;i<oldsize;i++) {
			if (oldstorage[i].used && oldstorage[i].in) {
				p = hashfn(oldstorage[i].i,true);
				storage[p].in = true;
			}
		}
		delete []oldstorage;
	}
};

#endif
