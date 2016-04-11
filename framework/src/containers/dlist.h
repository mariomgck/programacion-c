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
#ifndef _DLIST_H_
#define _DLIST_H_

#include <stdlib.h>
#include <iostream>
#include "malloc.h"

template <class T>
class dlist {

	public:

		dlist (int initsize = 10);
		dlist (T *s, int size, int realsize);
		~dlist();
		T nth(int n) const { return storage[n]; }
		
		inline void setnth(const T &item, int n) {
			if (n < length()) storage[n] = item;
			else { setlength(n+1); storage[n] = item; }
		}
		
		inline int length() const { return vsize; }
		inline void reducelength(int by)
		{ setlength(vsize-by); }

		inline void setlength(int l, bool hardcap=false) {
			if (hardcap) {
				if (l!=asize) {
					storage = (T *)realloc(storage,l*sizeof(T));
				}
				vsize = l;
				asize = l;
			} else {
				if (l<vsize) {
					if (l<(asize>>2)) {
						int lim = l<<1;
						for(int c=asize>>1;c>lim;asize=c,c>>=1)
							;
						storage = (T *)realloc(storage,asize*sizeof(T));
					}
					vsize = l;
				} else if (l!=vsize) {
					if (l>asize) {
						for(;asize<l;asize<<=1)
							;
						storage = (T *)realloc(storage,asize*sizeof(T));
					}
					vsize = l;
				}
			}
		}	

		// careful using lists of integers, as this function will *not*
		//  remove the specified argument from the list!
		// Furthermore, if you have a pointer to a list, using -= will
		//  adjust the pointer and not the list, so be sure to dereference
		//  first!
		inline void operator-=(int by) 
			{ reducelength(by); }
		// be similarly careful with += !
		inline void operator+=(T item) {
			setlength(vsize+1); storage[vsize-1] = item;
		}
		inline T &operator[](int i)
			{ return storage[i]; }
		// this is the same as nth -- some may find the notation for confusing
		//   and some may find it easier -- it is designed to get around the
		//   whole l-value vs. r-value problem without requiring a proxy class
		//   which can really screw up things like reference counting.
		inline T operator()(int i) const
			{ return storage[i]; }

		inline int find(const T &item) const { int i,m=length();
			for(i=0;i<m;i++) if (storage[i]==item) 
				return i;
			return -1;
		}

		inline void sort(int (*comp)(const void *,const void *)) {
			qsort(storage,length(),sizeof(T),comp);
		}

		inline void swap(dlist<T> &l) {
			T *ts;
			int tvs;
			ts = storage;
			storage = l.storage;
			l.storage = ts;
			tvs = vsize;
			vsize = l.vsize;
			l.vsize = tvs;
		}

			
		inline int remove(const T &item) {
			int i,j=0,m=length();
			for(i=0;i<m;i++)
				if (storage[i]==item) {
					storage[i] = storage[--m];
					reducelength(1);
					j++;
					i--;
				}
			return j;
		}

		inline void del(int pos) {
			storage[pos] = storage[length()-1];
			reducelength(1);
		}

		inline int replace(const T &x, const T &y) { int i,j=0,m=length();
			for(i=0;i<m;i++) if (storage[i]==x) {
				storage[i] = y; j++;
			}
			return j;
		}
		inline dlist<T> *intersection(dlist<T> *l) const {
			dlist<T> *ret = new dlist<T>;
			int i,m=length();
		
			for(i=0;i<m;i++)
				if (l->find(storage[i])) (*ret) += storage[i];
			return ret;
		}
		inline int subsetof(dlist<T> *l) const {
			int i,m=length();

			for(i=0;i<m;i++)
				if (l->find(storage[i])==-1) return 0;
			return 1;
		}

		inline int identical(dlist<T> *l) const {
			return l->length()==length() && subsetof(l);
		}
		
		inline void operator+=(const dlist<T> &l) {
			int m=l.length();
			int i=length();
			setlength(m+i);
			for(int j=0;j<m;i++,j++)
				storage[i] = l.storage[j];
		}
		
		inline void insert(const T &item, int pos) {
			int m=length();
			(*this) += storage[m-1];
			for(int i=m-2;i>=pos;i++) {
				storage[i] = storage[i+1];
			}
			storage[pos] = item;
		}
		
	private:
		dlist(const dlist &x) { }

		T *storage;
		int vsize,asize;
	
};
 

template <class T>
dlist<T>::dlist (int initsize) {
	storage = (T *)malloc(initsize*sizeof(T));
	vsize=0;asize=initsize;
}

template <class T>
dlist<T>::dlist (T *s, int size, int realsize) {
	storage = s;
	vsize = size; asize = realsize;
}

template <class T>
dlist<T>::~dlist(void) {
	if (storage) free(storage);
}

template <class T>
inline ostream &operator<<(ostream &s, const dlist<T> &l) {
	for(int i=0;i<l.length();i++) {
		s << l(i); if (i!=l.length()) s << ','; else s << ':' << endl;
	}
	return s;
}

template <class T>
inline istream &operator>>(istream &s, dlist<T> &l) {
	T t;
	char c;

	do {
		s >> t;
		s >> c;
		l += t;
	} while (c==',');
}

#endif
