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
#ifndef _ILIST_H_
#define _ILIST_H_

#include <stdlib.h>
#include <iostream>
#include <malloc.h>

using namespace std;

template <class T>
class ilist {
	public:

		ilist (int initsize = 10, int expandrate = 2);
		ilist (T *s, int size);
		~ilist();
		T nth(int n) const { return storage[n]; }
		/* { if (n<0 || n>=virtualsize) return *bl; 
				else return storage[n]; };*/
		inline void setnth(const T &item, int n) {
			if (n < length()) storage[n] = item;
			else difficultset(item,n);
		}
		
		inline int length() const { return vsize; } // return storage ? _msize(storage)/sizeof(T) : 0; };
		inline void reducelength(int by)
		{ if (storage==NULL && by==0) return;
		  storage = (T *)realloc(storage,(vsize-by)*sizeof(T)); 
			vsize-=by; } // _msize(storage)-(by*sizeof(T))); }

		inline void setlength(int l)
		{ storage = (T *)realloc(storage,l*sizeof(T)); vsize = l; }

		// careful using lists of integers, as this function will *not*
		//  remove the specified argument from the list!
		// Furthermore, if you have a pointer to a list, using -= will
		//  adjust the pointer and not the list, so be sure to dereference
		//  first!
		inline void operator-=(int by) 
			{ reducelength(by); }
		// be similarly careful with += !
		inline void operator+=(T item) {
			difficultset(item,length());
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

		inline void swap(ilist<T> &l) {
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
		inline ilist<T> *intersection(ilist<T> *l) const {
			ilist<T> *ret = new ilist<T>;
			int i,m=length();
		
			for(i=0;i<m;i++)
				if (l->find(storage[i])!=-1) (*ret) += storage[i];
			return ret;
		}
		inline int subsetof(ilist<T> *l) const {
			int i,m=length();

			for(i=0;i<m;i++)
				if (l->find(storage[i])==-1) return 0;
			return 1;
		}

		inline int identical(ilist<T> *l) const {
			return l->length()==length() && subsetof(l);
		}
		
		inline void operator+=(const ilist<T> &l) {
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
		// no copy constructor allow -- its not what you want
		// to do in all probability!
		ilist(const ilist &x) { }

		void difficultset(const T &item, int n);

		T *storage;
		int vsize;
	
};
 

template <class T>
ilist<T>::ilist (int initsize, int expandrate) {
	storage = (T *)NULL;
	vsize=0;
}

template <class T>
ilist<T>::ilist (T *s, int size) {
	storage = s;
	vsize = size;
}

template <class T>
ilist<T>::~ilist(void) {
	if (storage) free(storage);
}

/*
template <class T>
T ilist<T>::nth (int n) const {
	return storage[n];
}
*/

template <class T>
void ilist<T>::difficultset(const T &item, int n) {
	storage = (T *)realloc(storage,(n+1)*sizeof(T));
	vsize = n+1;
	storage[n] = item;
}	

template <class T>
inline ostream &operator<<(ostream &s, const ilist<T> &l) {
	for(int i=0;i<l.length();i++) {
		s << l(i); if (i!=l.length()) s << ','; else s << ':' << endl;
	}
	return s;
}

template <class T>
inline istream &operator>>(istream &s, ilist<T> &l) {
	T t;
	char c;

	do {
		s >> t;
		s >> c;
		l += t;
	} while (c==',');
}

#endif
