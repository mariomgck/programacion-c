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
#ifndef _PQUEUE_H_
#define _PQUEUE_H_


#include <globaldef.h>

#include <containers/ilist.h>
#include <containers/dlist.h>

#include <iostream>
#include <math.h>

using namespace std;

// the class T needs to respond to the rank method (takes nothing returns a FP)
// and it also needs to respond to the setposition method (which is called whenever it is
//  moved in the queue)
// Other than that, this is a simple priority queue built on a heap
template <class T>
class pqueue {
public:
	// a blank queue
	inline pqueue(void) {};
	// a queue initialized with s (note that the queue now owns the memory of s!)
	inline pqueue(T *s, int size, int realsize) : heap(s,size,realsize) { donebatch(); }
	inline ~pqueue(void) {};

	inline T head() const { return heap.nth(0); }
	void add(T t);

	// add batch adds a new member but does not redo the heap (you must call donebatch
	//  before using head, remove, or nth)
	inline void addbatch(T t) { heap += t; }
	// adds a whole list of entries at once
	inline void addbatch(const dlist<T> &t) { heap += t; }
	// after this the heap is once again correct
	inline void donebatch() {
		heap.sort(sortem);
		for(int i=0;i<heap.length();i++) heap[i].setposition(i);
	}

	static int sortem(const void * const t1, const void * const t2) {
		FP v1 = ((T *)t1)->rank();
		FP v2 = ((T *)t2)->rank();
		if (v1<v2) return -1;
		if (v1>v2) return 1;
		return 0;
	}

	void remove(int i);

	inline int size() const { return heap.length(); }
	inline T nth(int n) const { return heap.nth(n); }

	template <class U>
	friend istream &operator>>(istream &s, pqueue<U> &q);
	template <class U>
	friend ostream &operator<<(ostream &s, const pqueue<U> &q);

private:

	int heapify(int i,FP tempk);

	dlist<T> heap;
};


template <class T>
void pqueue<T>::add(T t) {

	int i;
	FP key = t.rank();

	heap+=t;
	for (i=heap.length()-1; i>0 && key > heap[(i-1)/2].rank(); i=(i-1)/2) {
		heap[i] = heap[(i-1)/2];
		heap[i].setposition(i);
	}
	heap[i] = t;
	heap[i].setposition(i);
}

template <class T>
void pqueue<T>::remove(int i) {

	T t;
	FP key;

	if (i<0 || i>=heap.length()) return;
	if (i==0) {
		heap[0].setposition(-1);
		heap[0] = heap[heap.length()-1];
		heap[0].setposition(0);
		heap -= 1;
		heapify(0,heap[0].rank());
		return;
	}
	i=heapify(i,-hugenum);
	heap[i].setposition(-1);
	if (i==heap.length()-1) {
		heap -= 1;
		return;
	}
	key = heap[heap.length()-1].rank();
	t = heap[heap.length()-1];
	heap -= 1;
	while(i>0 && key > heap[(i-1)/2].rank()) {
		heap[i] = heap[(i-1)/2];
		heap[i].setposition(i);
		i = (i-1)/2;
	}
	heap[i] = t;
	heap[i].setposition(i);
}
		
template <class T>
int pqueue<T>::heapify(int i, FP tempk) {

	int l,r,large;
	T temp;
	FP largek;

	do {
		l = i*2+1;
		r = i*2+2;

		if (l < heap.length() && heap[l].rank() > tempk) {
			large = l;
			largek = heap[l].rank();
		} else {
			large = i;
			largek = tempk;
		}
		if (r < heap.length() && heap[r].rank() > largek) {
			large = r;
			largek = heap[r].rank();
		}

		if (large==i) {
			heap[i].setposition(i);
			return i;
		}

		temp = heap[i];
		heap[i] = heap[large];
		heap[large] = temp;

		heap[i].setposition(i);

		i = large;
	} while (1);
}

template <class T>
inline ostream &operator<<(ostream &s, const pqueue<T> &q) {
	return s << q.heap;
}

template <class T>
inline istream &operator>>(istream &s, pqueue<T> &q) {
	return s >> q.heap;
}
		


#endif
