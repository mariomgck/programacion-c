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
#include "facehash.h"
		
facehash::facehash(vec &min, vec &max, vec &divsize, int l) {
	level = l;
	initialize(min,max,divsize);

}

facehash::~facehash() {
	if (!storage) return;
	for(int i=0;i<asize;i++) {
		if (storage[i].contents)
			delete storage[i].contents;
		if (storage[i].subhash)
			delete storage[i].subhash;
	}
	delete []storage;
}

facehash::facehash(int l) {
	level = l;
	storage = NULL;
	asize = 0;
	vsize = 0;
}

void facehash::initialize(vec &min, vec &max, vec &divsize) {
	//cout << "starting initialize (" << level << ')' << endl;
	asize = 16; // must be a power of 2!
	vsize = 0;
	storage = new cell[16]; // same number here!
	mincorner = min;
	maxcorner = max;
	for(int i=0;i<DIM;i++) 
		div[i] = (max[i]-min[i])/ncells[i];
}

#define truncit(f,s,i)  \
	if (f<0) s = 0; \
	else if (f>=ncells[i]) s = ncells[i]-1; \
	else s = (short)floor(f);

void facehash::dorange(face *f,void (*fn)(short *, face *)) {
	vec minv=f->min(),maxv=f->max();
	vec t1,t2;
	f->tangents(t1,t2);
	vec c = f->basepoint();
	c -= (c*t1)*t1 + (c*t2)*t2;
	FP currhigh,currlow;

	int d1,d2;
	vec diff = maxv-minv;
	int i,j;
	for(i=0;i<DIM;i++) diff[i] /= div[i];

	if (diff[0]>diff[1]) { d1 = 0; d2 = 1; }
	else { d1 = 1; d2 = 0; }
	for(i=2;i<DIM;i++) {
		if (diff[i]>diff[d1]) {
			d2 = d1;
			d1 = i;
		} else if (diff[i]>diff[d2]) {
			d2 = i;
		}
	}
	FP t11,t12,t21,t22;
	t11 = t1[d1];
	t12 = t2[d1];
	t21 = t1[d2];
	t22 = t2[d2];
	FP t = t11*t22-t12*t21;
	FP c1 = c[d1];
	FP c2 = c[d2];
	vec x0 = (t12*c2*t1 + t21*c1*t2 - t22*c1*t1 - t11*c2*t2)/t + c;
	vec x1 = (t22*t1-t21*t2)/t;
	vec x2 = (t11*t2-t12*t1)/t;
	x0 += x1*mincorner[d1];
	x0 += x2*mincorner[d2];
	x1 *= div[d1];
	x2 *= div[d2];
	short mn1,mn2,mx1,mx2;
	currlow = (minv[d1]-mincorner[d1])/div[d1];
	currhigh = (maxv[d1]-mincorner[d1])/div[d1];
	truncit(currlow,mn1,d1);
	truncit(currhigh,mx1,d1);
	currlow = (minv[d2]-mincorner[d2])/div[d2];
	currhigh = (maxv[d2]-mincorner[d2])/div[d2];
	truncit(currlow,mn2,d2);
	truncit(currhigh,mx2,d2);
	vec cnr1,cnr2,cnr3,cnr4,cmx,cmn;
	short low[DIM],high[DIM];
	short in[DIM];
	int d,k;
	for(i=mn1;i<mx1;i++) {
		low[d1] = i;
		high[d1] = i;
		in[d1] = i;
		for(j=mn2;j<mx2;j++) {
			low[d2] = j;
			high[d2]= j;
			in[d2] = j;
			cnr1 = x0 + i*x1 + j*x2;
			cnr2 = cnr1 + x1;
			cnr3 = cnr2 + x2;
			cnr4 = cnr1 + x2;
			cmx = max(max(cnr1,cnr2),max(cnr3,cnr4));
			cmn = min(min(cnr1,cnr2),min(cnr3,cnr4));
			for(k=0;k<DIM;k++) {
				if (k==d1 || k==d2) continue;
				currlow = (cmn[k]-mincorner[k])/div[k];
				currhigh = (cmx[k]-mincorner[k])/div[k];
				truncit(currlow,low[k],k);
				in[k] = low[k];
				truncit(currhigh,high[k],k);
			}
			k=0;
			while(k<DIM) {
				fn(in,f);
				for(k=0;k<DIM&&in[k]==high[k];k++)
					in[k] = low[k];
				if (k<DIM) in[k]++;
			}
		}
	}
}

void facehash::add(face *t) {
	if (!toaddset.ismember((long)t)) {
		toadd += t;
		toaddset.add((long)t);
	} else if (toremoveset.ismember((long)t)) {
		toremoveset.remove((long)t);
	}
}

void facehash::addall() {
	//cout << "ADDING ALL (" << level << ") " << toadd.length() << endl;
	for(int i=0;i<toadd.length();i++) {
		if (toremoveset.ismember((long)toadd.nth(i))) continue;
		dorange(toadd[i], 
		(void (*)(short *,face *))&facehash::addcell);
	}
	toadd.setlength(0);
	toaddset.clear();
	toremoveset.clear();
	
}

void facehash::remove(face *t) {
	if (toaddset.ismember((long)t)) {
		toremoveset.add((long)t);
	} else {
		dorange(t,(void (*)(short *,face *))&facehash::removecell);
	}
}

void facehash::addcell(short i[DIM], face *t) {
	unsigned int p = hashfn(i,true);
	if (storage[p].contents) {
		if (level<maxlevel &&
		    storage[p].contents->length()==maxcontents) {
			int d;
			vec mn = mincorner;
			for(d=0;d<DIM;d++)
				mn += div[d]*i[d];
			vec mx = mn+div;
			storage[p].subhash = new facehash(mn,mx,mx,level+1);
			for(d=0;d<storage[p].contents->length();d++)
				storage[p].subhash->
					add(storage[p].contents->nth(d));
			storage[p].subhash->add(t);
			delete storage[p].contents;
			storage[p].contents = NULL;
		} else *(storage[p].contents) += t;
	} else if (storage[p].subhash) {
		storage[p].subhash->add(t);
	} else {
		storage[p].contents = new ilist<face *>;
		*(storage[p].contents) += t;
	}
}

void facehash::removecell(short i[DIM],face *t) {
	unsigned int p = hashfn(i,false);
	if (p!=-1) {
		if (storage[p].contents) storage[p].contents->remove(t);
		else storage[p].subhash->remove(t);
	}
}

float facehash::castray(vec &x0, vec &v,
                        geoiterateface &itt) {
	cout << "castray not yet implemented with facehash" << endl;
	return 0;
}

float facehash::castray3d(vec &x0, vec &v,
                          geoiterateface &itt) {
	cout << "castray3d not yet implemented with facehash" << endl;
	return 0;
}

float facehash::search(vec &center, geoiterateface &itt)  {
	FP best = itt.initaldist();
	FP curr;
	qcell *next,*consider;
	pqueue<qcellptr> q;
	longset ls;
	q.add(next = tocell(center,ls));
	FP currdist = next->celldist;

	if (level==0) numq++;
	int p,i;

	// first take care of the "patches" in toadd
	int pos = hashfn(next->i,false);
	if (toadd.length()>0 &&
	    (pos == -1 || storage[pos].contents==NULL ||
	     storage[pos].contents->length()<toadd.length())) {
		addall();
	} else {
		for(i=0;i<toadd.length();i++) {
			nums++;
			curr = itt.process(toadd[i]);
			if (curr<best) best = curr;
		}
	}

	while(q.size()>0) {
		next = q.head().p;
		q.remove(0);
		currdist = next->celldist;
		if (currdist>best) { delete next; break; }
		p = hashfn(next->i,false);
		if (p!=-1) {
			if (storage[p].contents) {
				for(i=0;i<storage[p].contents->length();i++) {
					nums++;
					curr = itt.process((*(storage[p].contents))[i]);
					if (curr<best) best = curr;
					if (curr<currdist) break;
				}
			} else {
				curr = storage[p].subhash->search(center,itt);
				if (curr<best) best = curr;
			}
			if (best<currdist) break;
		}
		for(i=0;i<DIM;i++) {
			if ((consider=tocell(center,next->i,i,-1,ls))!=NULL) {
				q.add(consider);
			}
			if ((consider=tocell(center,next->i,i,1,ls))!=NULL) {
				q.add(consider);
			}
		}
		delete next;
	}
	for(i=0;i<q.size();i++)
		delete q.nth(i).p;
	return best;
}
		

	
