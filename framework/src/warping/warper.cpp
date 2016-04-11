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
#include <iostream>
#include <fstream>
#include <warping/warper.h>
#include <warping/warpmesh.h>
#include <warping/warpmesh2.h>
#include <warping/warpmesh3.h>
#include <warping/warpmesh4.h>

using namespace std;

warper::warper(const ilist<char *> &from, const ilist<char *> &to, 
	const char * const frompts, const char * const topts,
	const warperparams & params, int vernum) :
			p(params) {
	for(int i=0;i<from.length();i++)
		As += from(i);
	for(int i=0;i<to.length();i++)
		Bs += to(i);
	level = 0;
	ittnum = 0;
	A = NULL;
	Aold = NULL;
	B = NULL;
	affineA = matrix::eye;
	affineb = vec::zero;
	loadpts(frompts,topts);
	warpv = vernum;
	nextlevel();
}

warper::warper(const ilist<trimesh *> *from, const ilist<trimesh *> *to, 
	const ilist<vec> &frompts, const ilist<vec> &topts,
	const warperparams & params, int vernum) :
			p(params) {

	affineA = matrix::eye;
	affineb = vec::zero;
	loadpts(frompts,topts);
	Ams = (ilist<trimesh *> *)from;
	Bms = (ilist<trimesh *> *)to;
	level = 0;
	ittnum = 0;
	A = NULL;
	Aold = NULL;
	B = NULL;
	warpv = vernum;
	nextlevel();
}

warper::warper(const ilist<char *> &from, const ilist<trimesh *> *to, 
	const char * const frompts, const ilist<vec> &topts,
	const warperparams & params, int vernum) :
			p(params) {

	affineA = matrix::eye;
	affineb = vec::zero;
	loadpts(frompts,topts);
	for(int i=0;i<from.length();i++)
		As += from(i);
	Bms = (ilist<trimesh *> *)to;
	level = 0;
	ittnum = 0;
	A = NULL;
	Aold = NULL;
	B = NULL;
	warpv = vernum;
	nextlevel();
}

warper::warper(const ilist<trimesh *> *from, const ilist<char *> &to, 
	const ilist<vec> &frompts, const char * const topts,
	const warperparams & params, int vernum) :
			p(params) {

	affineA = matrix::eye;
	affineb = vec::zero;
	loadpts(frompts,topts);
	Ams = (ilist<trimesh *> *)from;
	for(int i=0;i<to.length();i++)
		Bs += to(i);
	level = 0;
	ittnum = 0;
	A = NULL;
	Aold = NULL;
	B = NULL;
	warpv = vernum;
	nextlevel();
}

warper::~warper()
{
	if (A!=NULL) delete A;
	if (Aold!=NULL) delete Aold;
	if (B!=NULL) delete B;
}

void warper::loadpts(const ilist<vec> &frompts, const ilist<vec> &topts) {

	fromvec += frompts;
	tovec += topts;
}

void warper::loadpts(const ilist<vec> &frompts, const char * const topts) {

	if (topts==NULL) return;
	ifstream fin;

	fromvec += frompts;

	int ttlnumber;
	vec curr;
	fin.open(topts);
	if (!fin.good()) {
		fromvec.setlength(0);
		return;
	}
	fin >> ttlnumber;
	for(int i=0;i<ttlnumber;i++) {
		fin >> curr;
		tovec += curr;
	}
	fin.close();
	completeload();
}

void warper::loadpts(const char * const frompts, const ilist<vec> &topts) {

	if (frompts==NULL) return;
	ifstream fin;

	int ttlnumber;
	fin.open(frompts);
	if (!fin.good()) return;
	fin >> ttlnumber;
	vec curr;
	for(int i=0;i<ttlnumber;i++) {
		fin >> curr;
		fromvec += curr;
	}
	fin.close();
	tovec += topts;
	completeload();
}

void warper::loadpts(const char * const frompts, const char * const topts) {

	if (frompts==NULL || topts==NULL) return;
	ifstream fin;

	int ttlnumber;
	fin.open(frompts);
	if (!fin.good()) return;
	fin >> ttlnumber;
	vec curr;
	for(int i=0;i<ttlnumber;i++) {
		fin >> curr;
		fromvec += curr;
	}
	fin.close();
	fin.open(topts);
	if (!fin.good()) {
		fromvec.setlength(0);
		return;
	}
	fin >> ttlnumber;
	for(int i=0;i<ttlnumber;i++) {
		fin >> curr;
		tovec += curr;
	}
	fin.close();
	completeload();
}

void warper::completeload(void) {

	if (fromvec.length()!=tovec.length()) {
		fromvec.setlength(0);
		tovec.setlength(0);
		return;
	}
	for(int i=0;i<tovec.length();) {
		if (!tovec[i].isvalid() || !fromvec[i].isvalid()) {
			tovec.del(i);
			fromvec.del(i);
		} else {
			for(int j=3;j<DIM;j++) {
				tovec[i][j] *= p.gamma;
				fromvec[i][j] *= p.gamma;
			}
			i++;
		}
	}
	findaffine();
}

// this finds a constrained affine transform (i.e. A must be block diagonal)
void warper::findaffine() {

	matrix C(0,false),D(0,false),F;
	vec s(0),t(0);

	for(int i=0;i<tovec.length();i++) {
		s += fromvec[i];
		t += tovec[i];
		C += matrix(fromvec[i],fromvec[i]);
		D += matrix(tovec[i],fromvec[i]);
	}
	D -= matrix(t,s)/tovec.length();
	C -= matrix(s,s)/tovec.length();
	
	// to make it block diagonal -- remove off block diagonal elements!
	// (yes -- this is the correct thing to do, not just a hack)
	for(int i=0;i<3;i++)
		for(int j=3;j<DIM;j++)
			C[i][j] = C[j][i] = D[i][j] = D[j][i] = 0;
	
	
	bool w;
	F = C.inv(w);
	if (!w) {
		cout << "having to use pseudo-inverse" << endl;
		F = (C.t()*C).inv(w);
		if (!w) return; // failure -- drat!
		F = F*C.t();
	}
	affineA = D*F;
	affineb = (t - affineA*s)/tovec.length();
	FP v=0,ov=0;
	for(int i=0;i<fromvec.length();i++) {
		ov += (fromvec[i]-tovec[i]).len2();
		fromvec[i] = affineA*fromvec[i] + affineb;
		v += (fromvec[i]-tovec[i]).len2();
	}
	cout << "total old squared error: " << ov << endl;
	cout << "total new squared error: " << v << endl;
}
	
bool warper::nextlevel() {

	//cout << "next level" << endl;
	if (level>=Ams->length() && level>=As.length() &&
	    level>=Bms->length() && level>=Bs.length()) return false;

	if (level<As.length()) {
		awarpmesh *next;
		trimesh *nextold;
		ifstream f(As[level]);
		if (!f.good()) {
			cout << "could not open file " << As[level] << endl;
			exit(1);
		}
		cout << "loading A: " << As[level] << endl;
		
		warpparms prms(p,level*p.t);
		cout << "warpmethod = " << warpv << endl;
		switch(warpv) {
			case 1:
				next = new warpmesh(f,prms,p.gamma,
					affineA,affineb);
			break;
			case 2:
				next = new warpmesh2(f,prms,p.gamma,
					affineA,affineb);
			break;
			case 3:
				next = new warpmesh3(f,prms,p.gamma,
					affineA,affineb);
			break;
			case 0:
			case 4:
				next = new warpmesh4(f,prms,p.gamma,
					affineA,affineb);
			break;
			default:
				cout << "Unknown minimization method (" <<
					 warpv << ") -- terminating" << endl;
				exit(1);
			break;
		}
		
		f.close();

		nextold = next->dup();

		if (A!=NULL) {
			next->warp(A->D(),Aold);
			delete Aold;
			delete A;
		}
		A = next;
		Aold = nextold;
	} else if (level<Ams->length()) {
		awarpmesh *next;
		trimesh *nextold;
		warpparms prms(p,level*p.t);
		switch(warpv) {
			case 1:
				next = new warpmesh(*Ams->nth(level),prms,
					p.gamma,affineA,affineb);
			break;
			case 2:
				next = new warpmesh2(*Ams->nth(level),prms,
				 	p.gamma,affineA,affineb);
			break;
			case 3:
				next = new warpmesh3(*Ams->nth(level),prms,
				 	p.gamma,affineA,affineb);
			break;
			case 0:
			case 4:
				next = new warpmesh4(*Ams->nth(level),prms,
				 	p.gamma,affineA,affineb);
			break;
			default:
				cout << "Unknown minimization method (" <<
					warpv << ") -- terminating" << endl;
				exit(1);
			break;
		}
		nextold = next->dup();
		if (A!=NULL) {
			next->warp(A->D(),Aold);
			delete Aold;
			delete A;
		}
		A = next;
		Aold = nextold;
	}
	if (level<Bs.length()) {
		if (B!=NULL) delete B;
		cout << "loading B: " << Bs[level] << endl;
		ifstream f2(Bs[level]);
		if (!f2.good()) {
			cout << "could not open file " << Bs[level] << endl;
			exit(1);
		}
		B = new trimesh(f2,p.gamma);
	} else if (level<Bms->length()) {
		if (B!=NULL) delete B;
		B = new trimesh(*Bms->nth(level),p.gamma);
	}
	A->removesprings();
	vec toplace;
	FP d2;
	face *f;
	vertex *v1,*v2;
	for(int i=0;i<fromvec.length();i++) {
		toplace = B->closestpoint(tovec[i],d2,f,v1,v2);
		A->adduserspring(fromvec[i],toplace);
	}
	level++;
	ittnum=0;
	return true;
}


bool warper::iterate() {

	if (level>=Ams->length() && level>=As.length() &&
	    level>=Bms->length() && level>=Bs.length()) {
		if (ittnum>=p.t0) {
			//A->printlevelstats();
			//B->printlevelstats();
			return false;
		}
	} else {
		if (ittnum>=p.t)
			nextlevel();
	}
	A->iterate(B);
/*
	if (ittnum%10==9 || ittnum%p.t==0 || ittnum==p.t-1 || ittnum==p.t0-1) {
		char buffer[100];
		sprintf(buffer,"output.%d.%d.tm",level,ittnum);
		ofstream of(buffer);
		A->saveit(of,trimesh::tmff1);
	}
*/
	ittnum++;
	return true;
}

void warper::snap() {

	vec p,newp;
	face *f;
	vertex *v1,*v2;
	FP d2;

	for(int i=0;i<A->numvertices();i++) {
		p = A->vertexpos(i);
		newp = B->closestpoint(p,d2,f,v1,v2);
		A->movevertex(i,newp);
	}
}
