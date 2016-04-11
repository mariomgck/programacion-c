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
#include <warping/warper2.h>

using namespace std;

warper2::warper2(const ilist<char *> &from, const ilist<char *> &to, 
			   const char *frompts, const char *topts, const warperparams & params) :
			p(params) {
	for(int i=0;i<from.length();i++)
		As += from(i);
	for(int i=0;i<to.length();i++)
		Bs += to(i);
	level = 0;
	ittnum = 0;
	p = params;
	A = NULL;
	Aold = NULL;
	B = NULL;
	affineA = matrix::eye;
	affineb = vec::zero;
	loadpts(frompts,topts);
	nextlevel();
}

warper2::~warper2()
{
	if (A!=NULL) delete A;
	if (Aold!=NULL) delete Aold;
	if (B!=NULL) delete B;
}

void warper2::loadpts(const char *frompts, const char *topts) {

	if (frompts==NULL || topts==NULL) return;
	ifstream fin;

	int ttlnumber;
	fin.open(frompts);
	if (!fin.good()) {
		return;
		//cout << "could not open " << frompts << endl;
		//exit(1);
	}
	fin >> ttlnumber;
	vec curr;
	int t, i;
	for(i = 0;i<ttlnumber;i++) {
		//fin >> t;
		//if (t) fin >> curr;
		//else curr = hugenum;
		fin >> curr;
		fromvec += curr;
	}
	fin.close();
	fin.open(topts);
	if (!fin.good()) {
		fromvec.setlength(0);
		return;
		//cout << "could not open " << frompts << endl;
		//exit(1);
	}
	fin >> i;
	if (i!=ttlnumber) {
		//cout << "disagreement on the number of user specified points" << endl;
		//exit(1);
		fromvec.setlength(0);
		tovec.setlength(0);
		return;
	}
	for(i=0;i<ttlnumber;i++) {
		//fin >> t;
		//if (t) fin >> curr;
		//else curr = hugenum;
		fin >> curr;
		tovec += curr;
	}
	for(i=0;i<tovec.length();) {
		if (!tovec[i].isvalid() || !fromvec[i].isvalid()) {
			//cout << "removing " << tovec[i] << " & " << fromvec[i] << endl;
			tovec.del(i);
			fromvec.del(i);
		} else {
			//cout << "keeping " << tovec[i] << " & " << fromvec[i] << endl;
			tovec[i][3] *= p.gamma;
			tovec[i][4] *= p.gamma;
			tovec[i][5] *= p.gamma;
			fromvec[i][3] *= p.gamma;
			fromvec[i][4] *= p.gamma;
			fromvec[i][5] *= p.gamma;
			i++;
		}
	}
	findaffine();
	fin.close();
}

// this finds a constrained affine transform (i.e. A must be block diagonal)
void warper2::findaffine() {

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
	//cout << "sanity check: " << endl;
	//(F*C).niceprint(cout) << endl;
	//(C*F).niceprint(cout) << endl;
	affineA = D*F;
	affineb = (t - affineA*s)/tovec.length();
	//cout << "found A: " << endl;
	//affineA.niceprint(cout) << endl;
	//cout << "found b: " << endl;
	//cout << affineb << endl;
	FP v=0,ov=0;
	//int j;
	for(int i=0;i<fromvec.length();i++) {
		//cout << "orig from: " << fromvec[i] << endl;
		//cout << "orig to:   " << tovec[i] << endl;
		//cout << "error:     " << (fromvec[i]-tovec[i]).len() << endl;
		ov += (fromvec[i]-tovec[i]).len2();
		fromvec[i] = affineA*fromvec[i] + affineb;
		//cout << "new from:  " << fromvec[i] << endl;
		//cout << "new error: " << (fromvec[i]-tovec[i]).len() << endl;
		v += (fromvec[i]-tovec[i]).len2();
		//cin >> j;
	}
	cout << "total old squared error: " << ov << endl;
	cout << "total new squared error: " << v << endl;
}
	
bool warper2::nextlevel() {

	if (level>=As.length() && level>=Bs.length()) return false;

	if (level<As.length()) {
		warpmesh2 *next;
		trimesh *nextold;
		ifstream f(As[level]);
		if (!f.good()) {
			cout << "could not open file " << As[level] << endl;
			exit(1);
		}
		cout << "loading A: " << As[level] << endl;
		
		warpparms prms(p,level*p.t);
		next = new warpmesh2(f,prms,p.gamma,affineA,affineb);
		
		f.close();

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


bool warper2::iterate() {

	if (level>=As.length() && level>=Bs.length()) {
		if (ittnum>=p.t0) return false;
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

void warper2::snap() {

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
