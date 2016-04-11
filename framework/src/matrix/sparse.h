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
#ifndef _SPARSE_H_
#define _SPARSE_H_

//#include <globaldef.h>
#include <containers/ilist.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>


using namespace std;

class svec{ // this actually isn't a sparse vector, it is just
		// used with the sparse matrix
public:


	svec(int size) { 
		x = new double[size]; 
		s=size; 
	}

	svec(const svec &v) {
		s = v.s; x = new double[s];

		for(int i=0; i<s; i++) 
			x[i] = v.x[i];
	}

	svec(const svec &v, double mu) {
		s = v.s; x = new double[s];
		for(int i=0;i<s;i++) 
		x[i] = v.x[i]*mu;
	}

	~svec() { delete []x; }

	bool isfinite(void) 
	{
		for(int i=0;i<s;i++) if (!finite(x[i])) return false;
		return true;
	}


	inline double len2() { 
		double ret=0; 
		int i;
		for(i=0;i<s;i++) 
			ret += x[i]*x[i];
		return ret;
	}

	inline svec &operator-=(const svec &y) 
	{

		for(int i=0; i<s;i++) 
			x[i] -= y.x[i];

		return *this;
	}



	inline void clear() {
		for(int i=0;i<s;i++) x[i] = 0;
	}

	inline void addmult(const svec &v, double mu) {
		for(int i=0;i<s;i++) x[i] += v.x[i]*mu;
	}

	inline ostream &print(ostream &st) const {
		for(int i=0;i<s;i++) st << x[i] << ' ';
		return st << endl;
	}

	inline bool isvalid() {
		for(int i=0;i<s;i++) if (!finite(x[i])) return false;
		return true;
	}
	
	int s;
	double *x;
};

class smatrix {

public:
	smatrix(int M, int N, int maxel) { 
		int i;
		n = N; m = M;
		x = new double *[m];
		in = new int *[m];
		me = maxel;
		for(i=0;i<m;i++) {
			x[i] = new double[maxel];
			in[i] = new int[maxel+1];
			in[i][0] = 0;
		}
	}

	smatrix(const smatrix &ma) {
		int i,j;
		n = ma.n; m = ma.m;
		me = ma.me;
		x = new double *[m];
		in = new int *[m];
		for(i=0;i<m;i++) {
			x[i] = new double[me];
			in[i] = new int[me+1];
			in[i][0] = ma.in[i][0];
			for(j=0;j<me;j++) {
				in[i][j+1] = ma.in[i][j+1];
				x[i][j] = ma.x[i][j];
			}
		}
	}

	smatrix(const smatrix &ma, double mu) {
		int i,j;
		n = ma.n; m = ma.m;
		me = ma.me;
		x = new double *[m];
		in = new int *[m];
		for(i=0;i<m;i++) {
			x[i] = new double[me];
			in[i] = new int[me+1];
			in[i][0] = ma.in[i][0];
			for(j=0;j<me;j++) {
				in[i][j+1] = ma.in[i][j+1];
				x[i][j] = ma.x[i][j]*mu;
			}
		}
	}

	~smatrix() {
		int i;
		for(i=0;i<m;i++) {
			delete []in[i];
			delete []x[i];
		}
		delete []in;
		delete []x;
	}

	inline void clear() {
		for(int i=0;i<m;i++)
			in[i][0] = 0;
	}

	inline void setrow(int rownum, int nel, int *indices, double *row) {
		int i;
		in[rownum][0] = nel;
		for(i=0;i<nel;i++) {
			in[rownum][i+1] = indices[i];
			x[rownum][i] = row[i];
		}
	}

	inline void set(int i, int j, double el) {
		int k;
		for(k=0;k<in[i][0] && in[i][k+1] != j; k++);
		if (k==in[i][0]) {
			if (++in[i][0] > me)
			   cout << "sparse matrix overflow" << endl;
			in[i][k+1] = j;
			x[i][k] = el;
		} else {
			x[i][k] = el;
		}
	}
	inline void add(int i, int j, double el) {
		int k;
		for(k=0;k<in[i][0] && in[i][k+1] != j; k++);
		if (k==in[i][0]) {
			if (++in[i][0] > me)
			   cout << "sparse matrix overflow" << endl;
			in[i][k+1] = j;
			x[i][k] = el;
		} else {
			x[i][k] += el;
		}
	}
	smatrix &operator+=(const smatrix &ma) {
		if (ma.n != n || ma.m != m) {
			cout << "improper smatrix addition!" << endl;
			exit(1);
		}
		for(int i=0;i<ma.m;i++)
			addrow(i,ma.in[i][0],ma.in[i]+1,ma.x[i]);
		return *this;
	}
		
	smatrix operator+(const smatrix &ma) {
		return smatrix(*this)+=ma;
	}

	
	inline void addrow(int rownum, int nel, int *indices, double *row) {
		int i,j;
		for(i=0;i<nel;i++) {
			for(j=0;j<in[rownum][0] &&
			        in[rownum][j+1] != indices[i]; j++);
			if (j==in[rownum][0]) {
				if (++in[rownum][0] > me)
				   cout << "sparse matrix overflow" << endl;
				in[rownum][j+1] = indices[i];
				x[rownum][j] = row[i];
			} else {
				x[rownum][j] += row[i];
			}
		}
	}
	inline void subrow(int rownum, int nel, int *indices, double *row) {
		int i,j;
		for(i=0;i<nel;i++) {
			for(j=0;j<in[rownum][0] &&
			        in[rownum][j+1] != indices[i]; j++);
			if (j==in[rownum][0]) {
				if (++in[rownum][0] > me)
				   cout << "sparse matrix overflow" << endl;
				in[rownum][j+1] = indices[i];
				x[rownum][j] = -row[i];
			} else {
				x[rownum][j] -= row[i];
			}
		}
	}

	// right multiply
	inline void rightmult(const svec &b, svec &ret, const svec &w) const {
		int i,j;

		for(i=0;i<m;i++) {
			ret.x[i] = 0;
			for(j=0;j<in[i][0];j++) {
				ret.x[i] += b.x[in[i][j+1]]*x[i][j]*w.x[i];
			}
		}
	}
	inline void rightmult(const svec &b, svec &ret) const {
		int i,j;

		for(i=0;i<m;i++) {
			ret.x[i] = 0;
			for(j=0;j<in[i][0];j++) {
				ret.x[i] += b.x[in[i][j+1]]*x[i][j];
			}
		}
	}

	inline void leftmult(const svec &b, svec &ret, const svec &w) const {
		int i,j;

		for(i=0;i<n;i++) ret.x[i] = 0;
		for(i=0;i<m;i++)
			for(j=0;j<in[i][0];j++) {
				ret.x[in[i][j+1]] += b.x[i]*x[i][j]*w.x[i];
			}
	}

	inline void solve(const svec &b, svec &x, const svec &w, double tol) const {
		svec rc(m), gc(n), gp(n), dc(n), tc(m);
		int i,k;

		rightmult(x,rc,w);
		for(i=0;i<m;i++) rc.x[i]-=b.x[i]*w.x[i];
		double rssb = rc.len2();
		leftmult(rc,gc,w);
		for(i=0;i<n;i++) gc.x[i] = -gc.x[i];
		double oldgm2 = 100000;
		for(k=0;k<n;k++) {
			double gm2 = gc.len2();
			if (k==0) tol = gm2*tol;
			else if (gm2<tol) {
				//cout << "error small (" << k << "," << gm2 << ")" << endl;
				break;
			}
			if (gm2>oldgm2 && k>40) {
				//cout << "error increase ("<<k<<") ["<<gm2<<","<<oldgm2<<"]"<<endl;
				break;
			}
			oldgm2 = gm2;
			if (k>0) {
				double bi = gm2/gp.len2();
				for (i=0;i<n;i++) dc.x[i] = gc.x[i]+bi*dc.x[i];
			} else {
				for (i=0;i<n;i++) dc.x[i] = gc.x[i];
			}
			rightmult(dc,tc,w);
			double ai=gm2/tc.len2();
			for(i=0;i<n;i++) x.x[i] += ai*dc.x[i];
			for(i=0;i<m;i++) rc.x[i] += ai*tc.x[i];
			for(i=0;i<n;i++) gp.x[i]=gc.x[i];
			leftmult(rc,gc,w);
			for(i=0;i<n;i++) gc.x[i] = -gc.x[i];
		}
		//if (k==n) cout << "error BIG" << endl;
	}

	inline int getm(void) { return m; }
	inline int getn(void) { return n; }

	inline bool isfinite() {
		for(int i=0;i<m;i++)
			for(int j=0;j<in[i][0];j++)
				if (!finite(x[i][j+1])) return false;
		
		return true;
	}

private:
	int n,m,me;
	double **x;
	int **in;
};

class smatrix2 {

public:
	smatrix2(int M, int N) { 
		int i;
		n = N; m = M;
		x = new ilist<double>[m];
		in = new ilist<int>[m];
	}

	smatrix2(const smatrix2 &ma) {
		int i;
		n = ma.n; m = ma.m;
		x = new ilist<double>[m];
		in = new ilist<int>[m];
		for(i=0;i<m;i++) {
			x[i] += ma.x[i];
			in[i] += ma.in[i];
		}
	}

	smatrix2(const smatrix2 &ma, double mu) {
		int i,j;
		n = ma.n; m = ma.m;
		x = new ilist<double>[m];
		in = new ilist<int>[m];
		for(i=0;i<m;i++) {
			in[i] += ma.in[i];
			x[i] += ma.x[i];
			for(j=0;j<x[i].length();j++)
				x[i][j] *= mu;
		}
	}

	~smatrix2() {
		int i;
		delete []in;
		delete []x;
	}

	inline void clear() {
		for(int i=0;i<m;i++) {
			in[i].setlength(0);
			x[i].setlength(0);
		}
	}

	inline void setrow(int rownum, int nel, int *indices, double *row) {
		int i;
		in[rownum].setlength(0);
		x[rownum].setlength(0);
		for(i=0;i<nel;i++) {
			if (row[i] == 0) continue;
			in[rownum] += indices[i];
			x[rownum] += row[i];
		}
	}

	inline void set(int i, int j, double el) {
		if (el==0) return;
		int k;
		for(k=0;k<in[i].length() && in[i][k] != j; k++);
		if (k==in[i].length()) {
			in[i] += j;
			x[i] += el;
		} else {
			x[i][k] = el;
		}
	}
	inline void add(int i, int j, double el) {
		if (el==0) return;
		int k;
		for(k=0;k<in[i].length() && in[i][k] != j; k++);
		if (k==in[i].length()) {
			in[i] += j;
			x[i] += el;
		} else {
			x[i][k] += el;
		}
	}

	inline void addmult(const smatrix2 &ma, double mu) {
		if (ma.n != n || ma.m != m) {
			cout << "improper smatrix2 addition!" << endl;
			exit(1);
		}
		for(int i=0;i<ma.m;i++)
			addrowmult(i,ma.in[i],ma.x[i],mu);
	}
		
	inline smatrix2 &operator+=(const smatrix2 &ma) {
		if (ma.n != n || ma.m != m) {
			cout << "improper smatrix2 addition!" << endl;
			exit(1);
		}
		for(int i=0;i<ma.m;i++)
			addrow(i,ma.in[i],ma.x[i]);
		return *this;
	}
		
	inline smatrix2 operator+(const smatrix2 &ma) {
		return smatrix2(*this)+=ma;
	}

	inline void addrow(int rownum, int nel, int *indices, double *row) {
		int i,j;
		for(i=0;i<nel;i++) {
			if (row[i]==0) continue;
			for(j=0;j<in[rownum].length() &&
			        in[rownum][j] != indices[i]; j++);
			if (j==in[rownum].length()) {
				in[rownum] += indices[i];
				x[rownum] += row[i];
			} else {
				x[rownum][j] += row[i];
			}
		}
	}
	inline void addrowmult(int rownum,
		const ilist<int> &indices, const ilist<double> &row, double mu) {
		int i,j;
		for(i=0;i<indices.length();i++) {
			if (row(i)==0) continue;
			for(j=0;j<in[rownum].length() &&
			        in[rownum][j] != indices(i); j++);
			if (j==in[rownum].length()) {
				in[rownum] += indices(i);
				x[rownum] += row(i)*mu;
			} else {
				x[rownum][j] += row(i)*mu;
			}
		}
	}
	inline void addrow(int rownum,
		const ilist<int> &indices, const ilist<double> &row) {
		int i,j;
		for(i=0;i<indices.length();i++) {
			if (row(i)==0) continue;
			for(j=0;j<in[rownum].length() &&
			        in[rownum][j] != indices(i); j++);
			if (j==in[rownum].length()) {
				in[rownum] += indices(i);
				x[rownum] += row(i);
			} else {
				x[rownum][j] += row(i);
			}
		}
	}
	inline void subrow(int rownum, int nel, int *indices, double *row) {
		int i,j;
		for(i=0;i<nel;i++) {
			if (row[i]==0) continue;
			for(j=0;j<in[rownum].length() &&
			        in[rownum][j] != indices[i]; j++);
			if (j==in[rownum].length()) {
				in[rownum] += indices[i];
				x[rownum] += -row[i];
			} else {
				x[rownum][j] -= row[i];
			}
		}
	}

	// right multiply
	inline void rightmult(const svec &b, svec &ret, const svec &w) const {
		int i,j;

		for(i=0;i<m;i++) {
			ret.x[i] = 0;
			for(j=0;j<in[i].length();j++) {
				ret.x[i] += b.x[in[i][j]]*x[i][j]*w.x[i];
			}
		}
	}
	inline void rightmult(const svec &b, svec &ret) const {
		int i,j;

		for(i=0;i<m;i++) {
			ret.x[i] = 0;
			for(j=0;j<in[i].length();j++) {
				ret.x[i] += b.x[in[i][j]]*x[i][j];
			}
		}
	}

	inline void leftmult(const svec &b, svec &ret, const svec &w) const {
		int i,j;

		for(i=0;i<n;i++) ret.x[i] = 0;
		for(i=0;i<m;i++)
			for(j=0;j<in[i].length();j++) {
				ret.x[in[i][j]] += b.x[i]*x[i][j]*w.x[i];
			}
	}
	inline void leftmult(const svec &b, svec &ret) const {
		int i,j;

		for(i=0;i<n;i++) ret.x[i] = 0;
		for(i=0;i<m;i++)
			for(j=0;j<in[i].length();j++) {
				ret.x[in[i][j]] += b.x[i]*x[i][j];
			}
	}

	inline void print(ostream &s) const {
		for(int i=0;i<m;i++) {
			for(int j=0;j<n;j++) {
				int k;
				for(k=0;k<in[i].length() &&
				            in[i][k]!=j;k++) ;
				if (k==in[i].length()) s << "0 ";
				else s << x[i][k] << ' ';
			}
			cout << endl;
		}
	}

	inline void solve(const svec &b, svec &x, double tol) const {
		svec rc(m), gc(n), dc(n), tc(m);
		int i,k;

		rightmult(x,rc);
		for(i=0;i<m;i++) rc.x[i]-=b.x[i];
		double move2,move2lim;
		double gc2save;
		for(k=0;k<n;k++) {
			leftmult(rc,gc);
			for(i=0;i<n;i++) gc.x[i] = -gc.x[i];
			double gm2 = gc.len2();
			if (k>0) {
				double bi = gm2/gc2save;
				for (i=0;i<n;i++) dc.x[i] = gc.x[i]+bi*dc.x[i];
			} else {
				for (i=0;i<n;i++) dc.x[i] = gc.x[i];
			}
			gc2save = gm2;
			rightmult(dc,tc);
			double ai=gm2/tc.len2();
			for(i=0;i<n;i++) x.x[i] += ai*dc.x[i];
			move2 = dc.len2()*ai*ai;
			if (k==0) {
				if (move2==0.0) break;
				move2lim = move2*tol;
			} else if (move2<move2lim) break;
			for(i=0;i<m;i++) rc.x[i] += ai*tc.x[i];
		}
	}

	// don't use this one!
	inline void solve(const svec &b, svec &x, const svec &w, double tol) const {
		svec rc(m), gc(n), gp(n), dc(n), tc(m);
		int i,k;

		rightmult(x,rc,w);
		for(i=0;i<m;i++) rc.x[i]-=b.x[i]*w.x[i];
		double rssb = rc.len2();
		leftmult(rc,gc,w);
		for(i=0;i<n;i++) gc.x[i] = -gc.x[i];
		double oldgm2 = 1000000;
		for(k=0;k<n;k++) {
			double gm2 = gc.len2();
			if (k==0) tol = gm2*tol;
			else if (gm2<tol) {
				//cout << "error small (" << k << "," << gm2 << ")" << endl;
				break;
			}
			if (gm2>oldgm2 && k>40) {
				//cout << "error increase ("<<k<<") ["<<gm2<<","<<oldgm2<<"]"<<endl;
				break;
			}
			oldgm2 = gm2;
			if (k>0) {
				double bi = gm2/gp.len2();
				for (i=0;i<n;i++) dc.x[i] = gc.x[i]+bi*dc.x[i];
			} else {
				for (i=0;i<n;i++) dc.x[i] = gc.x[i];
			}
			rightmult(dc,tc,w);
			double ai=gm2/tc.len2();
			for(i=0;i<n;i++) x.x[i] += ai*dc.x[i];
			for(i=0;i<m;i++) rc.x[i] += ai*tc.x[i];
			for(i=0;i<n;i++) gp.x[i]=gc.x[i];
			leftmult(rc,gc,w);
			for(i=0;i<n;i++) gc.x[i] = -gc.x[i];
		}
		//if (k==n) cout << "error BIG" << endl;
	}

	inline int getm(void) { return m; }
	inline int getn(void) { return n; }

	inline bool isfinite() {
		for(int i=0;i<m;i++)
			for(int j=0;j<in[i].length();j++)
				if (!finite(x[i][j])) return false;
		return true;
	}

private:
	int n,m;
	ilist<double> *x;
	ilist<int> *in;
};

#endif
