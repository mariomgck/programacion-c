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
#ifndef _GMATRIX_H_
#define _GMATRIX_H_

#include <math.h>
//#include <nan.h>
#include <iostream>
#include <fstream>
#include "globaldef.h"
using namespace std;

class gmatrix {
public:
	inline gmatrix(int m=1, int n=1) {
		this->m = m; this->n = n;
		s = m*n;
		x = new FP[s];
	}
	inline ~gmatrix() { delete []x; }
	inline gmatrix(const gmatrix &ma, bool transpose=false) {
		s = ma.m*ma.n;
		x = new FP[s];
		if (transpose) {
			int i,j,c;
			n = ma.m; m = ma.n;
			c = 0;
			for(i=0;i<m;i++) for(j=0;j<n;j++,c++)
				x[c] = ma.x[i+j*m];
		} else {
			n = ma.n; m = ma.m;
			int i;
			for(i=0;i<s;i++) x[i] = ma.x[i];
		}
	}
	inline gmatrix(int m, int n,const FP &a,bool diaonly=false) {
		this->m = m;
		this->n = n;
		s = m*n;
		x = new FP[s];
		if (diaonly) {
			int i,j=m*m;
			for(i=0;i<j;i++) x[i] = 0;
			if (m>=m)
				for(i=0;i<m;i++) x[i*n+i] = a;
			else for(i=0;i<n;i++) x[i*n+i] = a;
		} else {
			int i,j=m*n;
			for(i=0;i<j;i++) x[i] = a;
		}
	}
	inline gmatrix(FP *v,int m, int n) {
		this->m = m;
		this->n = n;
		s = m*n;
		int i;
		x = new FP[s];
		for(i=0;i<s;i++) x[i] = v[i];
	}
	// forms a matrix of the outer product (ie v1*v2') -- v2 is
	// "transposed" temporarily for this operation
	inline gmatrix(const gmatrix &v1, const gmatrix &v2) {
		if (v1.n!=v2.n) {
			s = 1;
			m=1; n=1; x = new FP[1];
			//x[0] = NAN;
			x[0] = 0.0/0.0;
		} else {
			n = v2.m; m = v1.m;
			s = n*m;
			x = new FP[s];
			int i,j,k,c=0;
			for(i=0;i<m;i++) for(j=0;j<n;j++,c++) {
				x[c] = 0;
				for(k=0;k<v1.n;k++)
					x[c] += v1.x[i*v1.n+k]*
						v2.x[j*v1.n+k];
			}
		}
	}

	inline gmatrix operator-() const {
		gmatrix ret(m,n);
		for(int i=0;i<s;i++) ret.x[i] = -x[i];
		return ret;
	}
	inline gmatrix& operator=(FP a) {
		for(int i=0;i<s;i++) x[i] = a;
		return *this;
	}
	inline gmatrix& operator=(const gmatrix &ma) {
		if (&ma==this) return *this;
		if (ma.n != n || ma.m != m) {
			s = ma.s; m = ma.m; n = ma.n;
			delete []x;
			x = new FP[s];
		}
		for(int i=0;i<s;i++) x[i] = ma.x[i];
		return *this;
	}

	inline bool isvalid() const {
		for(int i=0;i<s;i++) if(!_finite(x[i])) return false;
		return true;
	}

	inline gmatrix operator*(const gmatrix &ma) const {
		if (n!=ma.m) {
			cerr << "invalid matrix multiply" << endl;
			exit(1);
		}
		gmatrix ret(m,ma.n);
		int c=0;
		for(int i=0;i<m;i++) for(int j=0;j<ma.n;j++,c++) {
			ret.x[c] = 0;
			for(int k=0;k<n;k++)
				ret.x[c] += x[i*n+k] * ma.x[k*ma.n+j];
		}
		return ret;
	}
	inline void outer(const gmatrix &ma, gmatrix &ret) const {
		if (n!=ma.n || ret.m!=m || ret.n!=ma.m) {
			cerr << "invalid matrix outer multiply" << endl;
			exit(1);
		}
		int c=0;
		for(int i=0;i<m;i++) for(int j=0;j<ma.m;j++,c++) {
			ret.x[c] = 0;
			for(int k=0;k<n;k++)
				ret.x[c] += x[i*n+k] * ma.x[j*ma.m+k];
		}
	}
	inline void inner(const gmatrix &ma, gmatrix &ret) const {
		if (m!=ma.m || ret.m!=n || ret.n!=ma.n) {
			cerr << "invalid matrix inner multiply" << endl;
			exit(1);
		}
		int c=0;
		for(int i=0;i<n;i++) for(int j=0;j<ma.n;j++,c++) {
			ret.x[c] = 0;
			for(int k=0;k<m;k++)
				ret.x[c] += x[k*m+i] * ma.x[k*m+j];
		}
	}
		

	inline FP operator()(int i, int j) const {
		return x[i*n+j];
	}

	inline FP *operator[](int i) {
		return x+(i*n);
	}
	inline gmatrix t() const {
		return gmatrix(*this,true);
	}

	inline gmatrix &operator+=(const gmatrix &ma) {
		if (m!=ma.m||n!=ma.n) {
			cerr << "invalid matrix addition" << endl;
			exit(1);
		}
		for(int i=0;i<s;i++) x[i] += ma.x[i];
		return *this;
	}

	inline gmatrix &operator-=(const gmatrix &ma) {
		if (m!=ma.m||n!=ma.n) {
			cerr << "invalid gmatrix addition" << endl;
			exit(1);
		}
		for(int i=0;i<s;i++) x[i] -= ma.x[i];
		return *this;
	}

	inline gmatrix &operator*=(const gmatrix &ma) {
		gmatrix a = (*this)*ma;
		*this = a;
		return *this;
	}

	inline gmatrix &operator+=(const FP &a) {
		for(int i=0;i<s;i++) x[i] += a;
		return *this;
	}

	inline gmatrix &operator-=(const FP &a) {
		for(int i=0;i<s;i++) x[i] -= a;
		return *this;
	}

	inline gmatrix &operator*=(const FP &a) {
		for(int i=0;i<s;i++) x[i] *= a;
		return *this;
	}

	inline gmatrix &operator/=(const FP &a) {
		for(int i=0;i<s;i++) x[i] /= a;
		return *this;
	}

	inline FP max() const {
		FP t,ma = x[0];
		for(int i=1;i<s;i++) 
			if ((t=x[i])>ma) ma=t;
		return ma;
	}

	inline FP min() const {
		FP t,mi = x[0];
		for(int i=1;i<s;i++) 
			if ((t=x[i])<mi) mi=t;
		return mi;
	}

	inline FP absmin() const {
		FP t,mi = x[0]>0?x[0]:-x[0];
		for(int i=1;i<s;i++)
			if ((t=(x[i]>0?x[i]:-x[i]))<mi) mi=t;
		return mi;
	}

	inline FP absmax() const {
		FP t,ma = x[0]>0?x[0]:-x[0];
		for(int i=1;i<s;i++)
			if ((t=(x[i]>0?x[i]:-x[i]))>ma) ma=t;
		return ma;
	}

	inline bool operator==(const gmatrix &ma) const {
		if (ma.n!=n||ma.m!=m) return false;
		for(int i=0;i<s;i++) if (ma.x[i]!=x[i]) return false;
		return true;
	}
	inline bool operator==(const FP &a) const {
		for(int i=0;i<s;i++) if (x[i]!=a) return false;
		return true;
	}
	inline bool operator!=(const gmatrix &ma) const {
		if (ma.n!=n||ma.m!=m) return false;
		for(int i=0;i<s;i++) if (ma.x[i]!=x[i]) return true;
		return false;
	}
	inline bool operator!=(const FP &a) const {
		for(int i=0;i<s;i++) if (x[i]!=a) return true;
		return false;
	}

	inline FP norm2() const { // returns the square of the frobenius norm
		FP ret=x[0]*x[0];
		for(int i=1;i<s;i++) ret += x[i]*x[i];
		return ret;
	}
	inline FP norm() const { // returns the frobenius norm
		return sqrt(norm2());
	}

	friend ostream& operator<<(ostream& s, const gmatrix& ma);
	friend istream& operator>>(istream& s, gmatrix& ma);

	// LU decomposition -- ix is the row permutations
	int LUdecomp(gmatrix &LU, int *ix) const;
	// LU back substitution --
	//    ix from above fn call (this should be an LU combination)
	void LUbacksub(int *ix, FP *col) const;

	// solves equation Ax=b (A is this, x is the returned value)
	FP *solve(const FP *b, bool &worked) const;
	inline FP *solve(const FP *b) const { bool w; return solve(b,w); }

	gmatrix inv(bool &worked) const;
	inline gmatrix inv() const { bool w; return inv(w); }

	// w needs to points to an array of n fp numbers
	void svd(gmatrix &u, gmatrix &v, FP *w);

	inline ostream &niceprint(ostream& s) {
		for(int i=0;i<m;i++) {
			for(int j=0;j<n;j++)
				s << x[i*n+j] << ' ';
			s << endl;
		}
		return s;
	}

	inline int getm() const { return m; }
	inline int getn() const { return n; }

private:
	static FP pythag(FP a, FP b);

	int m,n,s;
	FP *x;
};

inline gmatrix operator+(const gmatrix &a, const gmatrix &b) {
	return gmatrix(a)+=b;
}
inline gmatrix operator-(const gmatrix &a, const gmatrix &b) {
	return gmatrix(a)-=b;
}
inline gmatrix operator+(const gmatrix &a, const FP &b) {
	return gmatrix(a)+=b;
}
inline gmatrix operator-(const gmatrix &a, const FP &b) {
	return gmatrix(a)-=b;
}
inline gmatrix operator+(const FP &a, const gmatrix &b) {
	return gmatrix(b)+=a;
}
inline gmatrix operator-(const FP &a, const gmatrix &b) {
	return gmatrix(b.getn(),b.getm(),a)-=b;
}
inline gmatrix operator*(const gmatrix &a, const FP &b) {
	return gmatrix(a)*=b;
}
inline gmatrix operator*(const FP &b, const gmatrix &a) {
	return gmatrix(a)*=b;
}
inline gmatrix operator/(const gmatrix &a, const FP &b) {
	return gmatrix(a)/=b;
}

inline ostream& operator<<(ostream& s, const gmatrix& ma) {
	s << ma.m << ' ' << ma.n << ' ';
	for(int i=0;i<ma.s;i++) { s << ma.x[i]; if (i!=ma.s) s << ' '; }
	return s;
}

inline istream& operator>>(istream& s, gmatrix& ma) {
	int tn,tm;
	s >> tm >> tn;
	if (tm!=ma.m || tn!=ma.n) {
		delete []ma.x;
		ma.s = tm*tn;
		ma.x = new FP[ma.s];
		ma.m = tm; ma.n = tn;
	}
	for(int i=0;i<ma.s;i++) { s >> ma.x[i]; }
	return s;
}
#endif
