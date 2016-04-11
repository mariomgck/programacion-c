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
#ifndef _VEC_H_
#define _VEC_H_

#include <math.h>
#include <iostream>
#include <fstream>
#include "globaldef.h"
#include "float.h"

using namespace std;

// a large set of functions for doing vec and a few matrix operators
// for vectors in a fixed dimensional space

const int DIM = 6;

class vec  
{
public:
	inline vec() {}
	inline vec(const FP v) {
		for(int i=0;i<DIM;i++) x[i] = v;
	}
	inline vec(const int v) {
		for(int i=0;i<DIM;i++) x[i] = v;
	}
	inline vec(const FP * const v) {
		for(int i=0;i<DIM;i++) x[i] = v[i];
	}

	inline ~vec() {}

	inline vec(const vec& a) {
		for(int i=0;i<DIM;i++) x[i] = a.x[i];
	}

	inline vec operator-() const {
		vec ret;
		for(int i=0;i<DIM;i++) ret.x[i] = -x[i];
		return ret;
	}
	inline vec& operator=(FP a) {
		for(int i=0;i<DIM;i++) x[i] = a;
		return *this;
	}
	inline vec& operator+=(const vec& a) {
		for(int i=0;i<DIM;i++) x[i] += a.x[i];
		return *this;
	}
	inline vec& operator-=(const vec& a) {
		for(int i=0;i<DIM;i++) x[i] -= a.x[i];
		return *this;
	}
	inline vec& operator+=(FP a) {
		for(int i=0;i<DIM;i++) x[i] += a;
		return *this;
	}
	inline vec& operator-=(FP a) {
		for(int i=0;i<DIM;i++) x[i] -= a;
		return *this;
	}
	inline vec& operator*=(FP a) {
		for(int i=0;i<DIM;i++) x[i] *= a;
		return *this;
	}
	inline vec& operator/=(FP a) {
		for(int i=0;i<DIM;i++) x[i] /= a;
		return *this;
	}
	inline vec& operator/=(const vec& a) {
		for(int i=0;i<DIM;i++) x[i] /= a.x[i];
		return *this;
	}
	inline FP operator*(const vec& a) const {
		FP d(x[0]*a.x[0]);
		for(int i=1;i<DIM;i++) d+=x[i]*a.x[i];
		return d;
	}
	inline bool operator==(const vec& a) const {
		for(int i=0;i<DIM;i++) if (x[i]!=a.x[i]) return false;
		return true;
	}
	inline bool operator!=(const vec& a) const {
		for(int i=0;i<DIM;i++) if (x[i]!=a.x[i]) return true;
		return false;
	}
	inline bool operator==(FP a) const {
		for(int i=0;i<DIM;i++) if (x[i]!=a) return false;
		return true;
	}
	inline bool operator!=(FP a) const {
		for(int i=0;i<DIM;i++) if (x[i]!=a) return true;
		return false;
	}
	inline const vec operator/ (FP b) const {
		return vec(*this) /= b;
	}

	inline bool isvalid() const {
		for(int i=0;i<DIM;i++) if (!_finite(x[i])) return false;
		return true;
	}
	inline bool visnan() const {
		for(int i=0;i<DIM;i++) if (_isnan(x[i])) return true;
		return false;
	}
	
	inline FP len2() const { FP d=(x[0]*x[0]); for(int i=1;i<DIM;i++) d+=x[i]*x[i]; return d; }
	inline FP len() const { return sqrt(len2()); }
	inline vec norm() const { return (*this)/len(); }

	inline FP max() const { FP m=x[0];
		for(int i=1;i<DIM;i++) if (m<x[i]) m=x[i];
		return m;
	}
	inline FP min() const { FP m=x[0];
		for(int i=1;i<DIM;i++) if (m>x[i]) m=x[i];
		return m;
	}

	inline FP absmax() const { FP m=x[0]<0 ? -x[0] : x[0];
		for(int i=1;i<DIM;i++)
			if (x[i]<0) {
				if (-x[i]>m) m = -x[i];
			} else {
				if (x[i]>m) m = x[i];
			}
		return m;
	}

	inline FP absmin() const { FP m=x[0]<0 ? -x[0] : x[0];
		for(int i=1;i<DIM;i++)
			if (x[i]<0) {
				if (-x[i]<m) m=-x[i];
			} else {
				if (x[i]<m) m=x[i];
			}
		return m;
	}

	inline void copyselect(vec &to,int start,int end) const {
		int i;
		for(i=0;i<start;i++) to.x[i] = 0;
		for(;i<=end;i++) to.x[i] = x[i];
		for(;i<DIM;i++) to.x[i] = 0;
	}

	// these two are for writable and non-writable access to the vector
	// some may find the two confusing, but they are necessary since there
	// is no way of getting the compiler to use the proper one via over-loading
	// (I don't want to use proxy classes just for this!)
	inline FP operator()(int i) const { return x[i]; }
	inline FP &operator[](int i)  { return x[i]; }

	friend ostream& operator<<(ostream& s, const vec& a);
	friend istream& operator>>(istream& s, vec& a);

	const static vec ones;
	const static vec zero;

private:
	FP x[DIM];

};


inline ostream& operator<<(ostream& s, const vec& a) {

	for(int i=0;i<DIM;i++) { s << a.x[i]; if (i!=DIM) s << ' '; }
	return s;
}

inline istream& operator>>(istream& s, vec& a) {
	for(int i=0;i<DIM;i++) { s >> a.x[i]; }
	return s;
}

inline const vec operator+ (const vec& a, const vec& b) {
	return vec(a) += b;
}

inline const vec operator- (const vec& a, const vec& b) {
	return vec(a) -= b;
}

inline const vec operator+ (const vec& a, FP b) {
	return vec(a) += b;
}

inline const vec operator+ (FP a, const vec& b) {
	return vec(b) += a;
}

inline const vec operator- (const vec& a, FP b) {
	return vec(a) -= b;
}

inline const vec operator* (const vec& a, FP b) {
	return vec(a) *= b;
}

inline const vec operator* (FP a, const vec& b) {
	return vec(b) *= a;
}

inline const vec operator/ (const vec& a, const vec &b) {
	return vec(a) /= b;
}

inline vec min(const vec &a, const vec &b) {
	vec r;
	for(int i=0;i<DIM;i++)
		r[i] = a(i)<b(i) ? a(i) : b(i);
	return r;
}

inline vec max(const vec &a, const vec &b) {
	vec r;
	for(int i=0;i<DIM;i++)
		r[i] = a(i)>b(i) ? a(i) : b(i);
	return r;
}

inline FP planeangle (const vec& t1, const vec& t2, const vec& s1, const vec& s2) {

	FP d,r;

	d = t1*s1;
	r = d*d;
	d = t1*s2;
	r += d*d;
	d = t2*s1;
	r += d*d;
	d = t2*s2;
	r += d*d;
	return 2-r;
}

inline vec lineintersection(const vec& p1, const vec& d1, const vec& p2, const vec& d2) {
	/*vec n,k,l;
	FP d,c;

	n = d1.norm();
	k = (d2-(n*d2)*n).norm();
	d = k*d2;
	if (d!=0) {
		l = p1-p2;
		c = k*l;
		return p2+d2*c/d;
	} 
	return p2;*/
	FP d,D,dd1;
	vec l;
	l=p2-p1;
	D = d1*d2;
	dd1 = d1*d1;
	d = dd1*(d2*d2) - D*D;
	if (d==0) return p2;
	d = (dd1*(d2*l) - D*(d1*l))/d;
	return p2-d*d2;
}



class matrix {
public:
	inline matrix() {}
	inline ~matrix() {}
	inline matrix(const matrix &m) {
		for(int i=0;i<DIM;i++) x[i] = m.x[i];
	}
	inline matrix(const FP &a,bool diaonly=false) {
		if (diaonly) {
			for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j] = 0;
			for(int i=0;i<DIM;i++) x[i][i] = a;
		} else {
			for(int i=0;i<DIM;i++) x[i] = a;
		}
	}
	inline matrix(FP v[DIM][DIM]) {
		for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j] = v[i][j];
	}
	inline matrix(const vec &v1, const vec &v2) {  // forms a matrix of the outer product
		for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) x[i][j] = v1(i)*v2(j);
	}

	inline matrix operator-() const {
		matrix ret;
		for(int i=0;i<DIM;i++) ret.x[i] = -x[i];
		return ret;
	}
	inline matrix& operator=(FP a) {
		for(int i=0;i<DIM;i++) x[i] = a;
		return *this;
	}


	inline bool isvalid() const {
		for(int i=0;i<DIM;i++) if(!x[i].isvalid()) return false;
		return true;
	}

	inline matrix operator*(const matrix &m) const {
		matrix ret;
		for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++) {
			ret.x[i][j] = 0;
			for(int k=0;k<DIM;k++)
				ret.x[i][j] += x[i](k) * m.x[k](j);
		}
		return ret;
	}
	inline vec operator*(const vec &v) const {
		vec ret;
		for(int i=0;i<DIM;i++) ret[i] = x[i]*v;
		return ret;
	}

	inline vec &operator[](int r) {
		return x[r];
	}
	inline vec operator[](int r) const {
		return x[r];
	}

	inline FP operator()(int m, int n) const {
		return x[m](n);
	}
	inline FP &operator()(int m, int n) {
		return x[m][n];
	}

	inline matrix t() const {
		matrix ret;
		for(int i=0;i<DIM;i++) for(int j=0;j<DIM;j++)
			ret.x[i][j] = x[j](i);
		return ret;
	}

	inline matrix &operator+=(const matrix &m) {
		for(int i=0;i<DIM;i++) x[i] += m.x[i];
		return *this;
	}

	inline matrix &operator-=(const matrix &m) {
		for(int i=0;i<DIM;i++) x[i] -= m.x[i];
		return *this;
	}

	inline matrix &operator*=(const matrix &m) {
		matrix a = (*this)*m;
		*this = a;
		return *this;
	}

	inline matrix &operator+=(const FP &a) {
		for(int i=0;i<DIM;i++) x[i] += a;
		return *this;
	}

	inline matrix &operator-=(const FP &a) {
		for(int i=0;i<DIM;i++) x[i] -= a;
		return *this;
	}

	inline matrix &operator*=(const FP &a) {
		for(int i=0;i<DIM;i++) x[i] *= a;
		return *this;
	}

	inline matrix &operator/=(const FP &a) {
		for(int i=0;i<DIM;i++) x[i] /= a;
		return *this;
	}

	inline FP max() const {
		FP t,m = x[0].max();
		for(int i=1;i<DIM;i++) 
			if ((t=x[i].max())>m) m=t;
		return m;
	}

	inline FP min() const {
		FP t,m = x[0].min();
		for(int i=1;i<DIM;i++) 
			if ((t=x[i].min())<m) m=t;
		return m;
	}

	inline FP absmin() const {
		FP t,m = x[0].absmin();
		for(int i=1;i<DIM;i++)
			if ((t=x[i].absmin())<m) m=t;
		return m;
	}

	inline FP absmax() const {
		FP t,m = x[0].absmax();
		for(int i=1;i<DIM;i++)
			if ((t=x[i].absmax())>m) m=t;
		return m;
	}

	inline bool operator==(const matrix &m) const {
		for(int i=0;i<DIM;i++) if (m.x[i]!=x[i]) return false;
		return true;
	}
	inline bool operator==(const FP &a) const {
		for(int i=0;i<DIM;i++) if (x[i]!=a) return false;
		return true;
	}
	inline bool operator!=(const matrix &m) const {
		for(int i=0;i<DIM;i++) if (m.x[i]!=x[i]) return true;
		return false;
	}
	inline bool operator!=(const FP &a) const {
		for(int i=0;i<DIM;i++) if (x[i]!=a) return true;
		return false;
	}

	inline FP norm2() const { // returns the square of the frobenius norm
		FP ret=x[0].len2();
		for(int i=1;i<DIM;i++) ret += x[i].len2();
		return ret;
	}
	inline FP norm() const { // returns the frobenius norm
		return sqrt(norm2());
	}

	friend ostream& operator<<(ostream& s, const matrix& a);
	friend istream& operator>>(istream& s, matrix& a);

	// LU decomposition -- ix is the row permutations
	int LUdecomp(matrix &LU, int *ix) const;
	// LU back substitution -- ix from above fn call (this should be an LU combination)
	void LUbacksub(int *ix, vec &col) const;

	// solves equation Ax=b (A is this, x is the returned value)
	vec solve(const vec &b, bool &worked) const;
	inline vec solve(const vec &b) const { bool w; return solve(b,w); }

	matrix inv(bool &worked) const;
	inline matrix inv() const { bool w; return inv(w); }

	const static matrix eye;
	const static matrix zero;

	inline ostream &niceprint(ostream& s) {
		for(int i=0;i<DIM;i++) {
			for(int j=0;j<DIM;j++)
				s << x[j][i] << ' ';
			s << endl;
		}
		return s;
	}


private:
	// these are the *rows* of the matrix
	vec x[DIM];
};

inline matrix operator+(const matrix &a, const matrix &b) {
	return matrix(a)+=b;
}
inline matrix operator-(const matrix &a, const matrix &b) {
	return matrix(a)-=b;
}
inline matrix operator+(const matrix &a, const FP &b) {
	return matrix(a)+=b;
}
inline matrix operator-(const matrix &a, const FP &b) {
	return matrix(a)-=b;
}
inline matrix operator+(const FP &a, const matrix &b) {
	return matrix(b)+=a;
}
inline matrix operator-(const FP &a, const matrix &b) {
	return matrix(a)-=b;
}
inline matrix operator*(const matrix &a, const FP &b) {
	return matrix(a)*=b;
}
inline matrix operator*(const FP &b, const matrix &a) {
	return matrix(a)*=b;
}
inline matrix operator/(const matrix &a, const FP &b) {
	return matrix(a)/=b;
}

inline ostream& operator<<(ostream& s, const matrix& m) {
	for(int i=0;i<DIM;i++) { s << m.x[i]; if (i!=DIM) s << ' '; }
	return s;
}

inline istream& operator>>(istream& s, matrix& m) {
	for(int i=0;i<DIM;i++) { s >> m.x[i]; }
	return s;
}

#endif
