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
#ifndef _WARPERPARAMS_H_
#define _WARPERPARAMS_H_

#include <iostream>
#include <globaldef.h>

using namespace std;

class warperparams  
{


public:
	// creates a new set of parameters from a file -- this is needed for the warper class
	warperparams(istream &s) {
		s >> *this;
	}
	
	warperparams() {
		momentum = maxmove2 = gamma = xi = alpha0 = eta = rho =
			epsilon = zeta = (FP)0.0;
		n = 0;
		t0 = t = maxsp = 0;
	}
	
	friend istream &operator>>(istream &is, warperparams &p);
	friend ostream &operator<<(ostream &os, const warperparams &p);

	~warperparams() {}

	inline FP colorscale() const { return gamma; }

public:
	FP alpha0, eta, rho, epsilon, zeta; // from thesis
	FP xi; // eta:alpha :: xi:zeta (analogy)
	FP n;
	int t0, t;
	int maxsp; // maximum number of springs from a vertex
	FP momentum, maxmove2; // momentum factor and maximum move of any point (squared)
	FP gamma;

};
inline istream &operator>>(istream &is, warperparams &p) {
	is >> p.gamma >> p.alpha0 >> p.zeta;
	is >> p.eta >> p.xi >> p.rho >> p.epsilon;
	is >> p.momentum >> p.maxmove2 >> p.n;
	is >> p.maxsp >> p.t0 >> p.t;
	return is;
}

inline ostream &operator<<(ostream &os, const warperparams &p) {
	os << p.gamma << ' ' << p.alpha0 << ' ' << p.zeta << endl;
	os << p.eta << ' ' << p.xi << ' ' << p.rho << ' ' << p.epsilon << endl;
	os << p.momentum << ' ' << p.maxmove2 << ' ' << p.n << endl;
	os << p.maxsp << ' ' << p.t0 << ' ' << p.t << endl;
	return os;
}

#endif
