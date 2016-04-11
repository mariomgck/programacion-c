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
#ifndef _WARPPARMS_H_
#define _WARPPARMS_H_


#include <iostream>
#include <globaldef.h>
#include <warping/warperparams.h>

// These are the parameters used at only one level of the matching process.
// I guess this class should really be warpparams (with the extra a) but now
// it is a little late.
class warpparms  
{

public:
	// make one from a warperparams taking into account that ittnum iterations
	// have already passed (so alpha is less than it was in the beginning)
	warpparms(const warperparams &wp, int ittnum) {
		int i;
		for(alpha0=wp.alpha0,i=0;i<ittnum;i++,alpha0*=wp.eta)
			;
		for(zeta=wp.zeta,i=0;i<ittnum;i++,zeta*=wp.xi)
			;
		eta = wp.eta; rho = wp.rho;
		xi = wp.xi;
		epsilon = wp.epsilon;
		n = wp.n;
		momentum = wp.momentum;
		maxmove2 = wp.maxmove2;
		maxsp = wp.maxsp;
	}

	warpparms(FP Alpha0=0, FP Eta=0, FP Xi=0, FP Rho=0, FP Epsilon=0, FP N=0, int Maxsp=0,
		FP Momentum=0, FP Maxmove2=0) {
		alpha0 = Alpha0;
		eta = Eta;
		xi = Xi;
		rho = Rho;
		epsilon = Epsilon;
		n = N;
		maxsp = Maxsp;
		momentum = Momentum;
		maxmove2 = Maxmove2;
	}

	warpparms(istream &s) {
		s >> alpha0 >> eta >> xi >> rho >> epsilon >> momentum >> maxmove2 >> n >> maxsp;
	}

	~warpparms() {}

public:
	FP alpha0, eta, rho, epsilon, zeta; // from thesis
	FP xi; // eta:alpha :: xi:zeta (analogy)
	FP n;
	int maxsp;
	FP momentum, maxmove2; // momentum factor and maximum move of any point (squared)

};

#endif
