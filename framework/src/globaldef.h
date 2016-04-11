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
#ifndef _GLOBALDEF_H_
#define _GLOBALDEF_H_

#include <stdlib.h>
#include <math.h>

// set FP to be whatever kind of floating point numbers you want to work with
// (they should be doubles really at this point)
typedef double FP;
//typedef float FP;

extern const double hugenum;

/*
typedef unsigned char bool;
#define true 1
#define false 0
*/

#define _finite(x) finite(x)
#define _isnan(x) isnan(x)
//#define HUGE_VAL 1000000

#endif
