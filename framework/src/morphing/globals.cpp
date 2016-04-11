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
#include <containers/geohash.h>
#include <geometry/vertex.h>
#include <geometry/face.h>
#include <geometry/edge.h>

#include <globaldef.h>

//const FP hugenum = (FP)HUGE_VAL;
template<>
short geohash<vertex *>::globallimit = 500;
template<>
short geohash<face *>::globallimit = 500;
template<>
int geohash<vertex *>::levellimit = 8;
template<>
int geohash<face *>::levellimit = 20;
template<>
int geohash<vertex *>::pushlimit = 0;
template<>
int geohash<face *>::pushlimit = 1;
template<> 
FP geohash<vertex *>::overlap = (float)0.0;
template<>
FP geohash<face *>::overlap = (float)0.2;
template<>
int *geohash<vertex *>::counts = new int[geohash<vertex *>::levellimit+1];
template<>
int *geohash<face *>::counts = new int[geohash<face *>::levellimit+1];
template<>
int geohash<vertex *>::totalnum = -1;
template<>
int geohash<face *>::totalnum = -1;
template<>
int geohash<vertex *>::numq = 0;
template<>
int geohash<face *>::numq = 0;
template<>
int geohash<face *>::nums = 0;
template<>
int geohash<vertex *>::nums = 0;

