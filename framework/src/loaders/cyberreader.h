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
#ifndef _CYBERREADER_H_
#define _CYBERREADER_H_

#include <iostream>
#include "globaldef.h"

#include <matrix/vec.h>
#include <image/image.h>

class cyberreader  
{
public:
	// scale is applied to spatial coordinates only.  gamma is applied to color values only
	// it is probably just best to leave them alone and adjust the color:shape ratio when 
	// loading up the mesh itself
	// Note that the streams need to be opens with ios::binary as the second argument!
	cyberreader(istream &cyberfile, istream &ppmfile, 
		FP scale = (FP)0.00001, FP gamma=1.0);

	// This loads a cyberscan with no color information
	cyberreader(istream &cyberfile, FP scale = (FP)0.00001);
	~cyberreader();

	// returns the scan size (height= number of latitidues, around=number of longitudes)
	inline void scansize(int &height, int &around) const { height = numlat; around = numlon; }

	// returns false if the point could not be found, otherwise leaves the 6D point
	// in ret
	inline bool point(int h, int theta, vec &ret) const { 
		int index = theta*numlat + h;

		if (radius[index]==hugenum) return false;
		FP th = theta*thetamult;
		
		ret[2] = -cos(th)*radius[index];
		ret[0] = sin(th)*radius[index];
		ret[1] = hmult*(h-numlat/2);
		ret[3] = r[index];
		ret[4] = g[index];
		ret[5] = b[index];
		return true;
	}
	// same as above but the index has already been computed (actually, this
	//  doesn't save any computational time, it just might be more convienient)
	inline bool point(int index, vec &ret) const {
		if (radius[index]==hugenum) return false;
		int theta = index/numlat;
		int h = index%numlat;
		FP th = theta*thetamult;

		ret[2] = -cos(th)*radius[index];
		ret[0] = sin(th)*radius[index];
		ret[1] = hmult*(h-numlat/2);
		ret[3] = r[index];
		ret[4] = g[index];
		ret[5] = b[index];
		return true;
	}
	// takes a 6D point and finds its back projection (it actually only cares
	//  about the spatial coordinates)  It assumes that the point is from
	//  the cyberware scan.
	inline void reversepoint(const vec &v, int &h, int &theta) {
		theta = (int)floor((atan2(v(0),v(2))+M_PI)/thetamult + 0.5);
		if (theta<0) theta = 0;
		if (theta>=numlon) theta = numlon-1;
		h = (int)floor(v(1)/hmult + numlat/2 + 0.5);
		if (h<0) h = 0;
		if (h>=numlat) h = numlat-1;
	}
	// just returns the radius at a point
	inline bool getradius(int h, int theta, FP &ret) const {
		int index = theta*numlat+h;
		if (radius[index]==hugenum) return false;
		ret = radius[index];
		return true;
	}
 
	// maximum radius
	inline FP maxr() const { return mr; }
	// total height of the scan -- this is the height of the cylinder
	inline FP height() const { return numlat*hmult; }
	// returns the gamma value passed in in the beginning
	inline FP getgamma() const { return gamma; }
	// this and the next function assume that all of the radii are of
	// the maximum radius (thus producing a bonding cylinder) and then
	// bound this cylinder with an axis aligned bounding box.  minpt 
	// returns the corner of the box with minimum coordinates and maxpt
	// returns the corner of the box with maximum coordinates.
	inline vec minpt() const { 
		vec v;
		v[3]=v[4]=v[5] = 0;
		v[0]=v[1]=-mr;
		v[2]=-numlat*hmult/2;
		return v;
	}
	inline vec maxpt() const {
		vec v;
		v[3]=v[4]=v[5] = gamma;
		v[0]=v[1]=mr;
		v[2]=numlat*hmult/2;
		return v;
	}

	void clip(int minheight, int maxheight);

	// the following two functions take into account the warp around nature of the image
	void reduce(); // reduces the data to the largest connected component
	void fillholes(); // fills in any holes in the data set

	// this returns an image (which cannot outlast the cyberreader) of 
	// the scan (the image planes are: red, green, blue, radius)
	// with the scale and gamma applied.
	image getimage() { 
		FP **t = new FP *[4]; t[0] = r; t[1] = g; t[2] = b; t[3] = radius;
		image ret(numlat,numlon,t,4,false,true); delete []t; return ret;
	}
		
private:
	void readcyberheader(istream &s, FP scale);
	void readcyberdata(istream &s, FP scale);
	void readnewcyberheader(istream &s, FP scale);

	int *findholes();

	// selected data from the cyberware header
	int numlat,numlon;
	int rshift;
	int totalheight;
	FP thetamult,hmult;
	FP mr;

	FP gamma,scale;

	// these have already had gamma and scale applied
	FP *radius;
	FP *r,*g,*b;
};

#endif
