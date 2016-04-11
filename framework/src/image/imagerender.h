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
#ifndef _IMAGERENDER_H_
#define _IMAGERENDER_H_

#include <geometry/trirender.h>
#include <image/image.h>
#include <float.h>

class imagerender : public trirender {
public:
	imagerender(int width=256, int height=256, FP backred=1,
	     FP backgreen=1, FP backblue=1, bool outline=false,
	     FP outlinered=0, FP outlinegreen=0, FP outlineblue = 0,
	     double minx=-HUGE_VAL, double maxx=HUGE_VAL,
	     double miny=-HUGE_VAL, double maxy=HUGE_VAL);
	virtual ~imagerender();

	virtual void start(double minx, double maxx,
		double miny, double maxy);
	virtual void rendertriangle(double x[3], double y[3], bool e[3],
		double cx[3], double cy[3],
		double r[3], double g[3], double b[3]);
	virtual void end();

	inline image getimage() { return im; }

private:
	image im;
	FP br,bg,bb;
	FP olr,olg,olb;

	double maxxs,maxys,minxs,minys;
	
	double left,right,up,down;
	double boxleft,boxright,boxup,boxdown;

	void drawhline(int x1, int x2, int y, FP r1, FP g1, FP b1,
		FP r2, FP g2, FP b2);
	void findline(float x1, float y1, float x2, float y2,
		float &m, float &b);
};

#endif
