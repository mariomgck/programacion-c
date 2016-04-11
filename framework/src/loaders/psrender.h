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
#ifndef _PSRENDER_H_
#define _PSRENDER_H_

#include <geometry/trirender.h>
#include <fstream>

using namespace std;

class psrender : public trirender {
public:
	psrender(ostream *is, double scale = 1,
		double resolution=200,
		bool outline = false,
		bool allsolid = false,
		double outlinered = 0,
		double outlinegreen = 0,
		double outlineblue = 0);
	virtual ~psrender();

	virtual void start(double minx, double maxx,
		double miny, double maxy);
	virtual void rendertriangle(double x[3], double y[3], bool e[3],
		double cx[3], double cy[3],
		double r[3], double g[3], double b[3]);
	virtual void end();

private:
	ostream *s;
	double res;

	double olr,olg,olb;
	double sc;
	bool ol,solid;
	
	void drawtriangle(double x[3], double y[3],
		double r, double g, double b,
		bool e[3]);
	void drawgradstriangle(double x[3], double y[3],
		double r[3], double g[3], double b[3], 
		double px[3], double py[3], bool e[3]);
};

#endif
