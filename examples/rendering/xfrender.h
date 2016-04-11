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
#ifndef _XFRENDER_H_
#define _XFRENDER_H_

#include <geometry/trirender.h>
#include <containers/ilist.h>
#include "forms.h"

// renders everything from (-1,-1) to (1,1) scaled to fit on the canvas
class xfrender : public trirender {
public:
	xfrender(FL_OBJECT *canvas, int red=255, int green=255, int blue=255,
		bool drawedges=false, int outlinered=0, int outlinegreen=0,
		int outlineblue=0);
	virtual ~xfrender();

	virtual void start(double minx, double maxx,
		double miny, double maxy);
	virtual void rendertriangle(double x[3], double y[3], bool e[3],
		double cx[3], double cy[3],
		double r[3], double g[3], double b[3]);
	virtual void end();

	void redraw();

	void changebackground(int red, int green, int blue);
	void changeoutline(bool drawedges, int outlinered, int outlinegreen,
		int outlintblue);

private:
	FL_OBJECT *c;

	class triangle {
	public:
		double x[3],y[3];
		bool e[3];
		int r,g,b;
	};
	int r,g,b;
	int olr,olg,olb;
	bool outline;

	ilist<triangle> t;

	void drawtriangle(const triangle &tr);
	void clearcanvas();
	
};

#endif
