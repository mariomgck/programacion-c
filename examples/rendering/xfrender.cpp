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
#include "xfrender.h"
#include <fstream>

using namespace std;

xfrender::xfrender(FL_OBJECT *canvas,int red,int green,int blue,
	bool drawedges,int outlinered,int outlinegreen,int outlineblue) {
	c = canvas;
	r = red;
	g = green;
	b = blue;
	outline = drawedges;
	olr = outlinered;
	olg = outlinegreen;
	olb = outlineblue;
}

xfrender::~xfrender() {
}

void xfrender::start(double minx, double maxx, double miny, double maxy) {
	t.setlength(0);
	//clearcanvas();
}

void xfrender::rendertriangle(double x[3], double y[3], bool e[3],
		double cx[3], double cy[3],
		double r[3], double g[3], double b[3]) {

	triangle tr;
	tr.r = tr.g = tr.b = 0;
	for(int i=0;i<3;i++) {
		tr.x[i] = x[i];
		tr.y[i] = y[i];
		tr.e[i] = e[i];	
		tr.r += (int)(r[i]*255);	
		tr.g += (int)(g[i]*255);
		tr.b += (int)(b[i]*255);
	}
	tr.r /= 3;
	tr.g /= 3;
	tr.b /= 3;
	if (tr.r<0) tr.r = 0; if (tr.r>255) tr.r = 255;
	if (tr.g<0) tr.g = 0; if (tr.g>255) tr.g = 255;
	if (tr.b<0) tr.b = 0; if (tr.b>255) tr.b = 255;
	//drawtriangle(tr);
	t += tr;
}

void xfrender::end() {
	fl_redraw_object(c);
}

void xfrender::changebackground(int red, int green, int blue) {
	r = red; g = green; b = blue;
	fl_redraw_object(c);
}

void xfrender::changeoutline(bool drawedges, int outlinered,
	int outlinegreen, int outlineblue) {

	outline = drawedges;
	olr = outlinered;
	olg = outlinegreen;
	olb = outlineblue;
	fl_redraw_object(c);
}

void xfrender::redraw() {
	clearcanvas();
	for(int i=0;i<t.length();i++)
		drawtriangle(t[i]);
}

void xfrender::drawtriangle(const triangle &tr) {
	double s = (c->w>c->h ? c->h/2 : c->w/2);
	int xoff,yoff;
	FL_POINT pts[3];
	
	if (c->w>c->h) {
		s = c->h/2;
		xoff = c->w/2-s;
		yoff = 0;
	} else {
		s = c->w/2;
		xoff = 0;
		yoff = c->h/2-s;
	}
	for(int i=0;i<3;i++) {
		pts[i].x = (short)((tr.x[i]+1)*s+c->x)+xoff;
		pts[i].y = (short)((1-tr.y[i])*s+c->y)+yoff;
	}
	fl_mapcolor(FL_FREE_COL1,tr.r,tr.g,tr.b);
	fl_polyf(pts,3,FL_FREE_COL1);
	if (outline) fl_mapcolor(FL_FREE_COL1,olr,olg,olb);
	for(int i=0;i<3;i++)
		if (tr.e[i]) fl_line(pts[i].x,pts[i].y,
			pts[(i+1)%3].x,pts[(i+1)%3].y,FL_FREE_COL1);
}

void xfrender::clearcanvas() {
	fl_mapcolor(FL_FREE_COL1,r,g,b);
	fl_rectf(c->x,c->y,c->w,c->h,FL_FREE_COL1);
}
