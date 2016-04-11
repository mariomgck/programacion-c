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
#ifndef _CULL_H_
#define _CULL_H_

#include <globaldef.h>
#include <geometry/trirender.h>

class cull : public trirender {
public:
	cull(trirender *renderer,
	     int res=160, int thresh = 0,bool liberal=false, bool clip=false,
	     double minx=-HUGE_VAL, double maxx=HUGE_VAL,
	     double miny=-HUGE_VAL, double maxy=HUGE_VAL);
	virtual ~cull();

	virtual void start(double minx, double maxx,
		double miny, double maxy);
	virtual void rendertriangle(double x[3], double y[3], bool e[3],
		double cx[3], double cy[3],
		double r[3], double g[3], double b[3]);
	virtual void end();

	inline bool getliberal() { return LIBERAL; }
	inline void setliberal(bool liberal) { LIBERAL = liberal; }

private:
	unsigned char *image;
	char numtable[256];

	double maxxs,maxys,minxs,minys;
	
	double left,right,up,down;
	double boxleft,boxright,boxup,boxdown;
	bool clipbox;

	trirender *rend;

	class tri {
	public:
		tri *next, *prev;
		double x[3],y[3],cx[3],cy[3],r[3],g[3],b[3];
		bool e[3];
	};

	int RES,LIBERAL,THRESH,RES8;

	tri *head;

	void makenumtable();
	void drawhline(int x1, int x2, int y);
	int checkhline(int x1, int x2, int y);
	void findline(float x1, float y1, float x2, float y2,
		float &m, float &b);
	bool drawtriangle(tri *curr);
	void removetriangles();
};

#endif
