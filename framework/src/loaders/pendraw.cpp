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
#include <loaders/pendraw.h>
#include <geometry/edge.h>


static void initps(ostream &os) {
	os << "%!PS-Adobe-3.0 EPSF-3.0" << endl 
	   << "%%BoundingBox: -2.0 -2.0 502.0 502.0" << endl
	   << "gsave" << endl << "newpath" << endl
	   << "250.0 250.0 scale" << endl
	   << "1.0 1.0 translate" << endl
	   << "/dl {" << endl
	   << "  moveto lineto 0 0 0 setrgbcolor stroke" << endl
	   << "} def" << endl
	   << "%%EndProlog" << endl << endl
	   << "%%Page: 1 1" << endl
	   << "gsave" << endl
	   << "0 setlinewidth" << endl;
}

static void drawline(ostream &os, float x1, float y1, float x2, float y2) {
	os << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << " dl" << endl;
}

static void finishps(ostream &os) {
	os << "grestore" << endl << "grestore" << endl << "showpage" << endl;
}
	
void pendraw(trimesh *penmesh, ostream &os, FP scale) {

	int i,j;
	edge *e;
	int i1,i2;
	vec v1,v2;

	initps(os);
	for(i=0;i<penmesh->numedges();i++) {
		e = penmesh->getedge(i);
		if (!e->isend()) continue;
		i1 = e->vertexindex(0);
		i2 = e->vertexindex(1);
		if (i1==i2-1 || i1==i2+1 ||
		    (i1%2==0 && i2%2==0)) {
			v1 = penmesh->vertexpos(i1);
			v2 = penmesh->vertexpos(i2);
			drawline(os,v1[0]*scale, v1[1]*scale, v2[0]*scale, v2[1]*scale);
		}
	}
	finishps(os);
}
