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
#include <math.h>
#include <image/imagerender.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

imagerender::imagerender(int width, int height, FP backred,
             FP backgreen, FP backblue, bool outline,
             FP outlinered, FP outlinegreen, FP outlineblue,
             double minx, double maxx, double miny, double maxy) :
	im(width,height) {

	im.addchannel();
	im.addchannel();
	im.addchannel();
	olr = outlinered;
	olg = outlinegreen;
	olb = outlineblue;
	br = backred;
	bg = backgreen;
	bb = backblue;
	up = miny;
	down = maxy;
	left = minx;
	right = maxx;
}

imagerender::~imagerender() {
}

void imagerender::start(double minx, double maxx, double miny, double maxy) {

	im.set(br,0);
	im.set(bg,1);
	im.set(bb,2);

	minxs = minx;
	minys = miny;
	maxxs = maxx;
	maxys = maxy;
	if (up==-HUGE_VAL) boxup = miny;
	else boxup = up;
	if (down==HUGE_VAL) boxdown = maxy;
	else boxdown = down;
	if (left==-HUGE_VAL) boxleft = minx;
	else boxleft = left;
	if (right==HUGE_VAL) boxright = maxx;
	else boxright = right;
}


void imagerender::end() {
}

#define clipcolor(c) (c>0.0 ? (c<1.0 ? c : 1.0) : 0.0)

void imagerender::drawhline(int x1, int x2, int y, FP r1, FP g1, FP b1,
	FP r2, FP g2, FP b2) {


	y = im.getheight()-y-1;
	if (x1>x2) { return; }
	// {int x = x1; x1 = x2; x2 = x; }
	if (y<0 || y>=im.getheight() || x1>=im.getwidth() || x2<0) return;
	if (x1<0) x1 = 0;
	if (x2>=im.getwidth()) x2 = im.getwidth()-1;

	int x,i;
	FP dr,dg,db;
	dr = (r2-r1)/(x2-x1);
	dg = (g2-g1)/(x2-x1);
	db = (b2-b1)/(x2-x1);
	for(i=x1+y*im.getwidth(),x=x1;x<=x2;x++,i++) {
		im.p(i,0) = clipcolor(r1);
		im.p(i,1) = clipcolor(g1);
		im.p(i,2) = clipcolor(b1);
		r1 += dr; g1 += dg; b1 += db;
	}
}


#define min3(a,b,c) ((a<b)?((a<c)?a:c):((b<c)?b:c))
#define max3(a,b,c) ((a>b)?((a>c)?a:c):((b>c)?b:c))

void imagerender::findline(float x1, float y1, float x2, float y2, 
		float &m, float &b) {

	m = (x1-x2)/(y1-y2);
	b = x1 - m*y1;
}

void imagerender::rendertriangle(double x[3], double y[3], bool e[3],
                double cx[3], double cy[3],
                double r[3], double g[3], double b[3]) {

	float topx,topy,leftx,lefty,rightx,righty;
	float rm,lm,rk,lk,plb,prb,clb,crb,trb,tlb;
	FP ra,rb,rc,ga,gb,gc,ba,bb,bc,den;
	int b1,b2;
	int lyf,ryf,fl=0;
	int ret = 0;
	
	den = -cx[1]*cy[2]+cx[0]*cy[2]-cx[0]*cy[1]-
		cx[2]*cy[0]+cx[2]*cy[1]+cx[1]*cy[0];

	ra = (-cy[0]*r[2]+cy[1]*r[2]-cy[1]*r[0]-
		r[1]*cy[2]+r[1]*cy[0]+r[0]*cy[2])/den;
	rb = (-cx[0]*r[1]-cx[1]*r[2]+cx[2]*r[1]+
		cx[1]*r[0]+cx[0]*r[2]-cx[2]*r[0])/den;
	rc = (cy[0]*cx[1]*r[2]-r[0]*cx[1]*cy[2]-cy[1]*cx[0]*r[2]+
		cy[1]*cx[2]*r[0]+r[1]*cx[0]*cy[2]-cy[0]*cx[2]*r[1])/den;

	ga = (-cy[0]*g[2]+cy[1]*g[2]-cy[1]*g[0]-
		g[1]*cy[2]+g[1]*cy[0]+g[0]*cy[2])/den;
	gb = (-cx[0]*g[1]-cx[1]*g[2]+cx[2]*g[1]+
		cx[1]*g[0]+cx[0]*g[2]-cx[2]*g[0])/den;
	gc = (cy[0]*cx[1]*g[2]-g[0]*cx[1]*cy[2]-cy[1]*cx[0]*g[2]+
		cy[1]*cx[2]*g[0]+g[1]*cx[0]*cy[2]-cy[0]*cx[2]*g[1])/den;

	ba = (-cy[0]*b[2]+cy[1]*b[2]-cy[1]*b[0]-
		b[1]*cy[2]+b[1]*cy[0]+b[0]*cy[2])/den;
	bb = (-cx[0]*b[1]-cx[1]*b[2]+cx[2]*b[1]+
		cx[1]*b[0]+cx[0]*b[2]-cx[2]*b[0])/den;
	bc = (cy[0]*cx[1]*b[2]-b[0]*cx[1]*cy[2]-cy[1]*cx[0]*b[2]+
		cy[1]*cx[2]*b[0]+b[1]*cx[0]*cy[2]-cy[0]*cx[2]*b[1])/den;

	rc += ra*boxleft + rb*boxup;
	gc += ga*boxleft + gb*boxup;
	bc += ba*boxleft + bb*boxup;
	ra = ra/im.getwidth()*(boxright-boxleft);
	ga = ga/im.getwidth()*(boxright-boxleft);
	ba = ba/im.getwidth()*(boxright-boxleft);
	rb = rb/im.getheight()*(boxdown-boxup);
	gb = gb/im.getheight()*(boxdown-boxup);
	bb = bb/im.getheight()*(boxdown-boxup);

	if (y[1]<y[0]) {
		if (y[2]<y[1]) {
			topx = (float)im.getwidth()*(x[2]-boxleft)/
				(boxright-boxleft);
			topy = (float)im.getheight()*(y[2]-boxup)/
				(boxdown-boxup);
			if ((x[0]-x[2])/(y[0]-y[2])<(x[1]-x[2])/(y[1]-y[2])) {
				leftx = (float)im.getwidth()*(x[0]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[1]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[0]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[1]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)im.getwidth()*(x[1]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[0]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[1]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[0]-boxup)/
					(boxdown-boxup);
			}
		} else {
			topx = (float)im.getwidth()*(x[1]-boxleft)/
				(boxright-boxleft);
			topy = (float)im.getheight()*(y[1]-boxup)/
				(boxdown-boxup);
			if ((x[0]-x[1])/(y[0]-y[1])<(x[2]-x[1])/(y[2]-y[1])) {
				leftx = (float)im.getwidth()*(x[0]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[2]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[0]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[2]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)im.getwidth()*(x[2]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[0]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[2]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[0]-boxup)/
					(boxdown-boxup);
			}
		}
	} else {
		if (y[2]<y[0]) {
			topx = (float)im.getwidth()*(x[2]-boxleft)/
				(boxright-boxleft);
			topy = (float)im.getheight()*(y[2]-boxup)/
				(boxdown-boxup);
			if ((x[0]-x[2])/(y[0]-y[2])<(x[1]-x[2])/(y[1]-y[2])) {
				leftx = (float)im.getwidth()*(x[0]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[1]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[0]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[1]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)im.getwidth()*(x[1]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[0]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[1]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[0]-boxup)/
					(boxdown-boxup);
			}
		} else {
			topx = (float)im.getwidth()*(x[0]-boxleft)/
				(boxright-boxleft);
			topy = (float)im.getheight()*(y[0]-boxup)/
				(boxdown-boxup);
			if ((x[1]-x[0])/(y[1]-y[0])<(x[2]-x[0])/(y[2]-y[0])) {
				leftx = (float)im.getwidth()*(x[1]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[2]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[1]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[2]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)im.getwidth()*(x[2]-boxleft)/
					(boxright-boxleft);
				rightx = (float)im.getwidth()*(x[1]-boxleft)/
					(boxright-boxleft);
				lefty = (float)im.getheight()*(y[2]-boxup)/
					(boxdown-boxup);
				righty = (float)im.getheight()*(y[1]-boxup)/
					(boxdown-boxup);
			}
		}
	}
	plb = prb = topx;
	findline(topx,topy,rightx,righty,rm,rk);
	findline(topx,topy,leftx,lefty,lm,lk);
	ryf = (int)(righty);
	lyf = (int)(lefty);
	int lc = 0;
	for(int cy=(int)(topy);cy<=ryf||cy<=lyf;cy++) {
		if (cy==ryf && cy==lyf) {
			int plx=(int)min3(leftx,plb,plb);
			int prx=(int)max3(rightx,prb,prb);
			lc += prx-plx+1;
			drawhline(plx,prx,cy,
				plx*ra+cy*rb+rc,
				plx*ga+cy*gb+gc,
				plx*ba+cy*bb+bc,
				prx*ra+cy*rb+rc,
				prx*ga+cy*gb+gc,
				prx*ba+cy*bb+bc);
		}
		if (cy==ryf) {
			trb = rightx;
			fl++;
			findline(rightx,righty,leftx,lefty,rm,rk);
		} else trb = prb;
		if (cy==lyf) {
			tlb = leftx;
			fl++;
			findline(rightx,righty,leftx,lefty,lm,lk);
		} else tlb = plb;
		if (fl<2) {
			clb = lm*(cy+1) + lk;
			crb = rm*(cy+1) + rk;
		} else {
			clb = tlb;
			crb = trb;
		}
		tlb = min3(clb,tlb,plb);
		trb = max3(crb,trb,prb);
		if (tlb<=trb) {
			lc += (int)trb - (int)tlb + 1;
			int plx = (int)tlb;
			int prx = (int)trb;
			drawhline(plx,prx,cy,
				plx*ra+cy*rb+rc,
				plx*ga+cy*gb+gc,
				plx*ba+cy*bb+bc,
				prx*ra+cy*rb+rc,
				prx*ga+cy*gb+gc,
				prx*ba+cy*bb+bc);
		}
		plb = clb; prb = crb;
	}
	if (lc==0)
		 cout << "problem right here" << endl;
		
}
