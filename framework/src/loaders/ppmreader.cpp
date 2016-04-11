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
#include <loaders/ppmreader.h>
#include <globaldef.h>
#include <iostream>
#include <float.h>
#include <math.h>

using namespace std;

ppmreader::ppmreader(istream &s,FP *r, FP *g, FP *b, FP mult) {
	owndata = false;
	this->r = r;
	this->g = g;
	this->b = b;
	loadit(s,mult);
}

ppmreader::ppmreader(istream &s, FP mult) {
	r = NULL;
	g = NULL;
	b = NULL;
	owndata = true;
	loadit(s,mult);
}

ppmreader::ppmreader(int width, int height, double initr, double initg, double initb) {
	w = width;
	h = height;
	int m=w*h;
	r = new FP[m];
	g = new FP[m];
	b = new FP[m];
	for(int i=0;i<m;i++) {
		r[i]=initr;
		g[i]=initg;
		b[i]=initb;
	}
	owndata=true;
}

ppmreader::ppmreader(int width, int height, FP *r, FP *g, FP *b) {
	w = width;
	h = height;
	owndata = false;
	this->r = r;
	this->g = g;
	this->b = b;
}

ppmreader::ppmreader(const image &im) {
	w = im.getwidth();
	h = im.getheight();
	owndata = true;
	int m = w*h;
	r = new FP[m];
	g = new FP[m];
	b = new FP[m];
	if (im.numchannels()==1) {
		for(int i=0;i<m;i++) {
			r[i] = g[i] = b[i] = im.v(i);
		}
	} else if (im.numchannels()==2) { // hmmm... what to do?
		for(int i=0;i<m;i++) {
			r[i] = g[i] = b[i] = (im.v(i,0) + im.v(i,1))/2;
		}
	} else {
		for(int i=0;i<m;i++) {
			r[i] = im.v(i,0);
			g[i] = im.v(i,1);
			b[i] = im.v(i,2);
		}
	}
}

ppmreader::~ppmreader()
{
	if (owndata) {
		delete[] r; delete[] g; delete[] b;
	}
}

static int readppmvalue(istream &s) {
	char ch;
	int ret = 0;

	do {
		do { s>>ch; } while(ch==' ' || ch=='\n' || ch=='\t' || ch=='\r');
		if (ch=='#') {
			do {
				s>>ch;
			} while(ch!='\n'&&ch!='\r');
		} else break;
	} while(1);
	while(ch>='0'&&ch<='9') { ret=ret*10+(ch-'0'); s>>ch; }
	s.putback(ch);
	return ret;
}

void ppmreader::loadit(istream &s, FP mult) {

	char c1,c2;
	
	cout << "starting loading ppm data..." << endl;
	s >> c1 >> c2;
	if (c1!='P') return;
	c2 -= '0';
	if (c2<=0 || c2>6) return;
	bool ascii = c2<4;
	bool color = c2%3==0;
	cout << "ascii: " << ascii << " and color: " << color << endl;
	std::ios_base::fmtflags fl = s.flags();

	s.unsetf(ios::skipws);
	w = readppmvalue(s);
	h = readppmvalue(s);

	cout << "width,height: " << w << "," << h << endl;
	if (w==0 || h==0) {
		s.flags(fl);
		return;
	}
	int m = w*h;
	if (owndata) {
		r = new FP[m];
		g = new FP[m];
		b = new FP[m];
	}
	int maxv;
	if (c2!=1 && c2!=4) maxv = readppmvalue(s);
	else maxv = 1;
	if (!ascii) s >> c1;

	int cr,cg,cb;
	unsigned char ccr,ccg,ccb;

	for(int i=0;i<m;i++) {
		if (ascii) {
			if (color) {
				cr = readppmvalue(s);
				cg = readppmvalue(s);
				cb = readppmvalue(s);
			} else cr=cg=cb=readppmvalue(s);
		} else {
			if (color) {
				s >> ccr >> ccg >> ccb;
				cr = ccr; cg = ccg; cb = ccb;
			} else {
				if (c2==4) {
					if (i%8==0) s >> c1;
					else c1<<=1;
					ccr = c1>>7;
				} else s >> ccr;
				cr = cg = cb = ccr;
			}
		}
		r[i] = mult*(FP)cr/(FP)maxv;
		g[i] = mult*(FP)cg/(FP)maxv;
		b[i] = mult*(FP)cb/(FP)maxv;
	}
	cout << "done reading ppm" << endl;
	s.flags(fl);
}

// the weights on the gray scale values shoudl probably not be all 1/3 for each
//   color channel.  When I have the time I'll look up the proper mixing ratio
//   (something like 0.3,0.6,0.1 I think).
void ppmreader::save(ostream &s, int range, bool ascii, bool color, FP min, FP max) {

	if (!ascii && range>255) range=255;
	if (ascii) {
		if (color) s << "P3";
		else s << "P2";
	} else {
		if (color) s << "P6";
		else s << "P5";
	}
	s << endl;
	s << w << ' ' << h << endl << range << endl;
	int m=w*h;
	int c = 0;
	FP rr,gg,bb;
	for(int i=0;i<m;i++) {
		rr=r[i]; gg=g[i]; bb=b[i];
		if (!_finite(rr) || !_finite(gg) || !_finite(bb))
			rr=gg=bb=0;
		else {
			if (rr<min) rr=min;
			if (rr>max) rr=max;
			if (gg<min) gg=min;
			if (gg>max) gg=max;
			if (bb<min) bb=min;
			if (bb>max) bb=max;
		}
		if (ascii) {
			if (color) {
				s << (int)((rr-min)/(max-min)*range) << ' ' << 
				     (int)((gg-min)/(max-min)*range) << ' ' <<
					 (int)((bb-min)/(max-min)*range);
				c += 3;
			} else {
				s << (int)((rr+gg+bb-3*min)/(max-min)*range/3);
				c++;
			}
			if (c<11) s << ' ';
			else { s << endl; c=0; }
		} else {
			if (color) {
				s << (unsigned char)((rr-min)/(max-min)*range) <<
				     (unsigned char)((gg-min)/(max-min)*range) <<
					 (unsigned char)((bb-min)/(max-min)*range);
			} else {
				s << (unsigned char)((rr+gg+bb-3*min)/(max-min)*range/3);
			}
		}
	}
}
