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

#include <stdio.h>
#include <fstream>
#include <containers/ilist.h>
#include <loaders/penreader.h>

#define DEPTH 0.001

using namespace std;

penreader::penreader(const char *filename, bool hashon) {
	FILE *fp; // yes, we are doing this the old-fashioned way!
	ilist<ilist<int> *> strokes;
	int cs = -1;
	char buffer[100];
	int x,y;
	int maxx,maxy,minx,miny;
	FP scale,dx,dy;

	fp = fopen(filename,"r");
	if (fp==NULL) { ret = NULL; return; }
	maxx = maxy= -9999;
	minx = miny = 9999;
	while(fscanf(fp,"%[^\r\n]%*[\r\n]",buffer)>0) {
		if (buffer[0]=='S') {
			strokes += new ilist<int>;
			cs++;
		} else {
			sscanf(buffer,"%d%d%*d",&x,&y);
			if (x > maxx) maxx = x;
			if (x < minx) minx = x;
			if (y > maxy) maxy = y;
			if (y < miny) miny = y;
			*(strokes[cs]) += x;
			*(strokes[cs]) += y;
		}
	}
	vec p1,p2,p3;
	p1 = 0;
	p1[0] = p1[1] = -0.6;
	p2 = 1;
	p2[0] = p2[1] = 0.6;
	p3 = 0.5;
	ret = new trimesh(p1,p2,p3,hashon);
	if (maxx-minx < maxy-miny) {
		scale = 1.0/(maxy-miny);
		dx = (1.0 - (maxx-minx)*scale)/2.0 - 0.5;
		dy = -0.5;
	} else {
		scale = 1.0/(maxx-minx);
		dx = -0.5;
		dy = (1.0 - (maxy-miny)*scale)/2.0 - 0.5;
	}
	p1 = 0.9;
	p2 = 0.9;
	p3 = 0.9;
	for(cs=0;cs<strokes.length();cs++) {
		for(int i=2;i<strokes[cs]->length();i+=2) {
			p1[0] = dx+scale*((*(strokes[cs]))[i-2]-minx);
			p1[1] = -(dy+scale*((*(strokes[cs]))[i-1]-miny));
			p1[2] = DEPTH;
			// for whatever reason, this next statement is
			// not properly optimized by gcc (it removes the
			// statement and replaces it with something completely
			// wrong -- I've checked... so I'm doing something
			// else:
			//p2 = p1;
			for(int kludge=0;kludge<DIM;kludge++)
				p2[kludge] = p1[kludge];
			p2[2] = 0;
			p3[0] = dx+scale*((*(strokes[cs]))[i]-minx);
			p3[1] = -(dy+scale*((*(strokes[cs]))[i+1]-miny));
			p3[2] = 0;
			face *mf;
			mf = ret->addface(p1,p2,p3);
			p2 = p3;
			p2[2] = DEPTH;
			mf = ret->addface(p2,p1,p3);
		}
		delete strokes[cs];
	}
	fclose(fp);
}


penreader::~penreader() {
}

