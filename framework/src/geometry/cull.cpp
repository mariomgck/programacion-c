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
#include <geometry/cull.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

cull::cull(trirender *renderer,
		int res, int thresh, bool liberal, bool clip,
		double minx, double maxx, double miny, double maxy) {
	image = NULL;
	makenumtable();
	RES = ((res+7)/8)*8;
	RES8 = RES/8;
	THRESH = thresh;
	LIBERAL = liberal;
	head = NULL;
	rend = renderer;
	clipbox = clip;
	up = miny;
	down = maxy;
	left = minx;
	right = maxx;
}

cull::~cull() {
	if (image) delete []image;
	if (head) {
		tri *curr,*next;
		for(curr=head;curr;curr=next) {
			next = curr->next;
			delete curr;
		}
	}
}

void cull::start(double minx, double maxx, double miny, double maxy) {
	if (head) {
		tri *curr,*next;
		for(curr=head;curr;curr=next) {
			next = curr->next;
			delete curr;
		}
		head = NULL;
	}
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

void cull::rendertriangle(double x[3], double y[3], bool e[3],
                double cx[3], double cy[3],
                double r[3], double g[3], double b[3]) {

	tri *newone;

	newone = new tri;
	newone->prev = NULL;
	newone->next = head;
	if (head!=NULL) head->prev = newone;
	head = newone;
	for(int i=0;i<3;i++) {
		head->x[i] = x[i];
		head->y[i] = y[i];
		head->cx[i] = cx[i];
		head->cy[i] = cy[i];
		head->r[i] = r[i];
		head->g[i] = g[i];
		head->b[i] = b[i];
		head->e[i] = e[i];
	}
}

void cull::end() {
	tri *curr;

	removetriangles();
	rend->start(minxs,maxxs,minys,maxys);
	for(curr=head;curr&&curr->next;curr=curr->next);
	for(;curr;curr=curr->prev)
		rend->rendertriangle(curr->x,curr->y,curr->e,
			curr->cx,curr->cy,curr->r,curr->g,curr->b);
	rend->end();
}

void cull::makenumtable(void) {

	int i;

	for(i=0;i<256;i++) {
		numtable[i] = 0;
		if (i&1) numtable[i]++;
		if (i&2) numtable[i]++;
		if (i&4) numtable[i]++;
		if (i&8) numtable[i]++;
		if (i&16) numtable[i]++;
		if (i&32) numtable[i]++;
		if (i&64) numtable[i]++;
		if (i&128) numtable[i]++;
	}
}

#define NUMON(m) (numtable[m])

int cull::checkhline(int x1, int x2, int y) {

	int xd1,xd2,xr1,xr2;
	int x;
	unsigned char mask;
	int ypos;
	int ret=0;

	//drawhline(x1,x2,y);
	if (x1>x2) { return 0; } //x = x1; x1 = x2; x2 = x; }
	if (y<0 || y>=RES || x1>=RES || x2<0) return clipbox?0:THRESH+1;
	if (x1<0) {
		if (clipbox) x1 = 0;
		else return THRESH+1;
	}
	if (x2>=RES) {
		if (clipbox) x2 = RES;
		else return THRESH+1;
	}
	ypos = y*RES8;
	xd1 = x1>>3;  xd2 = x2>>3;
	xr1 = x1&0x07; xr2 = x2&0x07;
	if (xd1==xd2) {
		mask = ((unsigned char)0xff>>(7-xr2))>>xr1<<xr1;
		return NUMON(mask) - NUMON(image[ypos+xd1] & mask);
	}
	mask = ((unsigned char)0xff>>xr1)<<xr1;
	ret += NUMON(mask) - NUMON(image[ypos+xd1] & mask);
	for(x=xd1+1;x<xd2;x++) {
		ret += 8-NUMON(image[ypos+x]);
	}
	mask = ((unsigned char)0xff)>>(7-xr2);
	ret += NUMON(mask) - NUMON(image[ypos+xd2] & mask);
	return ret;
}

void cull::drawhline(int x1, int x2, int y) {

	int xd1,xd2,xr1,xr2;
	int x;
	unsigned char mask;
	int ypos;

	if (x1>x2) { return; } // {x = x1; x1 = x2; x2 = x; }
	if (y<0 || y>=RES || x1>=RES || x2<0) return;
	if (x1<0) x1 = 0;
	if (x2>=RES) x2 = RES-1;
	ypos = y*RES8;
	xd1 = x1>>3;  xd2 = x2>>3;
	xr1 = x1&0x07; xr2 = x2&0x07;
	if (xd1==xd2) {
		mask = ((unsigned char)0xff>>(7-xr2))>>xr1<<xr1;
		image[ypos+xd1] |= mask;
		return;
	}
	mask = ((unsigned char)0xff>>xr1)<<xr1;
	//assert(mask!=0xff || xr1==0);
	image[ypos+xd1] |= mask;
	for(x=xd1+1;x<xd2;x++) {
		image[ypos+x] = 0xff;
	}
	mask = ((unsigned char)0xff)>>(7-xr2);
	image[ypos+xd2] |= mask;
/*
	for(x=x1;x<=x2;x++)
		image[y*RES8 + (x>>3)] |= 1<<(x&0x07);
*/
}


#define min3(a,b,c) ((a<b)?((a<c)?a:c):((b<c)?b:c))
#define max3(a,b,c) ((a>b)?((a>c)?a:c):((b>c)?b:c))

void cull::findline(float x1, float y1, float x2, float y2, 
		float &m, float &b) {

/*
	if (x1<x2) {
*/
		m = (x1-x2)/(y1-y2);
		b = x1 - m*y1;
/*
	} else if (x1>x2 || y1>y2) {	
		m = (x2-x1)/(y2-y1);
		b = x2 - m*y2;
	} else {
		m = (x1-x2)/(y1-y2);
		b = x1 - m*y1;
	}
*/
}

bool cull::drawtriangle(cull::tri *curr) {

	float topx,topy,leftx,lefty,rightx,righty;
	float rm,lm,rb,lb,plb,prb,clb,crb,trb,tlb;
	int b1,b2;
	int y,lyf,ryf,fl=0;
	int ret = 0;

	if (curr->y[1]<curr->y[0]) {
		if (curr->y[2]<curr->y[1]) {
			topx = (float)RES*(curr->x[2]-boxleft)/
				(boxright-boxleft);
			topy = (float)RES*(curr->y[2]-boxup)/
				(boxdown-boxup);
			if ((curr->x[0]-curr->x[2])/(curr->y[0]-curr->y[2])<
			    (curr->x[1]-curr->x[2])/(curr->y[1]-curr->y[2])) {
				leftx = (float)RES*(curr->x[0]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[1]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[0]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[1]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)RES*(curr->x[1]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[0]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[1]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[0]-boxup)/
					(boxdown-boxup);
			}
		} else {
			topx = (float)RES*(curr->x[1]-boxleft)/
				(boxright-boxleft);
			topy = (float)RES*(curr->y[1]-boxup)/
				(boxdown-boxup);
			if ((curr->x[0]-curr->x[1])/(curr->y[0]-curr->y[1])<
			    (curr->x[2]-curr->x[1])/(curr->y[2]-curr->y[1])) {
				leftx = (float)RES*(curr->x[0]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[2]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[0]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[2]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)RES*(curr->x[2]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[0]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[2]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[0]-boxup)/
					(boxdown-boxup);
			}
		}
	} else {
		if (curr->y[2]<curr->y[0]) {
			topx = (float)RES*(curr->x[2]-boxleft)/
				(boxright-boxleft);
			topy = (float)RES*(curr->y[2]-boxup)/
				(boxdown-boxup);
			if ((curr->x[0]-curr->x[2])/(curr->y[0]-curr->y[2])<
			    (curr->x[1]-curr->x[2])/(curr->y[1]-curr->y[2])) {
				leftx = (float)RES*(curr->x[0]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[1]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[0]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[1]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)RES*(curr->x[1]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[0]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[1]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[0]-boxup)/
					(boxdown-boxup);
			}
		} else {
			topx = (float)RES*(curr->x[0]-boxleft)/
				(boxright-boxleft);
			topy = (float)RES*(curr->y[0]-boxup)/
				(boxdown-boxup);
			if ((curr->x[1]-curr->x[0])/(curr->y[1]-curr->y[0])<
			    (curr->x[2]-curr->x[0])/(curr->y[2]-curr->y[0])) {
				leftx = (float)RES*(curr->x[1]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[2]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[1]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[2]-boxup)/
					(boxdown-boxup);
			} else {
				leftx = (float)RES*(curr->x[2]-boxleft)/
					(boxright-boxleft);
				rightx = (float)RES*(curr->x[1]-boxleft)/
					(boxright-boxleft);
				lefty = (float)RES*(curr->y[2]-boxup)/
					(boxdown-boxup);
				righty = (float)RES*(curr->y[1]-boxup)/
					(boxdown-boxup);
			}
		}
	}
	plb = prb = topx;
	findline(topx,topy,rightx,righty,rm,rb);
	findline(topx,topy,leftx,lefty,lm,lb);
	ryf = (int)(righty);
	lyf = (int)(lefty);
	for(y=(int)(topy);y<=ryf||y<=lyf;y++) {
		if (y==ryf && y==lyf) {
			if (ret<=THRESH)
				ret += checkhline((int)(min3(leftx,plb,plb)),
					(int)(max3(rightx,prb,prb)),y);
			if (LIBERAL) drawhline((int)(min3(leftx,plb,plb)),
					(int)(max3(rightx,prb,prb)),y);
			else drawhline((int)(max3(leftx,plb,plb)+1),
					(int)(min3(rightx,prb,prb)-1),y);
			break;
		}
		if (y==ryf) {
			trb = rightx;
			fl++;
			findline(rightx,righty,leftx,lefty,rm,rb);
		} else trb = prb;
		if (y==lyf) {
			tlb = leftx;
			fl++;
			findline(rightx,righty,leftx,lefty,lm,lb);
		} else tlb = plb;
		if (fl<2) {
			clb = lm*(y+1) + lb;
			crb = rm*(y+1) + rb;
		} else {
			clb = tlb;
			crb = trb;
		}
		b1 = (int)(min3(clb,tlb,plb));
		b2 = (int)(max3(crb,trb,prb));
		if (ret<=THRESH)
			ret += checkhline(b1,b2,y);
		if (LIBERAL)
			drawhline(b1,b2,y);
		else {
			tlb = max3(clb,tlb,plb)+1;
			trb = min3(crb,trb,prb)-1;
			if (tlb<=trb) drawhline((int)(tlb),
				(int)(trb),y);
		}
		plb = clb; prb = crb;
	}
	return ret>THRESH;
}
	
void cull::removetriangles() {

	tri *curr,*next;
	int i;

	//for(i=0,curr=head;curr;curr=curr->next,i++);
	//cout << i << endl;
	image = new unsigned char[RES*RES8];
	if (!image) {
		cerr << "ran out of memory allocating image" << endl;
		exit(1);
	}
	for(i=0;i<RES*RES/8;image[i++] = 0);
	//i = 0;
	for(curr=head;curr;curr=next) {
		next = curr->next;
		if (!drawtriangle(curr)) {
			//i++;
			if (curr==head) {
				head = head->next;
				head->prev = NULL;
			} else {
				if (curr->prev) curr->prev->next = curr->next;
				if (curr->next) curr->next->prev = curr->prev;
			}
			delete curr;
		}
	}
	//for(i=0,curr=head;curr;curr=curr->next,i++);
	delete []image;
	image = NULL;
	return;
}
