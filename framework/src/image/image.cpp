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
#include <image/image.h>
#include <math.h>
#include <float.h>
#include <loaders/ppmreader.h>
#include <fstream>

using namespace std;

image image::copy() const {
	FP *ch;
	image ret(width,height,xwrap,ywrap);
	int size = width*height;

	for(int c = 0;c<chs.length();c++) {
		ch = new FP[size];
		for(int i=0;i<size;i++)
			ch[i] = channel(c)[i];
		ret.chs += new chptr(ch,true);
	}
	return ret;
}

image image::reduce(int ntaps, FP *taps) const {
	FP *currblock,*tb;
	image ret((width+1)/2,(height+1)/2,xwrap,ywrap);

	for(int c = 0;c<chs.length();c++) {
		tb = filterrows(channel(c),width,height,ntaps,taps);
		currblock = halfrows(tb,width,height);
		delete []tb;
		tb = filtercolumns(currblock,width,(height+1)/2,ntaps,taps);
		delete []currblock;
		currblock = halfcolumns(tb,width,(height+1)/2);
		delete []tb;
		ret.chs += new chptr(currblock,true);
	}
	return ret;
}
	
image image::expand(int ntaps, FP *taps, int w, int h) const {

	FP *currblock,*tb;
	
	if (w==-1) w = width*2;
	if (h==-1) h = height*2;

	image ret(w,h,xwrap,ywrap);

	for(int c=0;c<chs.length();c++) {
		tb = doublerows(channel(c),width,height,h);
		currblock = filterrows(tb,width,h,ntaps,taps);
		delete []tb;
		//currblock = tb;
		tb = doublecolumns(currblock,width,h,w);
		delete []currblock;
		currblock = filtercolumns(tb,w,h,ntaps,taps);
		delete []tb;
		//currblock = tb;
		ret.chs += new chptr(currblock,true);
	}
	ret.mult(4);
	return ret;
}

FP *image::filterrows(FP *ch, int w, int h, int ntaps, FP *taps) const {
	FP *ret = new FP[w*h];
	int x,y,i,ti;

	int b=ntaps/2;
	int index = 0;
	if (true || ntaps>h) {
		for(y=0;y<h;y++) {
			for(x=0;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++) {
					if (ywrap) ret[index] += taps[ti]*ch[x+((i+y+h)%h)*w];
					else if (y<-i) ret[index] += taps[ti]*ch[x];
					else if (i+y>=h) ret[index] += taps[ti]*ch[w*h-w+x];
					else ret[index] += taps[ti]*ch[index+w*i];
				}
			}
		}
		return ret;
	}
	if (ywrap) {
		for(y=0;y<b;y++) {
			for(x=0;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++) 
					ret[index] += taps[ti]*ch[x+((i+y+h)%h)*w];
			}
		}
	} else {
		for(y=0;y<b;y++) {
			for(x=0;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++)
					ret[index] += y<-i ? taps[ti]*ch[x] : taps[ti]*ch[index+i*w];
			}
		}
	}
	for(;y<h-b;y++) {
		for(x=0;x<w;x++,index++) {
			ret[index]=0;
			for(ti=0,i=-b;i<=b;i++,ti++)
				ret[index] += taps[ti]*ch[index+i*w];
		}
	}
	if (ywrap) {
		for(;y<h;y++) {
			for(x=0;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++)
					ret[index] += taps[ti]*ch[x+((i+y+h)%h)*w];
			}
		}
	} else {
		int last = w*h-w;
		for(;y<h;y++) {
			for(x=0;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++)
					ret[index] += y+i>=h ? taps[ti]*ch[last+x] : taps[ti]*ch[index+i*w];
			}
		}
	}
	return ret;
}
FP *image::filtercolumns(FP *ch, int w, int h, int ntaps, FP *taps) const {
	FP *ret = new FP[w*h];
	int x,y,i,ti;

	int b=ntaps/2;
	int index = 0;
	if (true || ntaps>w) {
		for(y=0;y<h;y++) {
			for(x=0;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++) {
					if (xwrap) ret[index] += taps[ti]*ch[index-x+(i+x+w)%w];
					else if (x<-i) ret[index]+=0; //ret[index] += taps[ti]*ch[y*w];
					else if (i+x>=w) ret[index]+=0; //ret[index] += taps[ti]*ch[y*w+w-1];
					else ret[index] += taps[ti]*ch[index+i];
				}
			}
		}
		return ret;
	}
	for(y=0;y<h;y++) {
		if (xwrap) {
			for(x=0;x<b;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++) 
					ret[index] += taps[ti]*ch[index-x+(i+x+w)%w];
			}
		} else {
			for(x=0;x<b;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++)
					ret[index] += x<-i ? taps[ti]*ch[y*w] : taps[ti]*ch[index+i];
			}
		}
		for(;x<w-b;x++,index++) {
			ret[index]=0;
			for(ti=0,i=-b;i<=b;i++,ti++)
				ret[index] += taps[ti]*ch[index+i];
		}
		if (xwrap) {
			for(;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++)
					ret[index] += taps[ti]*ch[index-x+(i+x+w)%w];
			}
		} else {
			for(;x<w;x++,index++) {
				ret[index]=0;
				for(ti=0,i=-b;i<=b;i++,ti++)
					ret[index] += x+i>=w ? taps[ti]*ch[y*w+w-1] : taps[ti]*ch[index+i];
			}
		}
	}
	return ret;
}

FP *image::halfrows(FP *ch, int w, int h) const {
	FP *ret = new FP[w*(h+1)/2];
	int ri,ci;
	int x,y;

	ri=ci=0;
	for(y=0;y<(h+1)/2;y++) {
		for(x=0;x<w;x++,ri++,ci++)
			ret[ri] = ch[ci];
		ci+=w;
	}
	return ret;
}

FP *image::halfcolumns(FP *ch, int w, int h) const {
	FP *ret = new FP[(w+1)/2*h];
	int ri,ci;
	int x,y;
	const int fix = w%2 ? -1 : 0;

	ri=ci=0;
	for(y=0;y<h;y++) {
		for(x=0;x<(w+1)/2;x++,ri++,ci+=2)
			ret[ri] = ch[ci];
		ci += fix;
	}
	return ret;
}

FP *image::doublerows(FP *ch, int w, int h, int newh) const {
	FP *ret = new FP[w*newh];
	int ri,ci;
	int x,y,ry;

	ri=ci=0;
	for(y=ry=0;y<h && ry<newh;y++,ry++) {
		for(x=0;x<w;x++)
			ret[ri++] = ch[ci++];
		ry++;
		if (ry<newh)
			for(x=0;x<w;x++)
				ret[ri++] = 0;
	}
	for(;ry<newh;ry++)
		for(x=0;x<w;x++)
			ret[ri++] = 0;
	return ret;
}

FP *image::doublecolumns(FP *ch, int w, int h, int neww) const {
	FP *ret = new FP[neww*h];
	int ri,ci;
	int x,y,rx;

	ri=ci=0;
	for(y=0;y<h;y++) {
		for(x=rx=0;x<w && rx<neww;x++,rx++) {
			ret[ri++] = ch[ci++];
			rx++;
			if (rx<neww) ret[ri++] = 0;
		}
		for(;rx<neww;rx++)
			ret[ri++] = 0;
	}
	return ret;
}

image image::backwarp(const image &backflow, FP scale) const {
	if (backflow.chs.length()!=2 || backflow.width!=width || backflow.height!=height) 
		return copy();

	image ret(width,height,xwrap,ywrap);
	int x,y,x1,y1,x2,y2,c;
	int pos;
	FP xfrac,yfrac;

	for(c=0;c<chs.length();c++) ret.addchannel();

	pos = 0;
	for(y=0;y<height;y++)
		for(x=0;x<width;x++,pos++) {
			if (!_finite(backflow.channel(0)[pos]) ||
				!_finite(backflow.channel(1)[pos])) xfrac=yfrac=0;
			else {
				xfrac = (FP)x + backflow.channel(0)[pos]*scale;
				yfrac = (FP)y + backflow.channel(1)[pos]*scale;
			}
			x1 = (int)floor(xfrac);
			y1 = (int)floor(yfrac);
			xfrac -= x1;
			yfrac -= y1;
			x2 = x1+1;
			y2 = y1+1;
			if (x1<0) {
				if (xwrap) x1 += width;
				else x1 = 0;
			} else if (x1>=width) {
				if (xwrap) x1 -= width;
				else x1 = width-1;
			}
			if (x2<0) {
				if (xwrap) x2 += width;
				else x2 = 0;
			} else if (x2>=width) {
				if (xwrap) x2 -= width;
				else x2 = width-1;
			}
			if (y1<0) {
				if (ywrap) y1 += height;
				else y1 = 0;
			} else if (y1>=height) {
				if (ywrap) y1 -= height;
				else y1 = height-1;
			}
			if (y2<0) {
				if (ywrap) y2 += height;
				else y2 = 0;
			} else if (y2>=height) {
				if (ywrap) y2 -= height;
				else y2 = height-1;
			}
			for(c=0;c<chs.length();c++)
				ret.channel(c)[pos] =
					((FP)1.0-yfrac)*(value(x1,y1,c)*((FP)1.0-xfrac) + value(x2,y1,c)*xfrac) +
					yfrac*(value(x1,y2,c)*((FP)1.0-xfrac) + value(x2,y2,c)*xfrac);
		}

	return ret;
}

image image::dx() const {

	image ret(width,height,xwrap,ywrap);

	int c,x,y,pos;

	for(c=0;c<chs.length();c++) ret.addchannel();

	pos=0;
	for(y=0;y<height;y++)
		for(x=0;x<width;x++,pos++) {
			if (x==0) {
				if (xwrap) {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = (channel(c)[pos+1] - channel(c)[pos+width-1])/2;
				} else {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = 0;
				}
			} else if (x==width-1) {
				if (xwrap) {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = (channel(c)[pos-width+1] - channel(c)[pos-1])/2;
				} else {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = 0;
				}
			} else {
				for(c=0;c<chs.length();c++)
					ret.channel(c)[pos] = (channel(c)[pos+1] - channel(c)[pos-1])/2;
			}
		}
	return ret;
}

image image::dy() const {

	image ret(width,height,xwrap,ywrap);

	int c,x,y,pos,size=width*height;

	for(c=0;c<chs.length();c++) ret.addchannel();

	pos=0;
	for(y=0;y<height;y++)
		for(x=0;x<width;x++,pos++) {
			if (y==0) {
				if (ywrap) {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = (channel(c)[pos+width] - channel(c)[size-width+pos])/2;
				} else {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = 0;
				}
			} else if (y==height-1) {
				if (ywrap) {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = (channel(c)[y] - channel(c)[pos-width])/2;
				} else {
					for(c=0;c<chs.length();c++)
						ret.channel(c)[pos] = 0;
				}
			} else {
				for(c=0;c<chs.length();c++)
					ret.channel(c)[pos] = (channel(c)[pos+width] - channel(c)[pos-width])/2;
			}
		}
	return ret;
}

// I'm assuming for the moment that size>width and size>height
image image::squaresum(int size, bool onlyfinite) const {
	image ret(width,height,xwrap,ywrap);

	int c,x,y,pos,tx,ty,ttx;
	FP val,tv;
	int i,d=size/2;

	for(c=0;c<chs.length();c++) ret.addchannel();

	for(c=0;c<chs.length();c++) {
		pos = 0;
		for(y=0;y<height;y++) {
			val = 0;
			for(i=-d+y;i<=d+y;i++)
				for(x=-d;x<=d;x++) {
					if (x<0) {
						if (xwrap) tx = x+width;
						else continue;
					} else tx = x;
					if (y<0) {
						if (ywrap) ty = y+width;
						else continue;
					} else ty = y;
					if (!onlyfinite || _finite(tv=value(tx,ty,c)))
						val += tv;
				}
			ret.channel(c)[pos++] = val;
			for(x=1;x<width;x++,pos++) {
				tx = x-d-1;
				if (tx<0) {
					if (xwrap) tx += width;
					else tx = -1;
				}
				ttx = x+d;
				if (ttx>=width) {
					if (xwrap) ttx -= width;
					else ttx = -1;
				}
				for(i=-d+y;i<=d+y;i++) {
					if (i<0) {
						if (ywrap) ty = i+height;
						else continue;
					} else if (i>=height) {
						if (ywrap) ty = i-height;
						else continue;
					} else ty = i;
					if (tx!=-1 && (!onlyfinite || _finite(tv=value(tx,ty,c))))
						val -= tv;
					if (ttx!=-1 && (!onlyfinite || _finite(tv=value(ttx,ty,c))))
						val += tv;
				}
				ret.channel(c)[pos] = val;
			}
		}
	}
	return ret;
}

void image::save(const char * const name, int c1, int c2, int c3, FP min, FP max) {
	if (min==-hugenum) min = this->min();
	if (max==hugenum) max = this->max();
	ofstream of(name);
	//of.unsetf(ios::skipws);
	ppmreader p(width,height,channel(c1),channel(c2),channel(c3));
	p.save(of,255,false,c1!=c2||c2!=c3,min,max);
	of.close();
}

bool image::load(const char * const name) {
	ifstream f(name);
	f.unsetf(ios::skipws);
	ppmreader p(f);
	f.close();
	if ((width!=-1 || height!=-1) &&
		(width!=p.width() || height!=p.height())) return false;
	p.release();
	width = p.width();
	height = p.height();
	xwrap = false;
	ywrap = false;
	chs += new chptr(p.redchannel(),true);
	chs += new chptr(p.greenchannel(),true);
	chs += new chptr(p.bluechannel(),true);
	return true;
}
