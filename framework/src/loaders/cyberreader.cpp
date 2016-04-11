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

#include <string.h>
#include <iostream>
#include <fstream>
//#include <ios.h>

#include <containers/ilist.h>
#include <loaders/cyberreader.h>
#include <loaders/ppmreader.h>

using namespace std;

cyberreader::cyberreader(istream &cyberfile, istream &ppmfile, FP scale, FP gamma)
{
	ios_base::fmtflags fl = cyberfile.flags();
	cyberfile.unsetf(ios_base::skipws);
	
	readcyberheader(cyberfile,scale);
	readcyberdata(cyberfile,scale);
	
	cyberfile.flags(fl);

	r = new FP[numlat*numlon];
	g = new FP[numlat*numlon];
	b = new FP[numlat*numlon];
	ppmreader colorreader(ppmfile,r,g,b,gamma);

	// this is temporary to allow easy cropping of the data:
	for(int i=0;i<numlat*numlon;i++)
		if (r[i]==1.0 && g[i]==1.0 && b[i]==1.0) {
			radius[i] = hugenum;
			if (i>0) radius[i-1] = hugenum;
			if (i<numlat*numlon-1) radius[i+1] = hugenum;
			if (i>numlat-1) radius[i-numlat] = hugenum;
			
			if (i<numlat*numlon-numlat) radius[i+numlat] = hugenum;
		}

	this->gamma = gamma;
	this->scale = scale;

}

cyberreader::cyberreader(istream &cyberfile, FP scale) {

	ios_base::fmtflags fl = cyberfile.flags();
	cyberfile.unsetf(ios_base::skipws);
	
	readcyberheader(cyberfile,scale);
	readcyberdata(cyberfile,scale);
	
	cyberfile.flags(fl);
	r = new FP[numlat*numlon];
	g = new FP[numlat*numlon];
	b = new FP[numlat*numlon];
	for(int i=0;i<numlat*numlon;i++)
		r[i] = g[i] = b[i] = (FP)0.9;
	this->gamma = 1.0;
	this->scale = scale;
}

cyberreader::~cyberreader()
{
	delete[] radius;
	delete[] r;
	delete[] g;
	delete[] b;
}



// I'm doing this the slow way right now -- maybe I'll read into a buffer and
// access if I'm feeling like it later -- I just want to make sure I have it 
// right for the moment
void cyberreader::readcyberheader(istream &s, FP scale) {
	unsigned char ch;

	s.ignore(4); // skip 4 byte pointer
	s >> ch;
	if (ch=='r') {
		readnewcyberheader(s,scale);
		return;
	}

	int offset = ch;

	for(int i=1;i<4;i++) { s >> ch; offset=(offset<<8) | ch; } // read data offset
	s.ignore(40+4+2+2+1+1); // skip name, time, camera, setup code, saved char, valid char
	numlat = 0;
	for(int i=0;i<2;i++) { s >> ch; numlat=(numlat<<8) | ch; } // read number of lat. intervals
	numlon = 0;
	for(int i=0;i<2;i++) { s >> ch; numlon=(numlon<<8) | ch; } // read number of long. intervals
	rshift = 0;
	for(int i=0;i<2;i++) { s >> ch; rshift=(rshift<<8) | ch; } // read shift on radius values
	s.ignore(2+4+4+4); // ignore long. shift, flags, and both increments
	totalheight = 0;
	for(int i=0;i<4;i++) { s >> ch; totalheight<<=8; totalheight|=ch; } // read scan height
	s.ignore(offset-(4+4+40+4+2+2+1+1+2+2+2+2+4+4+4));
	thetamult = 2*M_PI/numlon;
	hmult = ((FP)totalheight/numlat)*scale;
}

void cyberreader::readnewcyberheader(istream &s, FP scale) {
	string arg;

	do {
		getline(s,arg,'\n');
		getline(s,arg,'=');
		if (arg=="NLG")
			s >> numlon;
		else if (arg=="NLT")
			s >> numlat;
		else if (arg=="RSHIFT")
			s >> rshift;
		else if (arg=="LTINCR")
			s >> totalheight; // just to hold until we can use numlat to compute
	} while (arg!="DATA");
	totalheight *= numlat; // compute total from inc*num
	getline(s,arg,'\n');
	thetamult = 2*M_PI/numlon;
	hmult = ((FP)totalheight/numlat)*scale;
}
 
void cyberreader::readcyberdata(istream &s, FP scale) {

	int lon,lat,a,b;
	int val;
	unsigned char *data = new unsigned char[numlon*numlat*2],ch;

	mr = 0;
	radius = new FP[numlon*numlat];
	a = numlon*numlat-numlat; // we must flip around the theta direction for reasons unknown
	b = 0;
	s.read((char *)data,numlon*numlat*2);
	for(lon=0;lon<numlon;lon++,a-=numlat*2)
		for(lat=0;lat<numlat;lat++,a++,b++) {
			val = data[b<<1];
			ch = data[(b<<1) + 1];
			val = (val << 8) | ch;
			if (val&0x8000) radius[a] = hugenum;
			else {
				val <<= rshift;
				radius[a] = val*scale;
				if (radius[a]>mr) mr=radius[a];
			}
		}
	delete[] data;
}

static inline int followref(const ilist<int> &r, int i) {
	while(r(i)<=0) i = -r(i);
	return i;
}

static inline int merge(ilist<int> &r, int i1, int i2) {
	int t1 = followref(r,i1);
	int t2 = followref(r,i2);
	if (t1!=t2) {
		r[t1] += r[t2];
		r[t2] = -t1;
	}
	return t1;
}

void cyberreader::reduce() {
	int *label;
	ilist<int> ref;

	int i=0,t,z;
	label = new int[numlat*numlon];
	for(t=0;t<numlon;t++) {
		for(z=0;z<numlat;z++,i++) {
			if (radius[i]==hugenum) label[i] = -1;
			else {
				if (t==0) {
					label[i] = ref.length();
					ref += 1;
				} else {
					if (z==0) {
						if (label[i-numlat]==-1 || label[i-numlat+1]==-1) {
							label[i] = ref.length();
							ref += 1;
						} else {
							label[i] = merge(ref,label[i-numlat],label[i-numlat+1]);
							ref[label[i]]++;
						}
					} else if (z==numlat-1) {
						if (label[i-1]==-1 || label[i-numlat]==-1) {
							label[i] = ref.length();
							ref += 1;
						} else {
							label[i] = merge(ref,label[i-numlat],label[i-1]);
							ref[label[i]]++;
						}
					} else {
						if ((label[i-1]==-1 && label[i-numlat+1]==-1) || label[i-numlat]==-1) {
							label[i] = ref.length();
							ref += 1;
						} else if (label[i-1]==-1) {
							label[i] = merge(ref,label[i-numlat],label[i-numlat+1]);
							ref[label[i]]++;
						} else if (label[i-numlat+1]==-1) {
							label[i] = merge(ref,label[i-numlat],label[i-1]);
							ref[label[i]]++;
						} else {
							label[i] = merge(ref,merge(ref,label[i-numlat],label[i-1]),
							                     label[i-numlat+1]);
							ref[label[i]]++;
						}
					}
				}
			}
		}
	}
	i = (numlon-2)*numlat;
	for(z=0;z<numlat;z++) {
		if (label[z]==-1 || label[i+z]==-1) continue;
		if (z==0) {
			if (label[1]!=-1) {
				merge(ref,merge(ref,label[0],label[i+1]),label[i]);
			}
		} else if (z==numlat-1) {
			if (label[z-1]!=-1) {
				merge(ref,merge(ref,label[z],label[z-1]),label[i+z]);
			}
		} else {
			if (label[z-1]!=-1 && label[i+z+1]!=-1) {
				merge(ref,merge(ref,label[z],label[z-1]),merge(ref,label[i+z],label[i+z+1]));
			} else if (label[z-1]!=-1) {
				merge(ref,merge(ref,label[z],label[i+z]),label[z-1]);
			} else if (label[i+z+1]!=-1) {
				merge(ref,merge(ref,label[z],label[i+z]),label[i+z+1]);
			}
		}
	}
	i=0;
	int j,p;
	p = 0;
	for(j=1;j<ref.length();j++)
		if (ref[j]>ref[p]) p=j;
	for(t=0;t<numlon;t++)
		for(z=0;z<numlat;z++,i++) 
			if (label[i]==-1 || followref(ref,label[i])!=p) radius[i]=hugenum;
	delete []label;

}

void cyberreader::fillholes() {

	int *h = findholes();
	FP change,old;
	int z,t,i;
	int ittnum=0;

	//ppmreader outf(numlon,numlat);

	i = 0;
	for(t=0;t<numlon;t++)
		for(z=0;z<numlat;z++,i++)
			if (h[i]==1) radius[i] = 0;
	//		else outf.red(t,z) = outf.green(t,z) = outf.blue(t,z) = 1.0;

	//ofstream of("check.ppm");
	//outf.save(of,true,false);
	//of.close();
	//exit(1);

	do {
		//cout << "itt number: " << ittnum << endl;
		change = 0.0;
		i=0;
		for(t=0;t<numlon;t++) {
			i++;
			for(z=1;z<numlat-1;z++,i++) {
				if (h[i]==0) continue;
				old = radius[i];
				if (t==0) {
					radius[i] = (radius[i-1]+radius[i+1]+radius[i+numlat]+radius[i+numlat*numlon-numlat])/(FP)4.0;
				} else if (t==numlon-1) {
					radius[i] = (radius[i-1]+radius[i+1]+radius[i-numlat*numlon+numlat]+radius[i-numlat])/(FP)4.0;
				} else {
					radius[i] = (radius[i-1]+radius[i+1]+radius[i+numlat]+radius[i-numlat])/(FP)4.0;
				}
				old -= radius[i];
				if (old<0.0) old = -old;
				if (change<old) change=old;
			}
			i++;
		}
		ittnum++;
	} while (change>0.0001);
	delete []h;
}

int *cyberreader::findholes() {
	
	ilist<int> outside;
	int *labels = new int[numlat*numlon];
	int i=0;

	outside += 1;
	for(int t=0;t<numlon;t++,i+=numlat)
		if (radius[i]==hugenum) labels[i] = 0;

	for(int z=1;z<numlat;z++) {
		i=z;
		for(int t=0;t<numlon;t++,i+=numlat) {
			if (radius[i]==hugenum) {
				if (t==0) {
					if (radius[i-1]!=hugenum) {
						labels[i] = outside.length();
						outside+=1;
					} else {
						labels[i] = labels[i-1];
					}
				} else if (t==numlon-1 && radius[z]==hugenum) {
					labels[i] = labels[z];
					if (radius[i-1]==hugenum) labels[i] = merge(outside,labels[i],labels[i-1]);
					if (radius[i-numlat]==hugenum) 
						labels[i] = merge(outside,labels[i],labels[i-numlat]);
				} else {
					if (radius[i-numlat]!=hugenum && radius[i-1]!=hugenum) {
						labels[i] = outside.length();
						outside += 1;
					} else if (radius[i-numlat]!=hugenum) {
						labels[i] = labels[i-1];
					} else if (radius[i-1]!=hugenum) {
						labels[i] = labels[i-numlat];
					} else {
						labels[i] = merge(outside,labels[i-numlat],labels[i-1]);
					}
				}
				if (z==numlat-1) merge(outside,0,labels[i]);
			}
		}
	}

	i=0;
	int out = followref(outside,0);

	for(int t=0;t<numlon;t++)
		for(int z=0;z<numlat;z++,i++)
			if (radius[i]==hugenum && followref(outside,labels[i])!=out) labels[i] = 1;
			else labels[i] = 0;

	return labels;
}

void cyberreader::clip(int minheight, int maxheight) {

	for(int x=0;x<minheight;x++)
		for(int y=0;y<numlat;y++)
			radius[x+y*numlon] = hugenum;
	for (int x=maxheight+1;x<numlon;x++)
		for(int y=0;y<numlat;y++)
			radius[x+y*numlon] = hugenum;
}
