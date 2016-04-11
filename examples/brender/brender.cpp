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
#include <X11/X.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <loaders/mload.h>
#include <loaders/load.h>
#include <matrix/vec.h>

#include <geometry/trimesh.h>
#include <geometry/cull.h>
#include <image/imagerender.h>

#include <iostream>
#include <fstream>
#include <string.h>

double zoom = 1;
double dirlight = 0.6;
double amblight = 0.1;
ilist<double> alpha;

model *mm;

using namespace std;

class affinevecconvert : public vecconvert {
public:
	affinevecconvert(const matrix &m=matrix::eye, const vec &v=vec::zero) {
		a = m; c = v;
	}
	virtual void convertpt(vec in, double &x, double &y, double &z,
		double &r, double &g, double &b) {
		vec out = a*in + c;
		x = out[0];
		y = out[1];
		z = out[2];
		r = out[3];
		g = out[4];
		b = out[5];
	}
	void setrot(double xrot,double yrot, double zrot, double scale) {
	
		a = matrix::eye;
		a[0][0] = scale;
		a[1][1] = scale;
		a[2][2] = scale;
		matrix rot(matrix::eye);
		rot[2][2] = rot[0][0] = cos(M_PI*yrot/180);
		rot[2][0] = -(rot[0][2] = sin(M_PI*yrot/180));
		a = rot*a;
		rot = matrix::eye;
		rot[2][2] = rot[1][1] = cos(M_PI*xrot/180);
		rot[1][2] = -(rot[2][1] = sin(M_PI*xrot/180));
		a = rot*a;
		rot = matrix::eye;
		rot[0][0] = rot[1][1] = cos(M_PI*zrot/180);
		rot[1][0] = -(rot[0][1] = sin(M_PI*zrot/180));
		a = rot*a;
	}
	matrix a;
	vec c;
};

affinevecconvert currconv;

void snapit(char *fn, int w, int h, int red, int green, int blue, bool outline, int olr, int olg, int olb, bool split) {

	ofstream outf(fn,ios::out|ios::binary);

	if (!outf.good()) {
		printf ("could not create file %s\n",fn);
		return;
	}

	imagerender imr(w,h,red,green,blue,outline,olr,olg,olb,
			-1,1,-1,1);
	mm->getmesh()->render(imr,split,&currconv,amblight,1,dirlight);
	ppmreader re(imr.getimage());
	re.save(outf,255,false,true);
	outf.close();
	printf ("rendered to %s\n",fn);
}

void renderline(ifstream &infile) {


	double td,rot;
	char filename[100];
	int h,w,i;

	infile >> filename;
	double xrot,yrot,zrot,scale,fr,b;
	int red,green,blue,olr,olg,olb;
	bool outline,split,colorize;
	char viewfile[100];
	infile >> w >> h >> viewfile;
	ifstream vf(viewfile);
	red = -10;
	if (vf.good()) {
		char buffer[100];
		vf >> buffer;
		if (strcmp(buffer,"VIEWFF1")) {
			printf ("%s is an invalid view file\n",viewfile);
			vf.close();
		} else {
			vf >> outline >> split >> fr >> b >> red >> green >> blue
				>> olr >> olg >> olb >> scale >> yrot >> xrot >> zrot
				>> colorize;
			scale = pow(2,scale);
			vf.close();
		}
	}
	if (red==-10) {
		printf ("could not open %s as a view file\n",viewfile);
		outline = false; split = false; red = green = blue = 255;
		olr = olg = olb = 0; scale = 1.0; yrot = zrot = xrot = 0.0;
		colorize = true;
	}
	currconv.setrot(xrot,yrot,zrot,scale);
	for(i=0;i<alpha.length();i++) {
		infile >> td;
		alpha.setnth(td,i);
	}
	mm->setparameters(alpha);
	printf ("rendering image %s\n",filename);
	printf ("at %dx%d from angle (%lf,%lf,%lf) at zoom level %lf and with "
		"colors %s\n",w,h,xrot,yrot,zrot,zoom,(colorize ? "on" : "off"));
		
	snapit(filename,w,h,red,green,blue,outline,olr,olg,olb,split);
}

const double hugenum = 100000;

int main(int argc, char **argv)
{
  int i,n;
  char buffer[100];

  if (argc==1) {
	printf ("usage: %s datafile\n",argv[0]);
	printf ("where datafile has:\n"
		"a single line with the filename of the morphable model,\n"
		"followed by the direct and then the ambient light levels (try 0.75 0.25),\n"
		"and then a another line with the number of images to be created.\n"
		"The file then has one line for each image with the line structure:\n"
		"filename width height viewfile alpha0 . . . alphan\n");
	exit(1);
  }
  ifstream infile(argv[1]);
  if (!infile.good()) {
	printf ("error: could not open file %s\n",argv[1]);
	exit(1);
  }
  infile >> buffer;
  infile >> dirlight >> amblight;
  infile >> n;

  warperparams wp;
  mm = loadmodel(buffer,wp);
  if (mm==NULL) {
	cerr << "could not open file " << buffer << endl;
	exit(1);
  }
  printf ("done loading model\n");

  for(i=0;i<mm->numdim();i++) {
	alpha += 0.0;
  }

  for(i=0;i<n;i++) {
	printf ("rendering image %d\n",i);
	renderline(infile);
  }
  infile.close();

  return 0;
}
