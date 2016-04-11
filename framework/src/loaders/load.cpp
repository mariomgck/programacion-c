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
#include <geometry/trimesh.h>
#include <warping/warpmesh.h>
#include <loaders/ppmreader.h>
#include <loaders/penreader.h>
#include <loaders/vtkreader.h>

#include <fstream>
#include <loaders/load.h>
#include <string.h>

trimesh *loadtrimesh(const char *filename, bool &spec, bool hashon) {

	trimesh *currmesh;

	currmesh = NULL;

	char e[10];
	int i;
	for(i=strlen(filename)-1;i>0&&filename[i]!='.';i--)
		;
	if (i!=0) {
		int j=0;
		for(i++;i<=strlen(filename)&&j<9;i++,j++) e[j]=filename[i];
	} else e[0] = 0;

	if (strcmp(e,"tm")==0) {
		ifstream in(filename);
		if (in.good()) {
			currmesh = new trimesh(in,1,hashon,matrix::eye,
				vec::zero,trimesh::tmff1);
			in.close();
		}
	} else if (strcmp(e,"ppm")==0 || strcmp(e,"pgm")==0) {
		ifstream in(filename);
		in.unsetf(ios::skipws);
		ppmreader read(in);
		currmesh = new trimesh(read,1,1,hashon);
		spec = false;
	} else if (strcmp(e,"dat")==0) {
		penreader read(filename,hashon);
		currmesh = read.mesh();
		spec = false;
	} else if (strcmp(e,"vtk")==0) {
		vtkreader read(filename,hashon);
		currmesh = read.mesh();
		spec = false;
	} else { // assume it is an old mesh format
		ifstream in(filename);
		if (in.good()) {
			currmesh = new trimesh(in,1,hashon,matrix::eye,
				vec::zero,trimesh::oldff);
			if (currmesh->numvertices()==0) {
				delete currmesh;
				currmesh = NULL;
			}
			in.close();
		}
	}
	return currmesh;
}
