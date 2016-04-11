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
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <warping/warper.h>
#include <warping/warpparms.h>
#include <warping/warperparams.h>
#include <warping/warpmesh.h>
#include <geometry/trimesh.h>
#include <loaders/load.h>
#include <image/imagerender.h>
#include <morphing/reducer.h>
#include <morphing/quadric.h>

using namespace std;
const double hugenum = 100000;

void buildlevels(trimesh *m, ilist<trimesh *> &tl) {

	tl += new trimesh(*m,1);
	trimesh *rm = new trimesh(*m);
	quadric q(rm);
	int n = q.complexity()/2;
	bool cont;
	int i=1;
	do {
		while((cont=q.reduce()) && n<q.complexity())
			;
		tl += new trimesh(*rm,1);
		n = q.complexity()/2;
		cout << "level " << i << ": " << m->maxdistance(tl[i])/m->width() << endl;
		i++;
	} while (cont);
	delete rm;
}

int main (int argc, char **argv) {

	if (argc==1) {
		cout << "usage:" << endl;
		cout << argv[0] << " sourcemesh outstem" << endl;
		exit(0);
	}
	trimesh *from;
	cout << "loading " << argv[1] << endl;
	bool d=false;
	from = loadtrimesh(argv[1],d);
	ilist<trimesh *> frommesh;
	buildlevels(from,frommesh);
	for(int i=0;i<frommesh.length();i++) {
		char buffer[100];
		sprintf(buffer,"%s.%d.tm",argv[2],i);
		ofstream ofs(buffer);
		if (ofs.good()) {
			frommesh[i]->saveit(ofs);
			ofs.close();
		} else {
			cout << "could not open " << buffer << " for writing" << endl;
		}
	}

return 0;
}
