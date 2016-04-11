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
#include <morphing/bootstrap.h>
#include <loaders/load.h>
#include <geometry/model.h>
#include <fstream>

using namespace std;
const double hugenum = 100000;
// command line is:
//    bootmodel parameterfile outputmodel [#proc [%left [inputmodel]]]

int main(int argc, char **argv) {

	if (argc==1) {
		cout << "command line is:" << endl;
		cout << argv[0] << " parameterfile outputmodel [#proc [%left [inputmodel]]]" << endl;
		exit(1);
	}
	ifstream inf(argv[1]);
	if (!inf.good()) exit(1);
	char wpfile[100];
	inf >> wpfile;
	int n;
	inf >> n;
	ilist<trimesh *> exs;
	char buffer[100];
	bool devnull;

	/*Carga una lista de mallas*/
	/*
	Run bootstrap with a data file that contains the warp parameters, the
	number of models, and a list of each model (with the base model first).
	data/carmodel.params is the file I used.  In this case, I chose to execute
	bootmodel as:
	./bootmodel data/carmodel.params data/carmodel.mm 1 0.0
	which uses the data/carmodel.params input file and produces the
	data/carmodel.mm morphable model.  It uses only one processor (instead of
	using pthreads to split up the computation) and on the first pass it takes
	all of the eigenvectors.
	*/
	for(int i=0;i<n;i++) {
		inf >> buffer;
		cout << "reading model " << buffer << endl;
		exs += loadtrimesh(buffer,devnull,true);
	}


	inf.close();
	inf.open(wpfile);
	if (!inf.good()) exit(1);
	warperparams wp(inf);
	FP cutlevel;
	inf >> cutlevel;
	inf.close();
	cout << "building bootstrapper" << endl;
	bootstrap *builder;

	if (argc>5) {
		inf.open(argv[5]);
		if (!inf.good()) {
			cerr << "could not open model " << argv[5] << endl;
			exit(1);
		}
		cout << "initial model loading from " << argv[5] << endl;
		builder = new bootstrap(wp,cutlevel,exs[0],inf);
		cout << builder->nummodeldim() << " initial dimensions" << endl;
		inf.close();
	} else {
		cout << "using blank initial model" << endl;
		builder = new bootstrap(wp,cutlevel,exs[0]);
	}

	for(int i=1;i<n;i++) {
		cout << "adding model " << i << endl;
		builder->addobject(exs[i]);
	}
	cout << "building model" << endl;

	while(builder->iterate(argc>4?atof(argv[4]):0.95,
	                       argc>3?atoi(argv[3]):1)) {
		cout << "saving after iteration" << endl;
		ofstream outf(argv[2]);
		if (!outf.good()) exit(1);
		builder->getmodel()->save(outf);
		outf.close();
	}

	cout << "done . . . saving" << endl;

	ofstream outf(argv[2]);
	if (!outf.good()) exit(1);
	builder->getmodel()->save(outf);
	outf.close();

	return 0;
}
		
		
