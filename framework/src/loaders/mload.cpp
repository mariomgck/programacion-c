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
#include <fstream>
#include <loaders/load.h>
#include <loaders/mload.h>
#include <string.h>

static void loadpts(char *filename, ilist<vec> &pts) {
	int pi = strlen(filename);
	char *pfn = new char[pi+2];
	strcpy(pfn,filename);
	pfn[pi-2] = 'p'; pfn[pi-1] = 't'; pfn[pi] = 's'; pfn[pi+1] = 0;
	ifstream is(pfn);
	if (is.good()) {
		int n;	
		vec v;
		is >> n;
		for(;n>0;n--) {
			is >> v;
			pts += v;
		}
		is.close();
	}
	delete []pfn;
}

model *loadmodel(const char *filename, char *params) {
	ifstream is (params);
	if (!is.good()) return NULL;
	warperparams p(is);
	is.close();
	return loadmodel(filename,p);
}

model *loadmodel(const char *filename, const warperparams &p) {

	model *m = NULL;
	int i;
	ifstream is;

	char e[10];
	for(i=strlen(filename)-1;i>0&&filename[i]!='.';i--)
		;
	if (i!=0) {
		int j=0;
		for(i++;i<=strlen(filename)&&j<9;i++,j++) e[j]=filename[i];
	} else e[0] = 0;

	if (strcmp(e,"mm")==0) {
		is.open(filename);
		if (!is.good()) return NULL;
		m = new model(is);
		is.close();
	} else { // maybe we have a single mesh file?
		ilist<vec> pts;
		bool spec;
		trimesh *currmesh = loadtrimesh(filename,spec,true);
		if (currmesh==NULL) { // maybe we have an old mm file?
			is.open(filename);
			if (!is.good()) return NULL;

			char buffer[80];
			is >> buffer;
			currmesh = loadtrimesh(buffer,spec,true);
			if (currmesh==NULL) { is.close(); return NULL; }
			loadpts(buffer,pts);

			m = new model(*currmesh,p,0.115,&pts);

			vec *vs = new vec[currmesh->numvertices()];
			ifstream cfp;
			while(!is.eof()) {
				is >> buffer;
				if (is.eof()) break;
				cfp.open(buffer);
				if (cfp.good()) {
					for(i=0;i<currmesh->numvertices();i++)
						cfp >> vs[i];
					cfp.close();
					m->addvector(vs,false);
				} else cout << "could not open " << buffer << endl;
			}
			delete []vs;
			is.close();
		} else {
			m = new model(*currmesh,p,0.115,&pts);
			delete currmesh;
		}
	}
	return m;
}
