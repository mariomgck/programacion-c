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
#include <morphing/quadric.h>
#include <pthread.h>

bootstrap::bootstrap(const warperparams &wp, FP cutlevel,
	trimesh *baseobj) {

	p = wp;
	cl = cutlevel;
	m = new model(*baseobj,p,cl);
	m->setcurrenttoorigin(true);
	m->removelast(1);
	base = baseobj;
}

bootstrap::bootstrap(const warperparams &wp, FP cutlevel,
	trimesh *baseobj, istream &basemodel) {

	p = wp;
	cl = cutlevel;
	m = new model(basemodel);
	base = baseobj;
}

bootstrap::~bootstrap() {
	delete m;
	for(int i=0;i<examples.length();i++)
		delete examples[i];
}

void bootstrap::addobject(trimesh *tm) {
	examples += tm;
}

int bootstrap::numobjects() {
	return examples.length();
}

int bootstrap::nummodeldim() {
	return m->numdim();
}

void buildlevels(trimesh *m, ilist<trimesh *> &l,
		FP cutlevel, FP cs) {

	ilist<trimesh *> tl;

	tl += new trimesh(*m,1);
	trimesh *rm = new trimesh(*m,cs);
	//rm->printlevelstats();
	quadric q(rm);
	int n = q.complexity()/2;
	bool cont;
	cutlevel *= m->width();
	do {
		while((cont=q.reduce()) && n<q.complexity())
			;
		tl += new trimesh(*rm,1);
		n = q.complexity()/2;
	} while (cont);
	FP dum;
	int i;
	for(i=tl.length()-1;i>0&&
		(dum=tl[0]->maxdistance(tl[i]))>cutlevel;i--) {
		delete tl[i];
	}
	for(;i>=0;i--)
		l += tl[i];
	delete rm;
}
void *bootstrap::match(void *voidmp) {
	bootstrap::matchprms *mp = (bootstrap::matchprms *)voidmp;
	ilist<trimesh *> frommesh,tomesh;
	ilist<vec> vlist;
	buildlevels(mp->from,frommesh,mp->cutlevel,mp->wp.colorscale());
	buildlevels(mp->to,tomesh,mp->cutlevel,mp->wp.colorscale());
	warper w(&frommesh,&tomesh,vlist,vlist,mp->wp,0);
	w.run();
	mp->ret = new trimesh(*(w.getA()));
	int i;
	for(i=0;i<frommesh.length();i++)
		delete frommesh[i];
	for(i=0;i<tomesh.length();i++)
		delete tomesh[i];
	return NULL;
}

bool bootstrap::iterate(FP percent, int parallel) {

	int nitems = examples.length();
	trimesh **matches = new trimesh*[nitems];
	volatile matchprms *mp = new matchprms[nitems];
	pthread_t *threads = new pthread_t[nitems];
	int i,j,k;

	int npts = m->veclength()/50;
	if (npts<m->numdim()*3) npts = m->numdim()*3;
	int nitt = m->numdim()*5;
	if (nitt>30) nitt = 30;
	for(i=0;i<nitems;i++) {
		if (m->numdim()>0) m->matchshape(examples(i),0.00001,nitt,npts);
		matches[i] = new trimesh(*(m->getmesh()));
	}
	for(i=0;i<nitems;i++) {
		mp[i].from = matches[i];
		mp[i].to = examples(i);
		const_cast<matchprms *>(mp)[i].wp = p;
		mp[i].cutlevel = cl;
		mp[i].ret = NULL;
	}
	for(i=0;i<nitems;i+=j) {
		if (parallel==1) {
			match(const_cast<matchprms *>(mp)+i);
			j = 1;
			cout << '.' << flush;
		} else {
			for(j=0;j<parallel && i+j<nitems;j++) {
				pthread_create(threads+i+j,NULL,
					bootstrap::match,(void *)(mp+i+j));
			}
			for(k=0;k<j;k++) {
				pthread_join(threads[i+k],NULL);
				cout << '.' << flush;
			}
		}
	}
	cout << ' ' << flush;
	model *newm = new model(*base,p,cl);
	vec *newp = new vec[base->numvertices()];
	for(j=0;j<base->numvertices();j++)
		newp[j] = base->vertexpos(j);
	newm->addvector(newp);
	for(i=0;i<nitems;i++) {
		if (mp[i].ret==NULL) continue;
		for(j=0;j<base->numvertices();j++)
			newp[j] = mp[i].ret->vertexpos(j);
		newm->addvector(newp);
	}
	ilist<FP> as;
	for(i=0;i<newm->numdim();i++)
		as += (FP)1.0/newm->numdim();
	newm->setparameters(as);
	newm->setcurrenttoorigin(true);
	bool retval;
	if (percent<=1.0) {
		FP *scales;
		newm->orthogonalize2(&scales);
		for(i=m->numdim()+1;i<newm->numdim();i++)
			if (scales[i]<scales[m->numdim()]*percent) break;
		cout << '(' << i-m->numdim() << ')' << endl;
		retval = newm->numdim()-i > 0;
		if (retval) newm->removelast(newm->numdim()-i);
		delete []scales;
	} else {
		cout << "(all)" << endl;
		retval = 0;
	}
	delete m;
	m = newm;
	for(i=0;i<nitems;i++)
		if (mp[i].ret != NULL) delete mp[i].ret;
	for(i=0;i<nitems;i++)
		delete matches[i];
	delete []matches;
	delete []mp;
	delete []threads;
	return retval;
}
