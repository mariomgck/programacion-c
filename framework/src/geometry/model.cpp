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
#include <geometry/model.h>
#include <warping/warpmesh.h>
#include <string.h>
#include <warping/warper.h>
#include <morphing/quadric.h>
#include <matrix/gmatrix.h>

model::model(const trimesh &basemesh, const warperparams &p, FP cutlevel,
		ilist<vec> *pts, bool split) {
	base = new trimesh(basemesh);
	m += new trimesh(basemesh);
	params = p;
	incomplete = true;
	if (pts!=NULL) upts += *pts;
	uptsname = NULL;
	basecutlevel = cutlevel;
	setbasevecs(split);
}

void model::setbasevecs(bool split) {
	origin = new vec[base->numvertices()];
	for(int j=0;j<base->numvertices();j++)
		origin[j] = 0;
	if (!split) {
		p += new vec[base->numvertices()];
		for(int i=0;i<base->numvertices();i++)
			p[0][i] = base->vertexpos(i);
		alpha += (FP)1;
	} else {
		int i;
		p += new vec[base->numvertices()];
		for(i=0;i<base->numvertices();i++)
			p[0][i] = splitvec(base->vertexpos(i),1);
		alpha += (FP)1;
		p += new vec[base->numvertices()];
		for(i=0;i<base->numvertices();i++)
			p[1][i] = splitvec(base->vertexpos(i),2);
		alpha += (FP)1;
	}
}

model::model(const ilist<char *> &basefiles, const warperparams &p,
		const char * const pts, bool split) {
	int i;
	for(i=0;i<basefiles.length();i++) {
		files += new char[strlen(basefiles(i))+1];
		strcpy(files[i],basefiles(i));
	}
	if (files.length()>0) {
		ifstream ins(files[0]);
		if (ins.good()) {
			base = new trimesh(ins);
			if (base->numvertices()==0) {
				delete base;
				base = NULL;
			}
			ins.close();
		}
		setbasevecs(split);
	} else base = NULL;
	params = p;
	incomplete = false;
	uptsname = strcpy(new char[strlen(pts)+1],pts);
	basecutlevel = (FP)-1.0;
}

model::model(istream &is, bool split) {

	char buffer[80];
	int i,n,n2,j;
	vec v;
	bool origininc=false;

	is >> buffer;
	if (strcmp(buffer,"MMFF")!=0)
		return;
	float vernum;
	is >> vernum;
	if (vernum!=(float)1) {
		if (vernum==(float)1.1) origininc=true;
		else {
			base = NULL;
			return;
		}
	}
	bool ref,complete;
	is >> ref >> complete;
	is >> params;
	is >> basecutlevel;
	if (ref) {
		is >> buffer;
		uptsname = new char[strlen(buffer)+1];
		strcpy(uptsname,buffer);
	} else {
		is >> n;
		for(i=0;i<n;i++) {
			is >> v;
			upts += v;
		}
		uptsname = NULL;
	}
	if (!complete) n = 1;
	else is >> n;
	for(i=0;i<n;i++) {
		if (ref) {
			is >> buffer;
			files += new char[strlen(buffer)+1];
			strcpy(files[files.length()-1],buffer);
		} else {
			m += new trimesh(is);
		}
	}
	is >> n >> n2;
	if (origininc) {
		origin = new vec[n2];
		for(j=0;j<n2;j++)
			is >> origin[j];
	}
	for(i=0;i<n;i++) {
		p += new vec[n2];
		for(j=0;j<n2;j++)
			is >> p[i][j];
	}
	FP a;
	for(i=0;i<n;i++) {
		is >> a;
		alpha += a;
	} 
	if (m.length()>0) {
		base = new trimesh(*(m[0]));
	} else if (files.length()>0) {
		ifstream ins(files[0]);
		if (ins.good()) {
			base = new trimesh(ins);
			if (base->numvertices()==0) {
				delete base;
				base = NULL;
			}
			ins.close();
		}
	} else base = NULL;
	warpbase();
	incomplete = !complete; // how profound!
}

model::~model() {
	int i;
	for(i=0;i<files.length();i++)
		delete []files[i];
	for(i=0;i<m.length();i++)
		delete m[i];
	for(i=0;i<p.length();i++)
		delete []p[i];
	if (base!=NULL) delete base;
	delete []origin;
}

bool model::save(ostream &os, bool ref, bool complete) {
	int i;

	if (complete && incomplete) {
		completelist();
	}
	if (files.length()==0 && ref) return false;
	if (complete && files.length()>m.length()) {
		ifstream ins;

		for(i=0;i<m.length();i++)
			delete m[i];
		m.setlength(0);
		for(i=0;i<files.length();i++) {
			ins.open(files[i]);
			if (!ins.good()) {
				int j;
				for(j=0;j<i;j++)
					delete m[j];
				m.setlength(0);
				return false;
			}
			m += new trimesh(ins);
			ins.close();
		}
	}
	os << "MMFF 1.1" << endl;
	os << ref << " " << complete << endl;
	os << params << endl;
	os << basecutlevel << endl;
	if (ref) {
		os << uptsname;
	} else {
		os << upts.length() << endl;
		for(i=0;i<upts.length();i++)
			os << upts[i] << endl;
	}
	if (!ref) {
		if (complete) {
			os << m.length() << endl;
			for(i=0;i<m.length()-1;i++)
				m[i]->saveit(os);
		}
		m[m.length()-1]->saveit(os);
	} else {
		if (complete) {
			os << files.length() << endl;
			for(i=0;i<files.length()-1;i++)
				os << files[i] << endl;
		}
		os << files[files.length()-1] << endl;
	}
	int j;
	os << p.length() << ' ' << base->numvertices() << endl;
	for(i=0;i<base->numvertices();i++)
		os << origin[i] << endl;
	for(i=0;i<p.length();i++)
		for(j=0;j<base->numvertices();j++)
			os << p[i][j] << endl;
	for(i=0;i<alpha.length();i++)
		os << alpha[i] << endl;
	os << endl;
}

void model::changeparams(const warperparams &p) {
	params = p;
}

void model::addobject(const ilist<char *> &filelist, const char * const pts,
			bool split, int vernum) {
	vec *pnew;
	if (m.length()>=files.length()) {
		if (incomplete) completelist();
		warper w(&m,filelist,upts,pts,params,vernum);

		w.run();
		int s;
		pnew = w.P(s);
	} else {
		warper w(files,filelist,uptsname,pts,params,vernum);
	
		w.run();
		int s;
		pnew = w.P(s);
	}
	addvector(pnew,split);
	delete []pnew;
}

void model::addobject(trimesh *inmesh, FP cutlevel, ilist<vec> *pts, bool split,
	int vernum) {

	ilist<trimesh *> ms;

	ms += inmesh;
	trimesh *rm = new trimesh(*inmesh);
	quadric q(rm);
	int n = q.complexity()/2;
	bool cont;
	do {
		while((cont=q.reduce()) && n<q.complexity())
			;
		if (ms[0]->maxdistance(rm)>cutlevel) break;
		ms += new trimesh(*rm);
		n = q.complexity()/2;
	} while (cont);
	delete rm;
	vec *pnew;
	if (m.length()>=files.length()) {
		if (incomplete) completelist();
		warper w(&m,&ms,upts,*pts,params,vernum);
		w.run();
		int s;
		pnew = w.P(s);
	} else {
		warper w(files,&ms,uptsname,*pts,params,vernum);
		w.run();
		int s;
		pnew = w.P(s);
	}
	addvector(pnew,split);
	delete []pnew;
	for(n=0;n<ms.length();n++)
		delete ms[n];
}
	
void model::completelist() {
	trimesh *rm = new trimesh(*(m[0]));
	quadric q(rm);
	int n = q.complexity()/2;
	bool cont;
	do {
		while((cont=q.reduce()) && n<q.complexity())
			;
		if (m[0]->maxdistance(rm)>basecutlevel) break;
		m += new trimesh(*rm);
		n = q.complexity()/2;
	} while (cont);
	delete rm;
	incomplete = false;
}

void model::addvector(const ilist<vec> &v, bool split) {

	if (!split) {
		p += new vec[v.length()];
		for(int i=0;i<v.length();i++)
			p[p.length()-1][i] = v(i)-origin[i];
		alpha += (FP)0;
	} else {
		p += new vec[v.length()];
		for(int i=0;i<v.length();i++)
			p[p.length()-1][i] = splitvec(v(i)-origin[i],1);
		alpha += (FP)0;
		p += new vec[v.length()];
		for(int i=0;i<v.length();i++)
			p[p.length()-1][i] = splitvec(v(i)-origin[i],2);
		alpha += (FP)0;
	}
}

void model::addvector(vec *v, bool split) {

	if (m.length()==0) return;
	if (!split) {
		p += new vec[base->numvertices()];
		for(int i=0;i<base->numvertices();i++) 
			p[p.length()-1][i] = v[i];
		alpha += (FP)0;
	} else {
		int i;
		p += new vec[base->numvertices()];
		for(i=0;i<base->numvertices();i++) 
			p[p.length()-1][i] = splitvec(v[i],1);
		alpha += (FP)0;
		p += new vec[base->numvertices()];
		for(i=0;i<base->numvertices();i++) 
			p[p.length()-1][i] = splitvec(v[i],2);
		alpha += (FP)0;
	}
}

int model::numdim() {
	return alpha.length();
}

int model::veclength() {
	if (base==NULL) return -1;
	return base->numvertices();
}

void model::setparameter(int i,FP v) {
	if (i>=0&&i<alpha.length()) alpha[i] = v;
	else cerr << "bad index, " << i << endl;
	warpbase();
}

void model::setparameters(const ilist<FP> &vs) {
	for(int i=0;i<vs.length()&&i<alpha.length();i++)
		alpha[i] = vs(i);
	warpbase();
}

FP model::getparameter(int i) {
	if (i>=0&&i<alpha.length()) return alpha[i];
	else return (FP)0;
}

void model::getparameters(ilist<FP> &vs,bool append) {
	if (append)
		for(int i=0;i<alpha.length();i++)
			vs += alpha[i];
	else
		for(int i=0;i<alpha.length()&&i<vs.length();i++)
			vs[i] = alpha[i];
}

void model::getvector(int i, ilist<vec> &v, bool append) {
	if (base==NULL) return;
	if (append)
		for(int j=0;j<base->numvertices();j++)
			v += p[i][j];
	else
		for(int j=0;j<base->numvertices() && j<v.length();j++)
			v[j] = p[i][j];
}

void model::getvector(int i, vec *v) {
	if (base==NULL) return;
	for(int j=0;j<base->numvertices();j++)
		v[j] = p[i][j];
}

void model::getorigin(ilist<vec> &v, bool append) {
	if (base==NULL) return;
	if (append)
		for(int j=0;j<base->numvertices();j++)
			v += origin[j];
	else
		for(int j=0;j<base->numvertices() && j<v.length();j++)
			v[j] = origin[j];
}

void model::getorigin(vec *v) {
	if (base==NULL) return;
	for(int j=0;j<base->numvertices();j++)
		v[j] = origin[j];
}

trimesh *model::getmesh() {
	return base;
}

void model::warpbase() {

	if (base==NULL) return;
	vec pos;

	for(int i=0;i<base->numvertices();i++) {
		pos = origin[i];
		for(int j=0;j<alpha.length() && j<p.length();j++)
			pos += p[j][i]*alpha[j];
		base->movevertex(i,pos);
	}
}

void model::matchshape(trimesh *shape, FP pchange, int nitt, int npts,
	FP lambda,
	void (*cbfn)(void *, int , int ), void *cbdata) {
	int i;
	FP *newalpha;

	for(i=0;i<alpha.length();i++)
		alpha[i] = (FP)1.0/alpha.length();
	warpbase();

	gmatrix P(alpha.length(),alpha.length());
	gmatrix Ptemp(alpha.length(),DIM);
	gmatrix Qtemp(DIM,1);
	gmatrix Q(alpha.length(),1);
	vec pt,q;
	FP realanswer;
	FP qconst;
	FP a[3],d2;
	face *f,*f2;
	vertex *v1,*v2;
	int j,k;

	FP change,newv;

	if (cbfn!=NULL) cbfn(cbdata,0,nitt);
	for(int itt=0;itt<nitt;itt++) {
		P = 0;
		Q = 0;
		qconst = (FP)npts*lambda/alpha.length();
		realanswer = 0;
		gmatrix ba(&(alpha[0]),1,alpha.length());
		for(i=0;i<npts;i++) {
			q = shape->samplepoint(f,npts);
			pt = base->closestpoint(q,d2,f,v1,v2);
			realanswer += (q-pt).len2();
			f->project(pt,a);
			for(j=0;j<alpha.length();j++)
				for(k=0;k<DIM;k++)
					Ptemp[j][k] =
					   a[0]*p[j][f->vertexindex(0)][k] +
					   a[1]*p[j][f->vertexindex(1)][k] +
					   a[2]*p[j][f->vertexindex(2)][k];
			P += gmatrix(Ptemp,Ptemp);
			for(k=0;k<DIM;k++) Qtemp[k][0] =
				q[k] - a[0]*origin[f->vertexindex(0)][k]
				     - a[1]*origin[f->vertexindex(1)][k]
				     - a[2]*origin[f->vertexindex(2)][k];
			Q += Ptemp*Qtemp;
			qconst += q.len2();
		}
		for(i=0;i<npts;i++) {
			pt = base->samplepoint(f,npts);
			q = shape->closestpoint(pt,d2,f2,v1,v2);
			realanswer += (q-pt).len2();
			f->project(pt,a);
			for(j=0;j<alpha.length();j++)
				for(k=0;k<DIM;k++)
					Ptemp[j][k] =
					   a[0]*p[j][f->vertexindex(0)][k] +
					   a[1]*p[j][f->vertexindex(1)][k] +
					   a[2]*p[j][f->vertexindex(2)][k];
			P += gmatrix(Ptemp,Ptemp);
			for(k=0;k<DIM;k++) Qtemp[k][0] =
				q[k] - a[0]*origin[f->vertexindex(0)][k]
				     - a[1]*origin[f->vertexindex(1)][k]
				     - a[2]*origin[f->vertexindex(2)][k];
			Q += Ptemp*Qtemp;
			qconst += q.len2();
		}
		P+=gmatrix(alpha.length(),alpha.length(),lambda*npts,true);
		Q+=gmatrix(alpha.length(),1,lambda*npts/alpha.length(),false);
		change =  qconst + (ba*gmatrix(P,ba) - 2*ba*Q)[0][0];
		//cout << "fit before: " << change << endl;
		newalpha = P.solve(Q[0]);
		gmatrix aa(newalpha,1,alpha.length());
		newv = qconst + (aa*gmatrix(P,aa) - 2*aa*Q)[0][0];
		change -= newv;
		//cout << "fit after: " << newv << endl;
		for(i=0;i<alpha.length();i++) {
			//cout << newalpha[i] << endl;
			alpha[i] = newalpha[i];
		}
		//cout << endl;
		warpbase();
		delete []newalpha;
		if (cbfn!=NULL) cbfn(cbdata,itt+1,nitt);
		if (change < newv*pchange) // just fake the rest!
			for(itt++;itt<nitt;itt++)
				if (cbfn!=NULL) cbfn(cbdata,itt+1,nitt);
	}
}

void model::setpmatrix(gmatrix &m) {
	int i,j,k;

	for(i=0;i<alpha.length();i++) {
		for(j=0;j<base->numvertices();j++) {
			for(k=0;k<DIM;k++) m[k+j*DIM][i] = p[i][j][k];
		}
	}
}

void model::pfrommatrix(const gmatrix &m, FP *scales) {
	
	int i,j,k;
	vec v;

	for(i=0;i<alpha.length();i++) {
		for(j=0;j<base->numvertices();j++) {
			for(k=0;k<DIM;k++) v[k] = m(k+j*DIM,i);
			p[i][j] = scales ? v*scales[i] : v;
		}
	}
}
	
void model::zeroorigin() {
	int i,j;
	gmatrix m(DIM*base->numvertices(),alpha.length());
	setpmatrix(m);
	gmatrix outp(alpha.length(),alpha.length());
	m.inner(m,outp);
	gmatrix c(alpha.length(),1);
	gmatrix orig(DIM*base->numvertices(),1);
	for(i=0;i<base->numvertices();i++)
		for(j=0;j<DIM;j++) orig[i*DIM+j][0] = origin[i][j];
	m.inner(orig,c);
	FP *ans;
	ans = outp.solve(c[0]);
	for(i=0;i<alpha.length();i++) alpha[i] += ans[i];
	delete []ans;
	for(i=0;i<base->numvertices();i++) origin[i] = 0;
	// in theory since origin was made as a linear combination
	// of the basis vectors, we don't have to rewarp the base mesh
}

void model::setcurrenttoorigin(bool changebasis) {
	int i;
	if (changebasis) {
		for(i=0;i<base->numvertices();i++)
			for(int j=0;j<alpha.length();j++)
				p[j][i] += origin[i] - base->vertexpos(i);
	}
	for(i=0;i<alpha.length();i++) alpha[i] = 0;
	for(i=0;i<base->numvertices();i++) origin[i] = base->vertexpos(i);
}

void model::removedim(int d) {
	delete [](p[d]);
	for(int i=d;i<p.length()-1;i++) {
		p[i] = p[i+1];
		alpha[i] = alpha[i+1];
	}
	alpha -= 1;
	p -= 1;
	warpbase();
}

void model::removelast(int d) {
	for(int i=p.length()-d;i<p.length();i++)
		delete [](p(i));
	alpha -= d;
	p -= d;
	warpbase();
}

void model::split(int v) {
	if (v==-1) {
		int d = alpha.length();
		for(int i=0;i<d;i++) split(i);
	} else {
		vec *newp = new vec[base->numvertices()];
		vec *oldp = new vec[base->numvertices()];
		bool change1 = false, change0 = false;
		for(int i=0;i<base->numvertices();i++) {
			newp[i] = splitvec(p[v][i],0);
			oldp[i] = splitvec(p[v][i],1);
			if (newp[i]!=0) change1=true;
			if (oldp[i]!=0) change0=true;
		}
		if (!change1 || !change0) {
			delete []newp;
			delete []oldp;
			return;
		}
		delete []p[v];
		p[v] = oldp;
		p += newp;
		alpha += alpha[v];
	}
}

void model::matchparms(vec *v) {
	FP *n = new FP[DIM*(base->numvertices()+base->numfaces())];
	FP newval;
	int i,j,k,c;

	for(i=0;i<base->numvertices();i++)
		for(j=0;j<DIM;j++) n[i*DIM+j] = (v[i]-origin[i])(j)*
						base->vertexnum(i)->sqrtarea();
	for(i=0;i<base->numfaces();i++)
		for(j=0;j<DIM;j++) n[(i+base->numvertices())*DIM+j] = 
			((v[base->facenum(i)->vertexindex(0)]-
			        origin[base->facenum(i)->vertexindex(0)])(j)+
			 (v[base->facenum(i)->vertexindex(1)]-
			        origin[base->facenum(i)->vertexindex(1)])(j)+
			 (v[base->facenum(i)->vertexindex(2)]-
			        origin[base->facenum(i)->vertexindex(2)])(j))*
				sqrt(base->facenum(i)->area());

	for(k=0;k<alpha.length();k++) {
		newval = 0;
		c = 0;
		for(i=0;i<base->numvertices();i++)
			for(j=0;j<DIM;j++) newval += n[c++]*p[k][i][j]*
						base->vertexnum(i)->sqrtarea();
		for(i=0;i<base->numfaces();i++)
			for(j=0;j<DIM;j++) newval += n[c++]*
				(p[k][base->facenum(i)->vertexindex(0)][j]+
				 p[k][base->facenum(i)->vertexindex(1)][j]+
				 p[k][base->facenum(i)->vertexindex(2)][j])*
					sqrt(base->facenum(i)->area());
		alpha[k] = newval;
	}
	delete []n;
	warpbase();
}


void model::orthogonalize2(FP **sret) {
	gmatrix m(DIM*(base->numvertices()+base->numfaces()),alpha.length());
	int i,j,k;
	const FP sqrt12 = sqrt(12.0);

	for(i=0;i<alpha.length();i++)
		for(j=0;j<base->numvertices();j++)
			for(k=0;k<DIM;k++)
				m[k+j*DIM][i] = p[i][j][k]*
					base->vertexnum(j)->sqrtarea()/sqrt12;
	for(i=0;i<alpha.length();i++)
		for(j=0;j<base->numfaces();j++)
			for(k=0;k<DIM;k++)
				m[(base->numvertices()+j)*DIM+k][i] = 
				 (p[i][base->facenum(j)->vertexindex(0)][k] +
				  p[i][base->facenum(j)->vertexindex(1)][k] +
				  p[i][base->facenum(j)->vertexindex(2)][k])*
					 sqrt(base->facenum(j)->area())/sqrt12;
	gmatrix amult(alpha.length(),alpha.length());
	FP *scales = new FP[alpha.length()];
	m.svd(m,amult,scales);

/*
	ofstream of("Data/v.mat");
	amult.niceprint(of);
	of.close();
	ofstream of2("Data/d.mat");
	for(int i=0;i<alpha.length();i++)
		of2 << scales[i] << endl;
	of2.close();
*/

	gmatrix mt(DIM*base->numvertices(),alpha.length());
	setpmatrix(mt);
	mt *= amult;
	mt *= 1.0/sqrt(alpha.length());
	pfrommatrix(mt);

/*	NOT NEEDED (as determined on 8/25/99)
	// scale axes by standard deviation
	// (no need to scale m since it is discarded after this)
	m *= amult;
	FP scalefact;
	for(i=0;i<alpha.length();i++) {
		scalefact = 0;
		for(j=0;j<DIM*(base->numvertices()+base->numfaces());j++)
			scalefact += m[j][i]*m[j][i];
		scalefact = sqrt(scalefact*alpha.length())/scales[i];
		for(j=0;j<base->numvertices();j++)
			p[i][j] /= scalefact;
	}
*/
	
	// spin alphas:
	gmatrix av(1,alpha.length());
	for(i=0;i<alpha.length();i++) av[0][i] = alpha[i];
	av = av*amult;
	for(i=0;i<alpha.length();i++) alpha[i] = av[0][i];

	// sort dimensions according to singular values
	// (simple selection sort)
	FP best; int bi;
	vec *tv;
	for(i=0;i<alpha.length()-1;i++) {
		bi = i; best = scales[i];
		for(j=i+1;j<alpha.length();j++) {
			if (scales[j]>best) {
				best = scales[j];
				bi = j;
			}
		}
		if (bi!=i) {
			scales[bi] = scales[i];
			scales[i] = best;
			best = alpha[bi];
			alpha[bi] = alpha[i];
			alpha[i] = best;
			tv = p[bi];
			p[bi] = p[i];
			p[i] = tv;
		}
	}
	//for(i=0;i<alpha.length();i++)
	//	cout << scales[i] << endl;
	if (sret!=NULL) *sret = scales;
	else delete []scales;
}

void model::orthogonalize(void) {
	gmatrix m(DIM*base->numvertices(),alpha.length());
	int i,j;

	// orthogonalize (and then scale by singular values):
	setpmatrix(m);
	gmatrix amult(alpha.length(),alpha.length());
	FP *scales = new FP[alpha.length()];
	m.svd(m,amult,scales);
	pfrommatrix(m,scales);

	// spin alphas:
	gmatrix av(1,alpha.length());
	for(i=0;i<alpha.length();i++) av[0][i] = alpha[i];
	av = av*amult;
	for(i=0;i<alpha.length();i++) alpha[i] = av[0][i];

	// sort dimensions according to singular values
	// (simple selection sort)
	FP best; int bi;
	vec *tv;
	for(i=0;i<alpha.length()-1;i++) {
		bi = i; best = scales[i];
		for(j=i+1;j<alpha.length();j++) {
			if (scales[j]>best) {
				best = scales[i];
				bi = j;
			}
		}
		if (bi!=i) {
			scales[bi] = scales[i];
			scales[i] = best;
			best = alpha[bi];
			alpha[bi] = alpha[i];
			alpha[i] = best;
			tv = p[bi];
			p[bi] = p[i];
			p[i] = tv;
		}
	}
	for(i=0;i<alpha.length();i++)
		cout << scales[i] << endl;
	delete []scales;
}
