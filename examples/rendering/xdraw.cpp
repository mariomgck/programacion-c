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
#include "xfrender.h"

#include <geometry/trimesh.h>
#include <geometry/cull.h>
#include <geometry/model.h>

#include <matrix/vec.h>
#include <morphing/quadric.h>
#include <loaders/mload.h>
#include <loaders/load.h>
#include <image/imagerender.h>
#include <loaders/psrender.h>
#include <warping/warper.h>

#include <fstream>
//#include <unistd.h>
#include <sys/param.h>
#include <unistd.h>
#include <pthread.h>

const double hugenum = 100000;

vec compose(vec minv, vec maxv, bool a, bool b, bool c) {
	vec ret;
	ret = minv;
	for(int i=0;i<DIM;i++) {
		switch(i%3) {
			case 0: if (a) ret[i] = maxv[i];
			break;
			case 1: if (b) ret[i] = maxv[i];
			break;
			case 2: if (c) ret[i] = maxv[i];
		}
	}
	return ret;
}

trimesh *makebox(trimesh *m) {

	trimesh *ret;
	vec minv,maxv;
	m->realextremes(minv,maxv);
	ret = new trimesh(minv,maxv,false,false);
	vec a,b,c;
	ret->addface(compose(minv,maxv,false,false,false),
	             compose(minv,maxv,true,false,false),
	             compose(minv,maxv,false,true,false));
	ret->addface(compose(minv,maxv,true,true,false),
	             compose(minv,maxv,true,false,false),
	             compose(minv,maxv,false,true,false));
	ret->addface(compose(minv,maxv,false,false,true),
	             compose(minv,maxv,true,false,true),
	             compose(minv,maxv,false,true,true));
	ret->addface(compose(minv,maxv,true,true,true),
	             compose(minv,maxv,true,false,true),
	             compose(minv,maxv,false,true,true));

	ret->addface(compose(minv,maxv,false,false,false),
	             compose(minv,maxv,false,false,true),
	             compose(minv,maxv,false,true,false));
	ret->addface(compose(minv,maxv,false,true,true),
	             compose(minv,maxv,false,false,true),
	             compose(minv,maxv,false,true,false));
	ret->addface(compose(minv,maxv,true,false,false),
	             compose(minv,maxv,true,false,true),
	             compose(minv,maxv,true,true,false));
	ret->addface(compose(minv,maxv,true,true,true),
	             compose(minv,maxv,true,false,true),
	             compose(minv,maxv,true,true,false));

	ret->addface(compose(minv,maxv,false,false,false),
	             compose(minv,maxv,true,false,false),
	             compose(minv,maxv,false,false,true));
	ret->addface(compose(minv,maxv,true,false,true),
	             compose(minv,maxv,true,false,false),
	             compose(minv,maxv,false,false,true));
	ret->addface(compose(minv,maxv,false,true,false),
	             compose(minv,maxv,true,true,false),
	             compose(minv,maxv,false,true,true));
	ret->addface(compose(minv,maxv,true,true,true),
	             compose(minv,maxv,true,true,false),
	             compose(minv,maxv,false,true,true));
	return ret;
}

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
	matrix a;
	vec c;
};

void nothing(FL_OBJECT *, long) {
}

class drawstate {
public:
	drawstate() {
		scale = -3; yrot=xrot=zrot=0;
		outline = false;
		split = false;
		color = true;
		red = green = blue = 255;
		olr = olg = olb = 0;
		double s = pow(2,scale);
		conv.a[0][0] = s;
		conv.a[1][1] = s;
		conv.a[2][2] = s;
		r = NULL;
		xfr = NULL;
		m = NULL;
		tm = NULL;
		c = NULL;
		box = NULL;
		boxuse = false;
		tomatch = NULL;
		sliderscale = 1;

		matchform = fl_bgn_form(FL_NO_BOX,500,130);
		fl_add_box(FL_UP_BOX,0,0,500,130,"");
		matchitt = fl_add_input(FL_NORMAL_INPUT,80,25,80,40,
				"iterations:");
			fl_set_object_callback(matchitt,nothing,0);
		fl_set_input(matchitt,"10");
		matchpts = fl_add_input(FL_NORMAL_INPUT,210,25,80,40,
				"points:");
			fl_set_object_callback(matchpts,nothing,0);
		fl_set_input(matchpts,"100");
		matchlambda = fl_add_input(FL_NORMAL_INPUT,340,25,80,40,
				"lambda:");
			fl_set_object_callback(matchlambda,nothing,0);
		fl_set_input(matchlambda,"0.0");
		FL_OBJECT *obj;
		obj = fl_add_button(FL_NORMAL_BUTTON,25,75,
					50,40,"Okay");
		fl_set_object_callback(obj,drawstate::approvematch,(long int)this);
		obj = fl_add_button(FL_NORMAL_BUTTON,425,75,
					50,40,"Cancel");
		fl_set_object_callback(obj,drawstate::cancelmatch,(long int)this);
		fl_end_form();

		ppmform = fl_bgn_form(FL_NO_BOX,350,130);
		fl_add_box(FL_UP_BOX,0,0,350,130,"");
		ppmwidth = fl_add_input(FL_NORMAL_INPUT,80,25,80,40,
				"width:");
		fl_set_input(ppmwidth,"256");
		fl_set_object_callback(ppmwidth,nothing,0);
		ppmheight = fl_add_input(FL_NORMAL_INPUT,210,25,80,40,
				"height:");
		fl_set_input(ppmheight,"256");
		fl_set_object_callback(ppmheight,nothing,0);
		obj = fl_add_button(FL_NORMAL_BUTTON,25,75,50,40,"Okay");
		fl_set_object_callback(obj,drawstate::ppmokay,(long int)this);
		obj = fl_add_button(FL_NORMAL_BUTTON,275,75,50,40,"Cancel");
		fl_set_object_callback(obj,drawstate::ppmcancel,(long int)this);
		fl_end_form();

		progressform = fl_bgn_form(FL_NO_BOX,250,100);
		fl_add_box(FL_UP_BOX,0,0,250,100,"");
		progressslider = fl_add_slider(FL_HOR_FILL_SLIDER,
				25,25,200,50,"");
		fl_deactivate_object(progressslider);
		fl_end_form();

	}
	double scale,yrot,xrot,zrot;
	bool outline;
	bool split,color;
	int red,green,blue,olr,olg,olb;
	affinevecconvert conv;
	trirender *r;
	xfrender *xfr;
	cull *c;
	model *m;
	int nitt,itt;
	int npts;
	float lambda;
	trimesh *tm;
	trimesh *box;
	bool boxuse;
	FL_FORM *matchform;
	FL_OBJECT *matchitt,*matchpts,*progressslider,*matchlambda;
	FL_OBJECT *ppmwidth,*ppmheight;
	FL_FORM *progressform,*ppmform;

	FL_OBJECT *canvas,*xrotslider,*yrotslider,*zrotslider,*scaleslider;
	FL_OBJECT *splitbutton,*outlinebutton,*fastrender,*bb,*colorbutton;
	FL_FORM *mainform;
	ilist<FL_OBJECT *> asliders;
	static bool firstselector[5];
	trimesh *tomatch;
	const char *ppmfn;
	FP sliderscale;

	void static cancelmatch(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->tomatch) delete t->tomatch;
		t->tomatch = NULL;
		fl_activate_form(t->mainform);
		fl_hide_form(t->matchform);
	}

	static void matchcallbackhelper(drawstate *t, int i, int m) {
		t->matchcallback(i,m);
	}

	void matchcallback(int i, int m) {
		itt = i;
	}

	static int warpupdate(XEvent *, void *d) {
		drawstate *t = (drawstate *)d;
		if (t->update!=t->myupdate) {
			cout << t->update << endl;
			t->redraw();
			t->myupdate = t->update;
		}
		if (t->update==-1) {
			fl_activate_form(t->mainform);
			fl_set_idle_callback(0,0);
		}
		return 0;
	}

	static int matchupdate(XEvent *, void *d) {
		drawstate *t = (drawstate *)d;
	
		if (t->itt!=fl_get_slider_value(t->progressslider)) {
			fl_set_slider_value(t->progressslider,t->itt);
			t->setsliders();
			t->redraw();
		}
		if (t->itt==t->nitt) {
			fl_hide_form(t->progressform);
			fl_activate_form(t->mainform);
			fl_set_idle_callback(0,0);
		}
		return 0;
	}

	static void matchrunhelper(drawstate *t) {
		t->matchrun();
	}

	void matchrun() {
		m->matchshape(tomatch,(FP)0.0,nitt,npts,lambda,
			(void (*)(void *,int,int))&drawstate::matchcallbackhelper,
			(void *)this);
		delete tomatch;
		pthread_exit(0);
	}

	static void approvematch(FL_OBJECT *ob, long d) {
		drawstate *t = (drawstate *)d;
		t->nitt = atol(fl_get_input(t->matchitt));
		t->npts = atol(fl_get_input(t->matchpts));
		t->lambda = atof(fl_get_input(t->matchlambda));
		fl_show_form(t->progressform,FL_PLACE_POSITION,FL_FULLBORDER,
			"Matching Progress");
		fl_set_slider_bounds(t->progressslider,0,t->nitt);
		fl_set_slider_value(t->progressslider,0);
		fl_set_idle_callback(drawstate::matchupdate,t);
		fl_hide_form(t->matchform);
		pthread_t th;
		pthread_create(&th,NULL,
		               (void * (*) (void *))&drawstate::matchrunhelper,t);
	}

	static void warp(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(2);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Warp...");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("Object File To Warp to",
			"",(t->firstselector[2]?"*":""),"");
		if (fn==NULL) return;
		t->firstselector[2] = false;
		bool temp;
		t->tomatch = loadtrimesh(fn,temp,true);
		if (t->tomatch==NULL) {
			fl_show_message("could not open",fn,"");
			return;
		}
		fl_use_fselector(4);
		ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Go");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn2 = fl_show_fselector("Parameter File",
			"",(t->firstselector[4]?"*":""),"");
		t->firstselector[4] = false;
		ifstream ps(fn2);
		if (!ps.good()) {
			fl_show_message("could not open",fn2,"");
			delete t->tomatch;
			return;
		}
		ps >> t->wp;
		t->cutlevel = (FP)(-1);
		ps >> t->cutlevel;
		if (t->cutlevel==(FP)-1) {
			fl_show_message("file",fn2,"was not complete through"
				"the cut level parameter");
			delete t->tomatch;
			return;
		}
		t->warpmethod = 0;
		ps >> t->warpmethod;
		ps.close();
		t->update = 0;
		t->myupdate = -1;
		fl_deactivate_form(t->mainform);
		fl_set_idle_callback(drawstate::warpupdate,t);
		pthread_t th;
		pthread_create(&th,NULL,(void *(*) (void *))drawstate::warprunhelper,t);
	}
	
	warperparams wp;
	FP cutlevel;
	int warpmethod;
	volatile int update,myupdate;

	static void buildlevels(trimesh *m, ilist<trimesh *> &l,
				FP mcutlevel, FP cs) {
		ilist<trimesh *> tl;

		tl += new trimesh(*m,1);
		trimesh *rm = new trimesh(*m,cs);
		quadric q(rm);
		int n = q.complexity()/2;
		bool cont;
		mcutlevel *= m->width();
		do {
			while((cont=q.reduce()) && n<q.complexity())
				;
			tl += new trimesh(*rm,1);
			n = q.complexity()/2;
		} while (cont);
		int i;
		for(i=tl.length()-1;i>0 &&
		                     tl[0]->maxdistance(tl[i])>mcutlevel;i--)
			;
		cout << i+1 << '/' << tl.length() << endl;
		for(;i>=0;i--)
			l += tl[i];
		delete rm;
	}

	static void *warprunhelper(drawstate *t) {
		t->warprun();
		return NULL;
	}

	void warprun() {
		ilist<trimesh *> frommesh,tomesh;
		ilist<vec> vlist;
		
		trimesh *holdtm = tm;
		
		cout << "building levels" << endl;
		buildlevels(tm,frommesh,cutlevel,wp.colorscale());
		buildlevels(tomatch,tomesh,cutlevel,wp.colorscale());
		cout << "using warp method: " << warpmethod << endl;
		warper w(&frommesh,&tomesh,vlist,vlist,wp,warpmethod);
		cout << "iterating" << endl;
		while(w.iterate()) {
			tm = w.getA();
			update++;
			while(myupdate!=update)
				;
		}
		cout << "done" << endl;
		delete m;
		m = new model(*tm,wp,cutlevel);
		for(int i=1;i<frommesh.length()-1;i++)
			delete frommesh[i];
		for(int j=0;j<tomesh.length();j++)
			delete tomesh[j];
		update = -1;
		pthread_exit(0);
	}

	static void match(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(2);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Match...");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("Object File to Match",
			"",(t->firstselector[2]?"*":""),"");
		t->firstselector[2] = false;
		bool temp;
		t->tomatch = loadtrimesh(fn,temp,true);
		if (t->tomatch==NULL) {
			fl_show_message("could not open",fn,"");
			return;
		}
		fl_deactivate_form(t->mainform);
		fl_show_form(t->matchform,FL_PLACE_POSITION,FL_FULLBORDER,
			"Matching Parameters");
	}

	static void loadmodel(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(1);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Load");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("Model File to Load",
			"",(t->firstselector[1]?"*":""),"");
		t->firstselector[1] = false;
		t->loadmodel(fn);
	}

	static void savemodel(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(1);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Save");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("File Name to Save",
			"",(t->firstselector[1]?"*":""),"");
		t->firstselector[1] = false;
		t->savemodel(fn);
	}

	static void saveimage(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(3);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Save");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("File to Save",
			"",(t->firstselector[3]?"*.ps *.ppm":""),"");
		t->firstselector[3] = false;
		if (fn) t->saveimage(fn);
	}

	void saveimage(const char *fn) {
		int l = strlen(fn);
		if (l>3 && fn[l-1] == 'm' && fn[l-2] == 'p' && fn[l-3] == 'p'
		    && fn[l-4] == '.') {
			ppmfn = fn;
			fl_deactivate_form(mainform);
			fl_show_form(ppmform,FL_PLACE_POSITION,
				FL_FULLBORDER,"Image Size");
		} else if (l>3 && fn[l-1] == 's' && fn[l-2] == 'p'
		               && fn[l-3] == '.') {
			ofstream outpsfile(fn);
			if (!outpsfile.good()) {
				fl_show_message("File ",fn," could not be opened for writing");
				return;
			}
			psrender *rend = new psrender(&outpsfile,1,200,outline,olr,olg,olb);
			tm->render(*rend,split,&conv,0.25,color);
			delete rend;
			outpsfile.close();
			fl_activate_form(mainform);
		} else {
			fl_show_message("File format for file",fn,
				"is not currenly supported");
		}
	}

	static void ppmcancel(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_activate_form(t->mainform);
		fl_hide_form(t->ppmform);
	}

	static void ppmokay(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_hide_form(t->ppmform);
		imagerender imr(atoi(fl_get_input(t->ppmwidth)),
		                atoi(fl_get_input(t->ppmheight)),
				t->red,t->green,t->blue,t->outline,
		                t->olr,t->olg,t->olb,-1,1,-1,1);
		t->tm->render(imr,t->split,&(t->conv),0.25,t->color);
		ppmreader re(imr.getimage());
		ofstream outf(t->ppmfn,ios::out|ios::binary);
		re.save(outf,255,false,true);
		outf.close();
		fl_activate_form(t->mainform);
	}

	static void flipculling(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t!=NULL && t->c!=NULL) {
			t->c->setliberal(!t->c->getliberal());
			t->redraw();
		}
	}

	static void flipbox(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t!=NULL) {
			t->boxuse = !t->boxuse;
			trimesh *temp = t->box;
			t->box = t->tm;
			t->tm = temp;
			t->redraw();
		}
	}

	void savemodel(const char *fn) {
		if (fn==NULL) return;
		ofstream outf(fn);
		if (!outf.good()) {
			fl_show_message("Could not save to",fn,
				"for unknown reasons");
		} else {
			m->save(outf);
			outf.close();
		}
	}

	void loadmodel(const char *fn) {
		if (fn==NULL) return;
		if (m) delete m;
		warperparams p;
		char olddir[200];
		if (getcwd(olddir,200)!=NULL) {
			char *newdir = new char[strlen(fn)+1];
			int i;
			for(i=strlen(fn)-1;i>=0 && fn[i]!='/';i--);
			if (i!=0) {
				int j;
				for(j=0;j<i;j++)
					newdir[j] = fn[j];
				newdir[j] = 0;
				chdir(newdir);
				m = ::loadmodel(fn+i+1,p);
				chdir(olddir);
			} else m = ::loadmodel(fn,p);
			delete []newdir;
		} else m = ::loadmodel(fn,p);
		if (m) {
			if (boxuse) {
				box = m->getmesh();
				if (tm) delete tm;
				tm = makebox(box);
			} else {	
				tm = m->getmesh();
				if (box) delete box;
				box = makebox(tm);
			}
		} else {
			if (boxuse) {
				if (box) delete box;
			} else {
				if (tm) delete tm;
			}
			tm = NULL;
			box = NULL;
		}
		setsliders();
		recalc();
	}

	void setsliders() {
		if (m==NULL) return;
		int i;
		if (m->numdim()<asliders.length()) {
			for(i=m->numdim();i<asliders.length();i++) {
				fl_delete_object(asliders[i]);
				fl_free_object(asliders[i]);
			}
			asliders.setlength(m->numdim());
		} else if (m->numdim()>asliders.length()) {
			int x = mainform->w-225;
			fl_addto_form(mainform);
			for(i=asliders.length();i<m->numdim();i++) {
				asliders +=
				 fl_add_valslider(FL_HOR_NICE_SLIDER,x,25+30*i,
					200,20,"");
				fl_set_slider_bounds(asliders[i],
					-sliderscale,sliderscale);
				fl_set_slider_step(asliders[i],0.01);
				fl_set_slider_precision(asliders[i],2);
				fl_set_slider_increment(asliders[i],0.01,0.1);
				fl_set_object_gravity(asliders[i],FL_NorthEast,
				     FL_NorthEast);
				fl_set_slider_return(asliders[i],
				     FL_RETURN_END_CHANGED);
				fl_set_object_callback(asliders[i],
				     drawstate::changeslider,(long int)this);
			}
			fl_end_form();
		}
		for(i=0;i<asliders.length();i++)
			fl_set_slider_value(asliders[i],m->getparameter(i));
	}

	static void doublescale(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		t->sliderscale *= 2;
		for(int i=0;i<t->asliders.length();i++)
			fl_set_slider_bounds(t->asliders[i],
				-t->sliderscale,t->sliderscale);
	}
	static void halfscale(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		t->sliderscale /= 2;
		for(int i=0;i<t->asliders.length();i++)
			fl_set_slider_bounds(t->asliders[i],
				-t->sliderscale,t->sliderscale);
	}
	static void splitvecs(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		t->m->split();
		t->setsliders();
	}
	static void zerocenter(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		t->m->zeroorigin();
		t->setsliders();
	}
	static void recenter(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		t->m->setcurrenttoorigin(true);
		t->setsliders();
	}
	static void setcenter(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		t->m->setcurrenttoorigin(false);
		t->setsliders();
	}
	static void even(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		ilist<FP> vs;
		for(int i=0;i<t->m->numdim();i++)
			vs += 1.0/t->m->numdim();
		t->m->setparameters(vs);
		t->setsliders();
		t->redraw();
	}
	static void zeroalpha(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		ilist<FP> vs;
		for(int i=0;i<t->m->numdim();i++)
			vs += 0;
		t->m->setparameters(vs);
		t->setsliders();
		t->redraw();
	}
	static void svd(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		t->m->orthogonalize();
		t->setsliders();
	}
	static void svd2(FL_OBJECT *,long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		t->m->orthogonalize2(NULL);
		t->setsliders();	
	}

	static void changeslider(FL_OBJECT *sl, long d) {
		drawstate *t = (drawstate *)d;
		if (t->m==NULL) return;
		int i;
		for(i=0;i<t->asliders.length();i++)
			if (t->asliders[i]==sl) break;
		if (i==t->asliders.length()) return;
		t->m->setparameter(i,fl_get_slider_value(t->asliders[i]));
		if (t->boxuse) {
			if (t->tm) delete t->tm;
			t->tm = makebox(t->box);
		} else {
			if (t->box) delete t->box;
			t->box = makebox(t->tm);
		}
		t->redraw();
	}

	static void save(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(0);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Save");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("View File to Save","",
			(t->firstselector[0]?"*.view":""),"");
		t->firstselector[0] = false;
		if (fn==NULL) return;
		ofstream fs(fn);
		if (!fs.good()) {
			fl_show_message("could not open",fn,"");
			return;
		}
		fs << "VIEWFF1" << endl;
		fs << t->outline << ' ' << t->split << endl;
		fs << t->c->getliberal() << ' ' << t->boxuse << endl;
		fs << t->red << ' ' << t->green << ' ' << t->blue << endl;
		fs << t->olr << ' ' << t->olg << ' ' << t->olb << endl;
		fs << t->scale << endl;
		fs << t->yrot << ' ' << t->xrot << ' ' << t->zrot << endl;
		fs << t->color << endl;
		fs.close();
	}
			
	static void load(FL_OBJECT *, long d) {
		drawstate *t = (drawstate *)d;
		fl_use_fselector(0);
		FD_FSELECTOR *ffs = fl_get_fselector_fdstruct();
		fl_set_object_label(ffs->ready,"Load");
		fl_fit_object_label(ffs->ready,1,1);
		const char *fn = fl_show_fselector("View File to Load","",
			(t->firstselector[0]?"*.view":""),"");
		t->firstselector[0] = false;
		t->load(fn);
	}

	void load(const char *fn) {
		if (fn==NULL) return;
		ifstream fs(fn);
		if (!fs.good()) {
			fl_show_message(fn,"could not be opened.","");
			return;
		}
		char buffer[80];
		fs >> buffer;
		if (strcmp(buffer,"VIEWFF1")!=0) {
			fl_show_message(fn,"is not a view file.","");
			fs.close();
			return;
		}
		color = true;
		bool fr,b;
		fs >> outline >> split >> fr >> b;
		fs >> red >> green >> blue;
		fs >> olr >> olg >> olb;
		fs >> scale;
		fs >> yrot >> xrot >> zrot;
		fs >> color;
		fs.close();
		fl_set_slider_value(xrotslider,xrot);
		fl_set_slider_value(yrotslider,yrot);
		fl_set_slider_value(zrotslider,zrot);
		fl_set_slider_value(scaleslider,scale);
		fl_set_button(colorbutton,color);
		fl_set_button(splitbutton,split);
		fl_set_button(outlinebutton,outline);
		fl_set_button(fastrender,fr);
		if (fr != c->getliberal()) flipculling(NULL,(long)this);
		if (b != boxuse) flipbox(NULL,(long)this);
		fl_set_button(bb,b);
		if (xfr!=NULL)
			xfr->changeoutline(outline,olr,olg,olb);
		recalc();
	}

	void recalc(void) {
		conv.a = matrix::eye;
		double s = pow(2,scale);
		conv.a[0][0] = s;	
		conv.a[1][1] = s;	
		conv.a[2][2] = s;	
		matrix rot(matrix::eye);
		rot[2][2] = rot[0][0] = cos(M_PI*yrot/180);
		rot[2][0] = -(rot[0][2] = sin(M_PI*yrot/180));
		conv.a = rot*conv.a;
		rot = matrix::eye;
		rot[2][2] = rot[1][1] = cos(M_PI*xrot/180);
		rot[1][2] = -(rot[2][1] = sin(M_PI*xrot/180));
		conv.a = rot*conv.a;
		rot = matrix::eye;
		rot[0][0] = rot[1][1] = cos(M_PI*zrot/180);
		rot[1][0] = -(rot[0][1] = sin(M_PI*zrot/180));
		conv.a = rot*conv.a;
		redraw();
	}

	void redraw() {
		if (tm) tm->render(*r,split,&conv,0.25,color);
	}
		
	static void changescale(FL_OBJECT *ob, long d) {
		drawstate *t = (drawstate *)d;
		t->scale = fl_get_slider_value(ob);
		t->recalc();
	}
	static void changexrot(FL_OBJECT *ob, long d) {
		drawstate *t = (drawstate *)d;
		t->xrot = fl_get_slider_value(ob);
		t->recalc();
	}
	static void changeyrot(FL_OBJECT *ob, long d) {
		drawstate *t = (drawstate *)d;
		t->yrot = fl_get_slider_value(ob);
		t->recalc();
	}
	static void changezrot(FL_OBJECT *ob, long d) {
		drawstate *t = (drawstate *)d;
		t->zrot = fl_get_slider_value(ob);
		t->recalc();
	}
	static void toggleoutline(FL_OBJECT *,long d) {
		drawstate *t = (drawstate *)d;
		t->outline = !(t->outline);
		//cout << t->outline << endl;
		if (t->xfr!=NULL)
			t->xfr->changeoutline(t->outline,t->olr,t->olg,t->olb);
	}
	static void togglesplit(FL_OBJECT *,long d) {
		drawstate *t = (drawstate*)d;
		t->split = !(t->split);
		//cout << t->split << endl;
		t->redraw();	
	}
	static void togglecolor(FL_OBJECT *, long d) {
		drawstate *t = (drawstate*)d;
		t->color = !(t->color);
		t->redraw();
	}
	static int redrawcanvas(FL_OBJECT *ob, int event, FL_Coord,
				FL_Coord, int, void *) {
		if (event==FL_DRAW && ob->u_vdata!=NULL) {
			drawstate *t = (drawstate *)ob->u_vdata;
			if (t!=NULL && t->xfr!=NULL)
				t->xfr->redraw();
		}
		return 0;
	}
};

bool drawstate::firstselector[5] = {true,true,true,true,true};

const char *scalefilter(FL_OBJECT *ob, double v, int prec) {
	static char buffer[20];
	sprintf(buffer,"%.*f",prec,pow(2,v));
	return buffer;
}

int main (int argc, char **argv) {
	fl_initialize(&argc,argv,0,0,0);

	xfrender *xfr;
	drawstate ds;
	int w=1000,h=750;
	int i;

	ds.mainform = fl_bgn_form(FL_NO_BOX,w,h);
	fl_add_box(FL_UP_BOX,0,0,w,h,"");
	FL_OBJECT *obj;
	obj = fl_add_frame(FL_DOWN_FRAME,50,50,w-450,h-130,"");
		fl_set_object_gravity(obj,FL_NorthWest,FL_SouthEast);
	FL_OBJECT *drawarea = fl_add_free(FL_NORMAL_FREE,50,50,w-450,h-130,
		"",ds.redrawcanvas);
		fl_set_object_gravity(drawarea,FL_NorthWest,FL_SouthEast);
	ds.canvas = drawarea;
	obj = fl_add_valslider(FL_VERT_NICE_SLIDER,10,65,30,h-160,"");
		fl_set_slider_bounds(obj,-10,10);
		fl_set_slider_value(obj,ds.scale);
		fl_set_slider_step(obj,0.01);
		fl_set_slider_precision(obj,3);
		fl_set_slider_increment(obj,0.05,1);
		fl_set_object_gravity(obj,FL_NorthWest,FL_SouthWest);
		fl_set_slider_return(obj,FL_RETURN_END_CHANGED);
		fl_set_object_callback(obj,ds.changescale,(long int)&ds);
		fl_set_slider_filter(obj,scalefilter);
		ds.scaleslider = obj;
	obj = fl_add_valslider(FL_VERT_NICE_SLIDER,w-390,65,30,h-160,"");
		fl_set_slider_bounds(obj,-180,180);
		fl_set_slider_value(obj,ds.xrot);
		fl_set_slider_step(obj,1);
		fl_set_slider_precision(obj,0);
		fl_set_slider_increment(obj,1,45);
		fl_set_object_gravity(obj,FL_NorthEast,FL_SouthEast);
		fl_set_slider_return(obj,FL_RETURN_END_CHANGED);
		fl_set_object_callback(obj,ds.changexrot,(long int)&ds);
		ds.xrotslider = obj;
	obj = fl_add_valslider(FL_HOR_NICE_SLIDER,65,15,w-480,20,"");
		fl_set_slider_bounds(obj,-180,180);
		fl_set_slider_value(obj,ds.zrot);
		fl_set_slider_step(obj,1);
		fl_set_slider_precision(obj,0);
		fl_set_slider_increment(obj,1,45);
		fl_set_object_gravity(obj,FL_NorthWest,FL_NorthEast);
		fl_set_slider_return(obj,FL_RETURN_END_CHANGED);
		fl_set_object_callback(obj,ds.changezrot,(long int)&ds);
		ds.zrotslider = obj;
	obj = fl_add_valslider(FL_HOR_NICE_SLIDER,65,h-65,w-480,20,"");
		fl_set_slider_bounds(obj,-180,180);
		fl_set_slider_value(obj,ds.yrot);
		fl_set_slider_step(obj,1);
		fl_set_slider_precision(obj,0);
		fl_set_slider_increment(obj,1,45);
		fl_set_object_gravity(obj,FL_SouthWest,FL_SouthEast);
		fl_set_slider_return(obj,FL_RETURN_END_CHANGED);
		fl_set_object_callback(obj,ds.changeyrot,(long int)&ds);
		ds.yrotslider = obj;
	obj = fl_add_checkbutton(FL_PUSH_BUTTON,w-350,50,100,30,"Fast Render");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.flipculling,(long int)&ds);
		ds.fastrender = obj;
	obj = fl_add_checkbutton(FL_PUSH_BUTTON,w-350,100,100,30,"Bounds Only");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.flipbox,(long int)&ds);
		ds.bb = obj;
	obj = fl_add_checkbutton(FL_PUSH_BUTTON,w-350,150,100,30,"Outlines");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.toggleoutline,(long int)&ds);
		ds.outlinebutton = obj;
	obj = fl_add_checkbutton(FL_PUSH_BUTTON,w-350,200,100,30,"Split");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.togglesplit,(long int)&ds);
		ds.splitbutton = obj;
	obj = fl_add_checkbutton(FL_PUSH_BUTTON,w-350,250,100,30,"Color");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.togglecolor,(long int)&ds);
		fl_set_button(obj,true);
		ds.colorbutton = obj;
	obj = fl_add_button(FL_NORMAL_BUTTON,w-350,300,100,30,"Load View");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.load,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-350,350,100,30,"Save View");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.save,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-350,400,100,30,"Load Model");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.loadmodel,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-350,450,100,30,"Save Image");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.saveimage,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-350,500,100,30,"Save Model");
		fl_set_object_gravity(obj,FL_NorthEast,FL_NorthEast);
		fl_set_object_callback(obj,ds.savemodel,(long int)&ds);


	obj = fl_add_button(FL_NORMAL_BUTTON,w-350,h-80,100,30,"Quit");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		
	obj = fl_add_button(FL_NORMAL_BUTTON,w-225,h-330,85,30,"Zeros");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.zeroalpha,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-125,h-330,85,30,"Even");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.even,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-225,h-280,85,30,"Split");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.splitvecs,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-125,h-280,85,30,"Center");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.recenter,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-225,h-230,85,30,"Set Center");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.setcenter,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-125,h-230,85,30,"Zero Center");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.zerocenter,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-225,h-180,85,30,"2x Scale");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.doublescale,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-125,h-180,85,30,"1/2 Scale");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.halfscale,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-225,h-130,85,30,"Match");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.match,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-125,h-130,85,30,"Warp");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.warp,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-225,h-80,85,30,"SVD");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.svd,(long int)&ds);
	obj = fl_add_button(FL_NORMAL_BUTTON,w-125,h-80,85,30,"SVD 2");
		fl_set_object_gravity(obj,FL_SouthEast,FL_SouthEast);
		fl_set_object_callback(obj,ds.svd2,(long int)&ds);
	
	drawarea->u_vdata = NULL;
	
	fl_end_form();
	
	xfr = new xfrender(drawarea,ds.red,ds.green,ds.blue,ds.outline,
		ds.olr,ds.olg,ds.olb);
	ds.xfr = xfr;
	drawarea->u_vdata = &ds;
/*
	cull *c[argc-2];
	for(i=2;i<argc;i++)
		c[i] = new cull((i==2?xfr:c[i-1]),strtol(argv[i],NULL,10));
	ds.r = (argc<3) ? xfr : c[argc-1];
*/
	ds.r = ds.c = new cull(xfr,2000);
	//ds.r = xfr;
	if (argc>1) ds.loadmodel(argv[1]);
	if (argc>2) ds.load(argv[2]);
	fl_set_form_position(ds.mainform,250,25);
	fl_show_form(ds.mainform,FL_PLACE_POSITION,FL_FULLBORDER,argv[1]);
	fl_do_forms();

/*
	for(i=2;i<argc;i++)
		delete c[i];

*/

	return 0;
}
