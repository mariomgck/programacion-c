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


#include <geometry/face.h>
#include <geometry/vertex.h>
#include <geometry/trimesh.h>
#include <geometry/edge.h>
#include <geometry/trirender.h>
#include <geometry/trirenderer.h>

#include <math.h>
#include "globaldef.h"



//static const straightvecconvert trimesh::svc;
 
trimesh::trimesh(const vec &minpt, const vec &maxpt, const vec &divsize, bool hashon) :
	faces(100), edges(100), vertices(100),
	facegrid(minpt,maxpt,divsize),
	vertexgrid(minpt,maxpt,divsize),
	areasvalid(false), areaface(0), arearem(0), colorscale(1.0), geohashing(hashon)
{
}

static vec ppmminpt(const ppmreader &pr, const FP &scale) {
	vec ret;

	ret[0] = pr.width()*-0.55*scale;
	ret[1] = pr.height()*-0.55*scale;
	ret[2] = (FP)-0.01;
	ret[3] = (FP)-0.1;
	ret[4] = (FP)-0.1;
	ret[5] = (FP)-0.1;
	return ret;
}

static vec ppmmaxpt(const ppmreader &pr, const FP &scale) {
	vec ret;

	ret[0] = pr.width()*0.55*scale;
	ret[1] = pr.height()*0.55*scale;
	ret[2] = (FP)0.01;
	ret[3] = (FP)1.1;
	ret[4] = (FP)1.1;
	ret[5] = (FP)1.1;
	return ret;
}

// this class is just for use for the constructors' median filter
class pixval {
public:
	FP r;
	int x,y;
};

static int pixvalcmp(const void *e1, const void *e2) {
	const pixval *v1=(const pixval *)e1,*v2=(const pixval *)e2;
	if (v1->r<v2->r) return -1;
	if (v1->r>v2->r) return 1;
	return 0;
}

trimesh::trimesh(const ppmreader &pr, int subsamplesize, const FP &scale, bool hashon) :
	faces(100), edges(100), vertices(100),
	facegrid(ppmminpt(pr,scale),ppmmaxpt(pr,scale),1),
	vertexgrid(ppmminpt(pr,scale),ppmmaxpt(pr,scale),1),
	areasvalid(false), areaface(0), arearem(0), geohashing(hashon) {

	int w,h;
	int v1,v2,v3;
	vec p,q;
	int n;
	int x,y;
	vertex *v;
	int *ind;

	colorscale = 1.0;
	w = pr.width()/subsamplesize;
	h = pr.height()/subsamplesize;
	ind = new int[h*w];
	int k=0;
	for(int i=0;i<w;i++) {
		for(int j=0;j<h;j++,k++) {
			p=0;n=0;
			for(int a=0;a<subsamplesize;a++)
				for(int b=0;b<subsamplesize;b++) {
					n++;
					x = i*subsamplesize+a;
					y = j*subsamplesize+b;
					p[0] += (x-pr.width()/2)*scale;
					p[1] += -(y-pr.height()/2)*scale;
					p[3] += pr.getr(x,y);
					p[4] += pr.getg(x,y);
					p[5] += pr.getb(x,y);
				}
				if (n>0) {
					ind[k] = addvertex(p/n);
				} else ind[k] = -1; // this shouldn't happen since I haven't allowed
				// holes in the image
		}
	}
	k=0;
	for(int i=0;i<w-1;i++,k++) {
		for(int j=0;j<h-1;j++,k++) {
			if ((v2=ind[k+h])==-1 || (v3=ind[k+1])==-1) continue;
			if ((v1=ind[k])!=-1) addface(v1,v2,v3);
			if ((v1=ind[k+h+1])!=-1) addface(v1,v2,v3);
		}
	}
	delete []ind;
}

int trimesh::addvertex(const vec &pos) {
	vertex *v = new vertex;
	v->p = pos;
	v->index = vertices.length();
	vertices += v;
	vertexgrid.add(v);
	return v->index;
} 

trimesh::trimesh(const cyberreader &cr,int subsamplesize, bool hashon) : 
	faces(100), edges(100), vertices(100),
	facegrid(cr.minpt()-(cr.maxpt()-cr.minpt())*(FP)0.1,cr.maxpt()+(cr.maxpt()-cr.minpt())*(FP)0.1,0),
	vertexgrid(cr.minpt()-(cr.maxpt()-cr.minpt())*(FP)0.1,cr.maxpt()+(cr.maxpt()-cr.minpt())*(FP)0.1,0),
	areasvalid(false), areaface(0), arearem(0), colorscale(1.0), geohashing(hashon) {

	int w,h;
	int v1,v2,v3;
	vec p,q;
	int n;
	vertex *v;
	int *ind;

	pixval *varr = new pixval[subsamplesize*subsamplesize];
	FP currr,ave;

	cr.scansize(h,w);
	w /= subsamplesize;
	h /= subsamplesize;

	ind = new int[h*w];
	int k =0;
	for(int i=0;i<w;i++) {
		for(int j=0;j<h;j++,k++) {
			p = 0; n=0; ave = 0;
			for(int a=0;a<subsamplesize;a++)
				for(int b=0;b<subsamplesize;b++) {
					/* median filter (on radius values): */
					if(cr.getradius(j*subsamplesize+a,i*subsamplesize+b,currr)) {
						varr[n].r = currr;
						varr[n].x = j*subsamplesize+a;
						varr[n].y = i*subsamplesize+b;
						n++;
						// following just for extremum filter
						ave+=currr;
					}
					/* mean filter:
					if(cr.point(j*subsamplesize+a,i*subsamplesize+b,q)) {
						n++; 
						p += q;
					}
					*/
				}

			/* median filter: */
			if (n==0) ind[k] = -1;
			else {
				qsort(varr,n,sizeof(pixval),pixvalcmp);
				v = new vertex;
				/* median filter: */
				cr.point(varr[n/2].x,varr[n/2].y,v->p);
				if (n&1 == 0) {
					cr.point(varr[n/2+1].x,varr[n/2+1].y,p);
					v->p += p;
					v->p /= 2;
				}
				
				/* extremum filter: 
				ave /= n;
				if (ave-varr[0].r > varr[n-1].r-ave) {
					cr.point(varr[0].x,varr[0].y,v->p);
				} else {
					cr.point(varr[n-1].x,varr[n-1].y,v->p);
				}
				*/

				v->index = vertices.length();
				vertices += v;
				vertexgrid.add(v);
				ind[k] = v->index;
			}

			/* mean filter:
			if (n>0) {
				v = new vertex;
				v->p = p/n;
				v->index = vertices.length();
				vertices += v;
				vertexgrid.add(v);
				ind[k] = v->index;
			} else ind[k] = -1;
			*/
		}
	}

	delete []varr;

	k=0;
	for(int i=0;i<w-1;i++,k++) {
		for(int j=0;j<h-1;j++,k++) {
			if ((v2=ind[k+h])==-1 || (v3=ind[k+1])==-1) continue;
			if ((v1=ind[k])!=-1) addface(v1,v2,v3);
			if ((v1=ind[k+h+1])!=-1) addface(v3,v2,v1);
		}
	}
	k = (w-1)*h; // which I think it should equal anyway, but what the heck.
	for(int j=0;j<h-1;j++) {
		if ((v2=ind[j])==-1 || (v3=ind[k+j+1])==-1) continue;
		if ((v1=ind[k+j])!=-1) addface(v1,v2,v3);
		if ((v1=ind[j+1])!=-1) addface(v3,v2,v1);
	}
	delete []ind;
}

trimesh::~trimesh()
{
	int i;

	for(i=0;i<vertices.length();i++)
		delete vertices[i];
	for(i=0;i<faces.length();i++)
		delete faces[i];
	for(i=0;i<edges.length();i++)
		delete edges[i];
}

ostream &operator<<(ostream &s, const trimesh * const m) {
	m->saveit(s);
	return s;
}

static inline FP boundit(const FP &a) {
	return a<(FP)0.0 ? (FP)0.0 : (a>(FP)1.0 ? (FP)1.0 : a);
}

void trimesh::saveit(ostream &s,fileformat ft) const {

	if (ft==tmff1) {
		int num;
		vertex *v;
		edge *e;
		face *f;
		vec tv;

		s << "TMFF1" << endl; // TMFF1 = TriMesh File Format 1
		s << vertexgrid.getmin() << endl << 
			vertexgrid.getmax() << endl << 
			vertexgrid.getdiv() << endl;
		s << vertices.length() << ' ' << edges.length() << ' ' << faces.length() << endl;
		num = vertices.length();
		for(int i=0;i<num;i++) {
			v = vertices(i);
			tv = v->p;
			for(int j=3;j<DIM;j++)
				tv[j] /= colorscale;
			s << tv << endl;
			s << v->faces.length() << ' ' << v->edges.length() << endl;
			for(int j=0;j<v->faces.length();j++) {
				s << v->faces[j]->index;
				if (j!=v->faces.length()-1) s << ' ';
			}
			s << endl;
			for(int j=0;j<v->edges.length();j++) {
				s << v->edges[j]->index;
				if (j!=v->edges.length()-1) s << ' ';
			}
			s << endl;
		}
		num = edges.length();
		for(int i=0;i<num;i++) {
			e = edges(i);
			s << e->faces.length() << endl;
			for(int j=0;j<e->faces.length();j++) {
				s << e->faces[j]->index;
				if (j!=e->faces.length()-1) s << ' ';
			}
			s << endl;
			s << e->vertices[0]->index << ' ' << e->vertices[1]->index << endl;
		}
		num = faces.length();
		for (int i=0;i<num;i++) {
			f = faces(i);
			s << f->edges[0]->index << ' '
				<< f->edges[1]->index << ' '
				<< f->edges[2]->index << endl;
			s << f->vertices[0]->index << ' '
				<< f->vertices[1]->index << ' '
				<< f->vertices[2]->index << endl;
		}
	} else if (ft==vrml1) { // VRML version 1.0
		int i;
		for(i=0;i<vertices.length()-1;i++) {
			if (vertices(i)->p(3)!=vertices(i+1)->p(3) ||
				vertices(i)->p(4)!=vertices(i+1)->p(4) ||
				vertices(i)->p(5)!=vertices(i+1)->p(5)) break;
		}
		bool solidcolor = (i==vertices.length()-1);

		s << "#VRML V1.0 ascii" << endl;
		s << "Coordinate3 {" << endl;
		s << "\tpoint [" << endl;
		for(i=0;i<vertices.length()-1;i++)
			s << "\t\t" << vertices(i)->p[0] << ' '
			  << vertices(i)->p[1] << ' '
			  << vertices(i)->p[2] << ',' << endl;
		s << "\t\t" << vertices(i)->p[0] << ' '
		  << vertices(i)->p[1] << ' '
		  << vertices(i)->p[2] << endl;
		s << "\t]" << endl;
		s << "}" << endl << endl;
		if (solidcolor) {
			s << "Material {" << endl;
			s << "\tdiffuseColor [" << boundit(vertices(0)->p[3]/colorscale) << ' '
				<< boundit(vertices(0)->p[4]/colorscale) << ' '
				<< boundit(vertices(0)->p[5]/colorscale) << ']' << endl;
			s << "}" << endl;
			s << "MaterialBinding {" << endl;
			s << "\tvalue OVERALL" << endl;
			s << '}' << endl;
		} else {
			s << "Material {" << endl;
			s << "\tdiffuseColor [" << endl;
			for(i=0;i<vertices.length()-1;i++)
				s << "\t\t" << boundit(vertices(i)->p[3]/colorscale) << ' '
				  << boundit(vertices(i)->p[4]/colorscale) << ' '
				  << boundit(vertices(i)->p[5]/colorscale) << ',' << endl;
			s << "\t\t" << boundit(vertices(i)->p[3]/colorscale) << ' '
				  << boundit(vertices(i)->p[4]/colorscale) << ' '
				  << boundit(vertices(i)->p[5]/colorscale) << endl;
			s << "\t]" << endl;
			s << "}" << endl << endl;
			s << "MaterialBinding {" << endl;
			s << "\tvalue PER_VERTEX_INDEXED" << endl;
			s << "}" << endl << endl;
		}
		s << "IndexedFaceSet {" << endl;
		s << "\tcoordIndex [" << endl;
		for(i=0;i<faces.length()-1;i++)
			s << "\t\t" << faces(i)->vertices[0]->index << ", "
			  << faces(i)->vertices[1]->index << ", "
			  << faces(i)->vertices[2]->index << ", "
			  << "-1," << endl;
		s << "\t\t" << faces(i)->vertices[0]->index << ", "
		  << faces(i)->vertices[1]->index << ", "
		  << faces(i)->vertices[2]->index << endl;
		s << "\t]" << endl;
		if (!solidcolor) {
			s << "\tmaterialIndex [" << endl;
			for (i=0;i<faces.length()-1;i++)
				s << "\t\t" << faces(i)->vertices[0]->index << ", "
				  << faces(i)->vertices[1]->index << ", "
				  << faces(i)->vertices[2]->index << ", "
				  << "0," << endl;
			s << "\t\t" << faces(i)->vertices[0]->index << ", "
			  << faces(i)->vertices[1]->index << ", "
			  << faces(i)->vertices[2]->index << endl;
			s << "\t]" << endl;
		}
		s << "}" << endl;
		/*
		s << "Material {" << endl;
		s << "\tdiffuseColor [ 0.0 0.0 1.0 ]" << endl;
		s << "}" << endl;
		s << "MaterialBinding {" << endl;
		s << "\tvalue OVERALL" << endl;
		s << "}" << endl;
		s << "IndexedLineSet {" << endl;
		s << "\tcoordIndex [" << endl;
		for(i=0;i<edges.length();i++)
			if (edges(i)->creaseangle()>0.3)
			s << "\t\t" << edges(i)->vertices[0]->index << ", "
				<< edges(i)->vertices[1]->index << ", -1," << endl;
		s << "\t]" << endl;
		s << "}" << endl;
		s << "Material {" << endl;
		s << "\tdiffuseColor [ 1.0 1.0 1.0 ]" << endl;
		s << "}" << endl;
		s << "MaterialBinding {" << endl;
		s << "\tvalue OVERALL" << endl;
		s << "}" << endl;
		s << "IndexedLineSet {" << endl;
		s << "\tcoordIndex [" << endl;
		for(i=0;i<edges.length();i++)
			if (edges(i)->creaseangle()<=0.3 &&
				edges(i)->creaseangle()>-0.99)
			s << "\t\t" << edges(i)->vertices[0]->index << ", "
				<< edges(i)->vertices[1]->index << ", -1," << endl;
		s << "\t]" << endl;
		s << "}" << endl;
		s << "Material {" << endl;
		s << "\tdiffuseColor [ 1.0 0.0 0.0 ]" << endl;
		s << "}" << endl;
		s << "MaterialBinding {" << endl;
		s << "\tvalue OVERALL" << endl;
		s << "}" << endl;
		s << "IndexedLineSet {" << endl;
		s << "\tcoordIndex [" << endl;
		for(i=0;i<edges.length();i++)
			if (edges(i)->creaseangle()<=-0.99)
			s << "\t\t" << edges(i)->vertices[0]->index << ", "
				<< edges(i)->vertices[1]->index << ", -1," << endl;
		s << "\t]" << endl;
		s << "}" << endl;
		*/
	} else if (ft==oldff) {
		vec tv;

		s << faces.length() << endl;
		for(int i=0;i<faces.length();i++) {
			s << 3 << endl;
			for(int j=0;j<3;j++) {
				tv = faces(i)->vertices[j]->p;
				for(int k=3;k<DIM;k++) tv[k] /= colorscale;
				s << "  " << tv << endl;
			}
			s << endl;
		}
	}
}

static void findbb(const vec &oldmin, const vec &oldmax, const matrix &A, const vec &b,
				   vec &newmin, vec &newmax) {

	int i,mask,j;
	vec curr;

	newmin = hugenum;
	newmax = -hugenum;
	for(i=0;i<1<<DIM;i++) {
		for(j=0,mask=1;j<DIM;j++,mask<<=1) {
			if (i&mask) curr[j] = oldmin(j);
			else curr[j] = oldmax(j);
		}
		curr = A*curr + b;
		for(j=0;j<DIM;j++) {
			if (newmin[j]>curr[j]) newmin[j] = curr[j];
			if (newmax[j]<curr[j]) newmax[j] = curr[j];
		}
	}
}

trimesh::trimesh(istream &s, FP cs, bool hashon, const matrix &A, const vec &b,
	const fileformat ft) :
	faces(100), edges(100), vertices(100),
		areasvalid(false), areaface(0), arearem(0), geohashing(false)
{
	char c1,c2,c3,c4,c5;
	vec themin,themax,thediv;
	vec newmin,newmax;

	if (ft==tmff1) {
		s >> c1 >> c2 >> c3 >> c4 >> c5;
		if (c1!='T' || c2!='M' || c3!='F' || c4 !='F' || c5!='1') {
			themin = 0; themax = 1; thediv = 0;
			// some kind of error thing should go here I suppose
			return;
		}
		s >> themin >> themax >> thediv;
	}

	colorscale = cs;
	
	finishloading(s,A,b,ft);
	geohashing = hashon; // was false from above
	if (hashon) {
		vec minp((FP)hugenum),maxp((FP)-hugenum);
		for(int i=0;i<vertices.length();i++)
			for(int j=0;j<6;j++) {
				if (vertices(i)->p(j)<minp(j)) minp[j]=vertices(i)->p(j);
				if (vertices(i)->p(j)>maxp(j)) maxp[j]=vertices(i)->p(j);
			}
		vertexgrid.initialize(minp,maxp,thediv);
		facegrid.initialize(minp,maxp,thediv);
		for(int i=0;i<faces.length();i++)
			facegrid.add(faces(i));
		for(int i=0;i<vertices.length();i++)
			vertexgrid.add(vertices(i));
	}
	
}

istream &operator>>(istream &s, trimesh *&m) {

	char c1,c2,c3,c4,c5;
	vec minv,maxv,divv;

	s >> c1 >> c2 >> c3 >> c4 >> c5;
	if (c1!='T' || c2!='M' || c3!='F' || c4!='F' || c5!='1')
		return s; // I guess I should implement some error thing -- I'll do that later
	s >> minv >> maxv >> divv;
	m = new trimesh(minv,maxv,divv);
	m->colorscale = 1.0;
	m->finishloading(s);
	return s;
}

void trimesh::finishloading(istream &s, const matrix &A, const vec &b, const fileformat ft) {
	
	int n,n1,n2,n3,i,j;
	vertex *v;
	edge *e;
	face *f;

	if (ft==tmff1) {
		s >> n;
		vertices.setlength(n,true);
		for(i=0;i<n;i++)
			vertices[i] = new vertex;
		s >> n;
		edges.setlength(n,true);
		for(i=0;i<n;i++)
			edges[i] = new edge;
		s >> n;
		faces.setlength(n,true);
		for(i=0;i<n;i++)
			faces[i] = new face;
		n = vertices.length();

		for(i=0;i<n;i++) {

			v = vertices[i];
			v->index = i;
			s >> v->p;

			for(j=3;j<DIM;j++)
				v->p[j] *= colorscale;

			v->p = A*v->p + b;
			s >> n1 >> n2;

			v->faces.setlength(n1);
			for(j=0;j<n1;j++) {
				s >> n3;
				v->faces[j] = faces[n3];
			}
			v->edges.setlength(n2);
			for(j=0;j<n2;j++) {
				s >> n3;
				v->edges[j] = edges[n3];
			}
		}
		n = edges.length();
		for(i=0;i<n;i++) {
			e = edges[i];
			e->index = i;
			s >> n1;
			e->faces.setlength(n1);
			for(j=0;j<n1;j++) {
				s >> n3;
				e->faces[j] = faces[n3];
			}
			s >> n1 >> n2;
			e->vertices[0] = vertices[n1];
			e->vertices[1] = vertices[n2];
		}
		n = faces.length();
		for(i=0;i<n;i++) {
			f = faces[i];
			f->index = i;
			s >> n1 >> n2 >> n3;
			f->edges[0] = edges[n1];
			f->edges[1] = edges[n2];
			f->edges[2] = edges[n3];
			s >> n1 >> n2 >> n3;
			f->vertices[0] = vertices[n1];
			f->vertices[1] = vertices[n2];
			f->vertices[2] = vertices[n3];
		}
	} else if (ft==oldff) {
		int n,t;
		vec v1,v2,v3;

		s >> n;
		for(;n>0;n--) {
			s >> t;
			if (t!=3) break;
			s >> v1 >> v2 >> v3;
			for(int i=3;i<DIM;i++) {
				v1[i] *= colorscale;
				v2[i] *= colorscale;
				v3[i] *= colorscale;
			}
			addface(A*v1+b,A*v2+b,A*v3+b);
		}
	} else return;
	if (geohashing) {
		n = vertices.length();
		for(i=0;i<n;i++)
			vertexgrid.add(vertices[i]);
	}
	n = faces.length();
	for(i=0;i<n;i++) {
		faces[i]->updatebounds();
		if (geohashing)	facegrid.add(faces[i]);
	}
}


static vec colorcorrect(const vec &in, FP cs) {
	vec ret = in;
	for(int i=3;i<DIM;i++) ret[i] *= cs;
	return ret;
}

trimesh::trimesh(const trimesh &from, FP cs, bool hashon, const matrix &A,
	const vec &b) : 
	faces(100), edges(100), vertices(100),
	facegrid(colorcorrect(from.facegrid.getmin(),cs/from.colorscale),
		colorcorrect(from.facegrid.getmax(),cs/from.colorscale),
		colorcorrect(from.facegrid.getdiv(),cs/from.colorscale)),
	vertexgrid(colorcorrect(from.vertexgrid.getmin(),cs/from.colorscale),
		colorcorrect(from.vertexgrid.getmax(),cs/from.colorscale),
		colorcorrect(from.vertexgrid.getdiv(),cs/from.colorscale)),
	areasvalid(false), areaface(0), arearem(0),
	colorscale(cs), geohashing(hashon)
{
	int i,j;
	vertex *v,*ov;
	face *f,*of;
	edge *e,*oe;

	for (i=0;i<from.vertices.length();i++) {
		v = new vertex;
		v->p = A*from.vertices(i)->p+b;
		for(j=3;j<DIM;j++)
			v->p[j] *= cs/from.colorscale;
		v->index = i;
		vertices += v;
	}
	for (i=0;i<from.faces.length();i++) {
		f = new face;
		f->index = i;
		faces += f;
	}
	for (i=0;i<from.edges.length();i++) {
		e = new edge;
		e->index = i;
		edges += e;
	}
	for (i=0;i<from.vertices.length();i++) {
		v = vertices(i);
		ov = from.vertices(i);
		for(j=0;j<ov->faces.length();j++) 
			v->faces += faces[ov->faces(j)->index];
		for(j=0;j<ov->edges.length();j++)
			v->edges += edges[ov->edges(j)->index];
	}
	for(i=0;i<from.faces.length();i++) {
		f = faces(i);
		of = from.faces(i);
		for(j=0;j<3;j++) 
			f->vertices[j] = vertices[of->vertices[j]->index];
		for(j=0;j<3;j++)
			f->edges[j] = edges[of->edges[j]->index];
	}
	for(i=0;i<from.edges.length();i++) {
		e = edges[i];
		oe = from.edges(i);
		for(j=0;j<2;j++)
			e->vertices[j] = vertices[oe->vertices[j]->index];
		for(j=0;j<oe->faces.length();j++)
			e->faces += faces[oe->faces(j)->index];
	}
	geohashing = from.geohashing;
	for(i=0;i<faces.length();i++) {
		faces[i]->updatebounds();
		if (geohashing)	facegrid.add(faces[i]);
	}
	if (geohashing) {
		for(i=0;i<vertices.length();i++)
			vertexgrid.add(vertices[i]);
	}
}

trimesh *trimesh::dup() const {
	trimesh *ret = new trimesh(facegrid.getmin(),facegrid.getmax(),facegrid.getdiv());

	int i,j;
	vertex *v,*ov;
	face *f,*of;
	edge *e,*oe;

	for (i=0;i<vertices.length();i++) {
		v = new vertex;
		v->p = vertices(i)->p;
		v->index = i;
		ret->vertices += v;
	}
	for (i=0;i<faces.length();i++) {
		f = new face;
		f->index = i;
		ret->faces += f;
	}
	for (i=0;i<edges.length();i++) {
		e = new edge;
		e->index = i;
		ret->edges += e;
	}
	for (i=0;i<vertices.length();i++) {
		v = ret->vertices(i);
		ov = vertices(i);
		for(j=0;j<ov->faces.length();j++) 
			v->faces += ret->faces[ov->faces(j)->index];
		for(j=0;j<ov->edges.length();j++)
			v->edges += ret->edges[ov->edges(j)->index];
	}
	for(i=0;i<faces.length();i++) {
		f = ret->faces(i);
		of = faces(i);
		for(j=0;j<3;j++) 
			f->vertices[j] = ret->vertices[of->vertices[j]->index];
		for(j=0;j<3;j++)
			f->edges[j] = ret->edges[of->edges[j]->index];
	}
	for(i=0;i<edges.length();i++) {
		e = ret->edges[i];
		oe = edges(i);
		for(j=0;j<2;j++)
			e->vertices[j] = ret->vertices[oe->vertices[j]->index];
		for(j=0;j<oe->faces.length();j++)
			e->faces += ret->faces[oe->faces(j)->index];
	}
	for(i=0;i<ret->faces.length();i++) {
		ret->faces[i]->updatebounds();
		if (geohashing)	ret->facegrid.add(ret->faces[i]);
	}
	if (geohashing) {
		for(i=0;i<ret->vertices.length();i++)
			ret->vertexgrid.add(ret->vertices[i]);
	}
	ret->geohashing = geohashing;
	ret->colorscale = colorscale;
	return ret;
}


face *trimesh::addface(int v1, int v2, int v3) {

	face *newf;
	vertex *cv, *pv;
	edge *ce;
	int j;

	newf = new face;

	newf->vertices[0] = vertices[v1];
	newf->vertices[1] = vertices[v2];
	newf->vertices[2] = vertices[v3];
	newf->updatebounds();
	vertices[v1]->faces += newf;
	vertices[v2]->faces += newf;
	vertices[v3]->faces += newf;

	for(int i=0;i<3;i++) {
		if (i==2) pv = newf->vertices[0];
		else pv = newf->vertices[i+1];
		cv = newf->vertices[i];
		j = cv->findedge(pv);
		if (j==-1) { // If this edge has never appeared before, we'd better make it
			ce = new edge;
			ce->vertices[0] = cv;
			ce->vertices[1] = pv;
			cv->edges += ce;
			pv->edges += ce;
			ce->index = edges.length();
			edges += ce;
		} else {
			ce = cv->edges[j];
		}
		ce->faces += newf;
		newf->edges[i] = ce;
	}
	newf->index = faces.length();
	faces += newf;
	if (geohashing)	facegrid.add(newf);
	return newf;
}

void trimesh::removeface(int v1, int v2, int v3,
		bool removeedge, bool removept) {

	ilist<face *> *toremove = vertices[v1]->faces.intersection(&(vertices[v2]->faces));
	ilist<face *> *toremove2 = toremove->intersection(&(vertices[v3]->faces));
	delete toremove;
	for(int i=0;i<toremove2->length();i++) {
		removeface(toremove2->nth(i)->index,removeedge,removept);
	}
	delete toremove2;
}
	
void trimesh::removeface(int fnum, bool removeedge, bool removept) {
	face *oldf = faces[fnum];
	ilist<edge *> elist;
	ilist<vertex *> vlist;
	int i,j;

	for(i=0;i<3;i++) {
		for(j=0;j<oldf->edges[i]->faces.length();j++)
			if (oldf->edges[i]->faces[j]==oldf)
				oldf->edges[i]->faces.del(j--);
		if (oldf->edges[i]->faces.length()==0 && removeedge)
			elist += oldf->edges[i];
		for(j=0;j<oldf->vertices[i]->faces.length();j++)
			if (oldf->vertices[i]->faces[j]==oldf)
				oldf->vertices[i]->faces.del(j--);
		if (oldf->vertices[i]->faces.length()==0 && removept)
			vlist += oldf->vertices[i];
	}
	for(i=0;i<elist.length();i++) {
		for(j=0;j<elist[i]->vertices[0]->edges.length();j++)
			if (elist[i]->vertices[0]->edges[j]==elist[i])
				elist[i]->vertices[0]->edges.del(j--);
		for(j=0;j<elist[i]->vertices[1]->edges.length();j++)
			if (elist[i]->vertices[1]->edges[j]==elist[i])
				elist[i]->vertices[1]->edges.del(j--);
		edges.del(elist[i]->index);
		edges[elist[i]->index]->index = elist[i]->index;
		elist[i]->index = -1;
		delete elist[i];
	}
	for(i=0;i<vlist.length();i++) {
		vertices.del(vlist[i]->index);
		vertices[vlist[i]->index]->index = vlist[i]->index;
		vlist[i]->index = -1;
		delete vlist[i];
	}
	faces.del(fnum);
	faces[fnum]->index = fnum;
	oldf->index = -1;
	delete oldf;
}

static int count =0;

void trimesh::removecomponent(int fnum, bool removeedge, bool removept) {
	face *oldf = faces[fnum];
	int i,j;

	if (oldf->gone) return;
	cout << ++count << '/' << faces.length() << endl;
	oldf->gone = true;
	while(1) {
		for(i=0;i<3;i++) {
			for(j=0;j<oldf->edges[i]->faces.length();j++)
				if (oldf->edges[i]->faces[j]!=oldf &&
				    !oldf->edges[i]->faces[j]->gone) break;
			if (j!=oldf->edges[i]->faces.length()) break;
		}
		if (i!=3) {
			removecomponent(oldf->edges[i]->faces[j]->index,removeedge,removept);
		} else break;
	}
	removeface(oldf->index,removeedge,removept);
}
			
static const FP eq_dist2 = (FP)0.00000000000001;

face *trimesh::addface(const vec &p1, const vec &p2, const vec &p3) {

	FP d;
	vertex *cv;
	int v1,v2,v3;
	
	// okay... so duplicating this code three times is pretty lame...
	//   maybe I'll fix the calling parameters later

	cv = closestvertex(p1,d,eq_dist2);
	if (cv==NULL || d>eq_dist2) {
		cv = new vertex;
		cv->p = p1;
		cv->index = vertices.length();
		vertices += cv;
		if (geohashing) vertexgrid.add(cv);
	}
	v1 = cv->index;

	cv = closestvertex(p2,d,eq_dist2);
	if (cv==NULL || d>eq_dist2) {
		cv = new vertex;
		cv->p = p2;
		cv->index = vertices.length();
		vertices += cv;
		if (geohashing) vertexgrid.add(cv);
	}
	v2 = cv->index;

	cv = closestvertex(p3,d,eq_dist2);
	if (cv==NULL || d>eq_dist2) {
		cv = new vertex;
		cv->p = p3;
		cv->index = vertices.length();
		vertices += cv;
		if (geohashing) vertexgrid.add(cv);
	}
	v3 = cv->index;
	return addface(v1,v2,v3);
}

class findcp : public geoiterate<face*> {
public:
	findcp(const vec &sp, const vec &st1, const vec &st2) : p(sp), t1(st1), t2(st2)
		{ fret = NULL; d2 = hugenum; v1=v2=NULL; }

	const vec &p;
	const vec &t1,&t2;
	vec pret;
	FP d2,anglefact;
	face *fret;
	vertex *v1,*v2;

	virtual FP process(face *f) {
		vec tp;
		vertex *tv1,*tv2;
		FP td;

		if (f->closestpoint(p,tp,d2,td,tv1,tv2,t1,t2,anglefact) && td<d2) {
			d2 = td;
			v1 = tv1; v2 = tv2;
			pret = tp;
			fret = f;
		}
		//count++;
		return d2;
	}

	virtual FP initaldist() { return d2; }
};

class findcv : public geoiterate<vertex*> {
public:
	findcv() { ret = NULL; d2 = hugenum; }

	vec p;
	FP d2;
	vertex *ret;

	virtual FP process(vertex *v) {
		FP dis;

		dis = (v->position() - p).len2();
		if (dis<d2) {
			d2 = dis; ret = v;
		}
		return d2;
	}

	virtual FP initaldist() { return d2; }
};

class findcast : public geoiterate<face *> {
public:
	inline findcast(const vec &p, const vec &v, const FP &maxd) {
		this->p = p; this->v = v; this->maxd = maxd; result = false;
	}
	inline ~findcast() {}
	FP maxd;
	vec p,v,pt;
	face *f;
	bool result;

	virtual FP process(face *f) {
		vec tp;
		FP td;

		if (f->rayintersection(p,v,maxd,td,tp) && td<maxd) {
			maxd = td;
			this->f = f;
			result = true;
			pt = tp;
		}
		return maxd;
	}

	virtual FP initaldist() { return maxd; }
};

class findcast3d : public geoiterate<face *> {
public:
	inline findcast3d(const vec &p, const vec &v, const FP &maxd) {
		this->p = p; this->v = v; this->maxd = maxd; result = false; //count=0;
	}
	inline ~findcast3d() {}
	FP maxd;
	vec p,v,pt;
	face *f;
	bool result;
	//int count;

	virtual FP process(face *f) {
		vec tp;
		FP td;

		//count++;
		if (f->rayintersection3d(p,v,maxd,td,tp) && td<maxd) {
			maxd = td;
			this->f = f;
			result = true;
			pt = tp;
		}
		return maxd;
	}

	virtual FP initaldist() { return maxd; }
};

vec trimesh::closestpoint(const vec &p, FP& dist2, face*& fret,
		   vertex* &v1, vertex*& v2, const vec &t1, const vec &t2,
		   const FP &anglefact, const FP &maxdist2) const {

	if (!geohashing) return vec(0);


	findcp fn(p,t1,t2);

	fn.anglefact = anglefact;
	fn.d2 = maxdist2;
	vertex *v = closestvertex(p,dist2,maxdist2);
	if (v) {
		for(int i=0;i<v->faces.length();i++)
			fn.process(v->faces[i]);
	}
	facegrid.search(p,fn);
	//printlevelstats();
	dist2 = fn.d2;
	v1 = fn.v1;
	v2 = fn.v2;
	fret = fn.fret;
	return fn.pret;
}

vertex *trimesh::closestvertex(const vec &p, FP &dist2,
		FP maxdist2) const {

	if (!geohashing) {
		vertex *ret;
		FP r;
		
		dist2 = maxdist2;
		for(int i=0;i<vertices.length();i++) {
			r = (vertices(i)->p-p).len2();
			if (r <= dist2) {
				dist2 = r;
				ret = vertices(i);
			}
		}
		if (dist2 != maxdist2) return ret;
		return NULL;
	}
	findcv fn;

	fn.p = p;
	fn.d2 = maxdist2;
	vertexgrid.search(p,fn);
	dist2 = fn.d2;
	return fn.ret;
}

bool trimesh::castray(const vec &p, const vec &v, vec &ret,
					  face *&fret, FP &retd, const FP &maxd) const {

	if (!geohashing) return false;
	findcast fn(p,v,maxd);
	facegrid.castray(p,v,fn);
	if (fn.result) {
		fret = fn.f;
		ret = fn.pt;
		retd = fn.maxd;
		return true;
	}
	return false;
}

bool trimesh::castray3d(const vec &p, const vec &v, vec &ret,
						face *&fret, FP &retd, const FP &maxd,
						face *hint) const {
	findcast3d fn(p,v,maxd);

	if (hint) {
		for(int i=0;i<3;i++)
			for(int j=0;j<hint->vertices[i]->faces.length();j++)
				fn.process(hint->vertices[i]->faces(j));
	}

	facegrid.castray3d(p,v,fn);
	if (fn.result) {
		fret = fn.f;
		ret = fn.pt;
		retd = fn.maxd;
		return true;
	}
	return false;
}
	

// DO NOT change the order in which edges and faces are removed without also changing
// (and *understanding*) the code in reducer::reduce which relies on the ordering here
// [yes, they probably shouldn't rely on each other that way -- but many things aren't
//  perfect and this code is one of those things]
void trimesh::removeedge(int edgenum, const vec &newpos) {

	edge *ee,*oe,*ne;
	face *f,*cf;
	vertex *tv;
	int i,j,k;
	
	ee = edges[edgenum];
	edges[edges.length()-1]->index = edgenum;
	ee->index = -1;
	edges[edgenum] = edges[edges.length()-1];
	edges-=1;
	if (geohashing) {
		vertexgrid.remove(ee->vertices[0]);
		vertexgrid.remove(ee->vertices[1]);
	}
	tv = ee->vertices[0];
	tv->pointmoved();
	if (geohashing) {
		for(i=0;i<tv->faces.length();i++)
			facegrid.remove(tv->faces[i]);
	}
	tv = ee->vertices[1];
	tv->pointmoved();
	if (geohashing) {
		for(i=0;i<tv->faces.length();i++)
			facegrid.remove(tv->faces[i]);
	}
	for(i=0;i<ee->faces.length();i++) {
		f = ee->faces[i];
		j = f->index;
		faces[j]->index = -1;
		faces[faces.length()-1]->index = j;
		faces[j] = faces[faces.length()-1];
		faces -= 1;
		if (geohashing)	facegrid.remove(f);
		for(j=0;j<3&&f->edges[j]!=ee;j++);
		if (j==0) { // this should *always* be true for normal topologies
			f->edges[0] = f->edges[1];
			f->vertices[0] = f->vertices[1];
			f->edges[1] = f->edges[2];
			f->vertices[1] = f->vertices[2];
		} else if (j==1) {
			f->edges[1] = f->edges[2];
			f->vertices[1] = f->vertices[2];
		}
		for(j=0;j<3;j++) if (f->vertices[j]==ee->vertices[1]) f->vertices[j]=ee->vertices[0];

		oe = f->edges[0];
		ne = f->edges[1];

		tv = oe->vertices[1];
		j = tv->edges.find(oe);
		if (j!=-1) {
			tv->edges[j] = tv->edges[tv->edges.length()-1];
			tv->edges -= 1;
		}
		tv = oe->vertices[0];
		j = tv->edges.find(oe);
		if (j!=-1) {
			tv->edges[j] = tv->edges[tv->edges.length()-1];
			tv->edges -= 1;
		}
		for(j=0;j<oe->faces.length();j++) {
			cf = oe->faces[j];
			if (cf==f) continue;
			ne->faces += cf;
			for (k=0;k<3&&cf->edges[k]!=oe;k++);
			if (k!=3) cf->edges[k] = ne;
		}
		j = ne->faces.find(f);
		if (j!=-1) {
			ne->faces[j] = ne->faces[ne->faces.length()-1];
			ne->faces -= 1;
		}
		j = oe->index;
		edges[j]->index = -1;
		edges[edges.length()-1]->index = j;
		edges[j] = edges[edges.length()-1];
		edges-=1;
		for(j=0;j<3;j++) {
			tv = f->vertices[j];
			k = tv->faces.find(f);
			if (k!=-1) {
				tv->faces[k] = tv->faces[tv->faces.length()-1];
				tv->faces -=1;
			}
		}
		delete oe;
		delete f;
	}
	for(i=0;i<ee->vertices[1]->faces.length();i++) {
		f = ee->vertices[1]->faces[i];
		if (ee->faces.find(f)!=-1) continue;
		for(j=0;j<3&&f->vertices[j]!=ee->vertices[1];j++);
		if (j!=3) f->vertices[j] = ee->vertices[0];
		ee->vertices[0]->faces += f;
	}
	for(i=0;i<ee->vertices[1]->edges.length();i++) {
		ne = ee->vertices[1]->edges[i];
		if (ne==ee) continue;
		if (ne->vertices[0]==ee->vertices[1]) ne->vertices[0] = ee->vertices[0];
		else ne->vertices[1] = ee->vertices[0];
		ee->vertices[0]->edges += ne;
	}
	tv = ee->vertices[0];
	i = tv->edges.find(ee);
	if (i!=-1) {
		tv->edges[i] = tv->edges[tv->edges.length()-1];
		tv->edges -= 1;
	}
	tv->p = newpos;
	if (geohashing) {
		vertexgrid.add(tv);
		for(i=0;i<tv->faces.length();i++)
			facegrid.add(tv->faces[i]);
	}
	tv = ee->vertices[1];
	vertices[vertices.length()-1]->index=tv->index;
	vertices[tv->index] = vertices[vertices.length()-1];
	vertices -= 1;
	tv->index = -1;
	delete tv;
	delete ee;
}


void trimesh::movevertex(int vertexnum, const vec &newpos) {
	vertex *v = vertices[vertexnum];
	face *f;
	vec pmin,pmax;

	if (geohashing) vertexgrid.move(v,v->p,v->p,newpos,newpos);
	vertices[vertexnum]->p = newpos;
	int mi = v->faces.length();
	for(int i=0;i<mi;i++) {
		f = v->faces[i];
		f->tvalid = false;
		pmin = f->min();
		pmax = f->max();
		f->updatebounds();
		if (geohashing)	facegrid.move(f,pmin,pmax,f->min(),f->max());
	}
}

vec trimesh::samplepoint(face *& fret, int numpts, bool checkvalid) {

	if (numpts < 1) numpts=1;
	if ((checkvalid && !areasvalid) ||
		faces.length() != areas.length()) {
		FP tot;
		int i,j,k = faces.length();
		int c;

		totalarea = 0;
		areas.setlength(faces.length());
		for(i=0;i<faces.length();i++) {
			tot = faces[i]->area();
			areas[i] = tot;
			totalarea += tot;
		}
		areas10.setlength((areas.length()+9)/10);
		for(c=i=0;i<areas.length();i+=10,c++) {
			tot = 0;
			for(j=0;j<10;j++) tot += areas[(i+j)%k];
			areas10[c] = tot;
		}
		areas100.setlength(areas10.length()/10);
		for(i=c=0;i<areas10.length()-9;i+=10,c++) {
			tot = 0;
			for(j=0;j<10;j++) tot += areas10[i+j];
			areas100[c] = tot;
		}
		areaface = 0;
		arearem = 0;
		areasvalid = true;
	}

	FP frand,s;
	int size = areas.length();

	frand = drand48()*((FP)3.0*(FP)totalarea/(FP)numpts);
	frand += arearem;
	while(1) {
		if (areas.length()>100 && areaface%100 == 0 &&
		    areaface/100 < areas100.length() &&
		    frand>(s=areas100[areaface/100])) {
			frand -= s;
			areaface += 100;
			if (areaface>=size) areaface -= size;
			continue;
		}
		if (areas.length()>10 && areaface%10 == 0 &&
		    areaface/10 < areas10.length() &&
		    frand>(s=areas10[areaface/10])) {
			frand -= s;
			areaface += 10;
			if (areaface>=size) areaface -= size;
			continue;
		}
		frand -= areas.nth(areaface);
		if (frand < 0) {
			arearem = areas[areaface]+frand;
			fret = faces[areaface];
			return fret->randompoint();
		}
		areaface++;
		if (areaface==areas.length()) areaface = 0;
	}
}

vec trimesh::vertexpos(int vertexnum) const { return vertices(vertexnum)->p; }

image trimesh::cyberwarescan(FP miny, FP maxy, int numlat,int numlon) const {

	image ret(numlat,numlon,false,true);
	if (!geohashing) return ret;
	ret.addchannel();
	ret.addchannel();
	ret.addchannel();
	ret.addchannel();
	int pos=0;
	vec p,v,pt;
	FP d,maxd;
	face *f,*nextf;
	FP maxr=0;

	for(int i=0;i<vertices.length();i++) {
		d = sqrt(vertices(i)->p[2]*vertices(i)->p[2]+
			     vertices(i)->p[0]*vertices(i)->p[0]);
		if (d>maxr) maxr = d;
	}
	maxr *= (FP)1.1;
	nextf=NULL;
	for(int y=0;y<numlon;y++) {
		for(int x=0;x<numlat;x++,pos++) {
			p=0;
			p[1] = (maxy-miny)*x/(FP)(numlat-1) + miny;
			v=0;
			v[2] = cos(2*M_PI*(FP)y/(FP)numlon);
			v[0] = -sin(2*M_PI*(FP)y/(FP)numlon);
			p[2] = -maxr*v[2];
			p[0] = -maxr*v[0];
			maxd = maxr;
			if (castray3d(p,v,pt,f,d,maxd,nextf)) {
				ret.channel(3)[pos] = maxr-d;
				pt = f->addcolor(pt);
				ret.channel(0)[pos] = pt[3];
				ret.channel(1)[pos] = pt[4];
				ret.channel(2)[pos] = pt[5];
				nextf=f;
			} else {
				nextf = NULL;
				ret.channel(0)[pos] = 0;
				ret.channel(1)[pos] = 0;
				ret.channel(2)[pos] = 0;
				ret.channel(3)[pos] = hugenum;
			}
		}
	}
	return ret;
}

void trimesh::colorcomplete() {

	edge *e,*e2;

	int numedges = edges.length();
	bool flip;
	vec d1,d2;
	for(int i=0;i<numedges;i++) {
		e = edges(i);
		if (e->faces.length() > 1) continue;
		int j;
		for(j=i+1;j<numedges;j++) {
			e2 = edges[j];
			flip = false;
			d1 = e->vertices[0]->p - e2->vertices[0]->p;
			d2 = e->vertices[1]->p - e2->vertices[1]->p;
			d1[3] = d1[4] = d1[5] = 0;
			d2[3] = d2[4] = d2[5] = 0;
			if (d1.len2()<0.000001 && d2.len2()<0.000001)
				break;
			flip = true;
			d1 = e->vertices[1]->p - e2->vertices[0]->p;
			d2 = e->vertices[0]->p - e2->vertices[1]->p;
			d1[3] = d1[4] = d1[5] = 0;
			d2[3] = d2[4] = d2[5] = 0;
			if (d1.len2()<0.000001 && d2.len2()<0.000001)
				break;
		}
		if (j==numedges) continue;
		addface(e->vertices[0]->index,e->vertices[1]->index,
				e2->vertices[1]->index);
		if (flip)
			addface(e2->vertices[0]->index,e2->vertices[1]->index,
				e->vertices[0]->index);
		else addface(e2->vertices[0]->index,e2->vertices[1]->index,
				e->vertices[1]->index);
	}
	dlist<vertex *> vl(5);
	for(int i = 0;i<vertices.length();i++) {
		vl.setlength(0);
		int j;
		for(j=0;j<vertices.length();j++) {
			d1 = vertices[i]->p - vertices[j]->p;
			d1[3]=d1[4]=d1[5] = 0;
			if (d1.len2()<0.00001) {
				if (j<i) break;
				vl += vertices[j];
			}
		}
		if (j<i) continue;
		if (vl.length()<3) continue;
		if (vl.length()>3) {
			cout << "CAN'T HANDLE THIS CASE" << endl;
			continue;
		}
		addface(vl[0]->index,vl[1]->index,vl[2]->index);
	}
}

void trimesh::realextremes(vec &minp, vec &maxp) {
	if (vertices.length()==0) return;
	minp = maxp = vertices[0]->p;
	vec v;
	for(int i=1;i<vertices.length();i++) {
		v = vertices[i]->p;
		for(int j=0;j<DIM;j++) {
			if (v[j]>maxp[j]) maxp[j] = v[j];
			if (v[j]<minp[j]) minp[j] = v[j];
		}
	}
}
FP trimesh::width() {
	vec minp,maxp;
	realextremes(minp,maxp);
	maxp -= minp;
	FP r = maxp[0];
	for(int j=1;j<DIM;j++)
		if (maxp[j]>r) r=maxp[j];
	return r;
}
FP trimesh::maxdistance(trimesh *m) {

	FP d,ret = 0;
	int i;
	face *f;
	vertex *v1,*v2;

	for(i=0;i<vertices.length();i++) {
		m->closestpoint(vertices[i]->p,d,f,v1,v2);
		if (d>ret) ret=d;
	}
	for(i=0;i<m->vertices.length();i++) {
		closestpoint(m->vertices[i]->p,d,f,v1,v2);
		if (d>ret) ret=d;
	}
	return sqrt(ret);
}
FP trimesh::mindistance(trimesh *m) {

	FP d,ret = HUGE_VAL;
	int i;
	face *f;
	vertex *v1,*v2;

	for(i=0;i<vertices.length();i++) {
		m->closestpoint(vertices[i]->p,d,f,v1,v2,vec::zero,
			vec::zero,0,ret);
		if (d==0) return 0;
		if (d<ret) ret=d;
	}
	for(i=0;i<m->vertices.length();i++) {
		closestpoint(m->vertices[i]->p,d,f,v1,v2,vec::zero,
			vec::zero,0,ret);
		if (d==0) return 0;
		if (d<ret) ret=d;
	}
	return sqrt(ret);
}
FP trimesh::distance(trimesh *m, int npts) {

	FP d,ret = 0;
	int i;
	face *f;
	vertex *v1,*v2;

	for(i=0;i<npts;i++) {
		m->closestpoint(samplepoint(f,npts/2),d,f,v1,v2);
		ret += sqrt(d);
	}
	return d/(FP)npts;
}

void trimesh::render(trirender &renderer, bool split, vecconvert *vc,
	FP ambient, bool colorized, FP direct) {

	ilist<double> x,y,z,r,g,b;
	ilist<int> f;
	double cx,cy,cz,cr,cg,cb;
	vec pt;

	int i;
	for(i=0;i<vertices.length();i++) {
		pt = vertices[i]->p;
		for(int j=3;j<DIM;j++)
			pt[j] /= colorscale;
		vc->convertpt(pt,cx,cy,cz,cr,cg,cb);
		x += cx;
		y += cy;
		z += cz;
		if (colorized) {
			r += cr;
			g += cg;
			b += cb;
		} else {
			r += 0.8;
			g += 0.8;
			b += 0.8;
		}
	}
	for(i = 0;i<faces.length();i++) {
		f += faces[i]->vertices[0]->index;
		f += faces[i]->vertices[1]->index;
		f += faces[i]->vertices[2]->index;
	}
	trirenderer(renderer,split,x,y,z,r,g,b,f,ambient,direct);
}
		

	
