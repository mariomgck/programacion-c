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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <loaders/psrender.h>

#define P 0.00001

#define E(a,b) ((b)-(a) < P && (b)-(a) > -P)

psrender::psrender(ostream *is, double scale, double resolution, bool outline,
	bool allsolid,
	double outlinered, double outlinegreen, double outlineblue) {
	res = resolution;
	s = is;
	olr = outlinered;
	olg = outlinegreen;
	olb = outlineblue;
	ol = outline;
	sc = scale;
	solid = allsolid;
}

psrender::~psrender() {
}

void psrender::start(double minx, double maxx, double miny, double maxy) {

	*s << "%!PS-Adobe-3.0 EPSF-3.0" << endl
	  << "%%BoundingBox: -2.0 -2.0 502.0 502.0 " << endl
	  << "gsave\nnewpath" << endl
	  << "250.0 250.0 scale" << endl
	  << "1.0 1.0 translate" << endl;
	*s << "/drawtri {" << endl
	  << "% draw a triangle (r g b x y x y x y)" << endl
	  << "  moveto lineto lineto" << endl
	  << "  setrgbcolor" << endl
	  << "  fill" << endl
	  << "} def" << endl;
	*s << "/drawtrio {" << endl
	  << "% outline a triangle (x y x y x y)" << endl
	  << "  moveto lineto lineto closepath" << endl
	  << "  0 0 0 setrgbcolor" << endl
	  << "  stroke" << endl
	  << "} def" << endl;
	*s << "/drawlinesegm {" << endl
	  << "% draw a line (x y x y r g b)" << endl
	  << "  setrgbcolor" << endl
	  << "  moveto lineto" << endl
	  << "  stroke" << endl
	  << "} def" << endl;
	*s << "/drawtrig {" << endl
	  << "% draw a triangle with linear interpolation of color" << endl
	  << "% (rx ry rc gx gy gc bx by bc "
	  << "minx miny maxx maxy x y x y x y)" << endl
	  << "  gsave" << endl
	  << "  newpath moveto lineto lineto closepath clip" << endl
	  << "  /maxy exch def" << endl
	  << "  /maxx exch def" << endl
	  << "  /miny exch def" << endl
	  << "  /minx exch def" << endl
	  << "  /bc exch def" << endl
	  << "  /by exch def" << endl
	  << "  /bx exch def" << endl
	  << "  /gc exch def" << endl
	  << "  /gy exch def" << endl
	  << "  /gx exch def" << endl
	  << "  /rc exch def" << endl
	  << "  /ry exch def" << endl
	  << "  /rx exch def" << endl
	  << "  /x minx def" << endl
	  << "  /y miny def" << endl
	  << "  /origx " << -res << " minx mul def" << endl
	  << "  /origy " << -res << " miny mul def" << endl
	  << "  /xdis maxx minx sub def" << endl
	  << "  /ydis maxy miny sub def" << endl
	  << "  /count 0 def" << endl
	  << "  /maxcount xdis " << res << " mul cvi 1 add def" << endl
	  << "  /tempy ydis " << res << " mul cvi 1 add def" << endl
	  << "  maxcount  tempy  8" << endl
	  << "  [" << res << " 0 0 " << res << " origx origy]" << endl
	  << "  {drawgradpix} false 3 colorimage" << endl
	  << "  grestore" << endl
	  << "} def" << endl
	  << "" << endl
	  << "/drawgradpix {" << endl
	  << "  /str3 3 string def" << endl
	  << "  str3 0" << endl
	  << "   rx x mul  ry y mul add rc add" << endl
	  << "   255 mul cvi dup 255 gt {pop 255} if" << endl
	  << "   dup 0 lt {pop 0} if" << endl
	  << "  put" << endl
	  << "  str3 1" << endl
	  << "   gx x mul  gy y mul add gc add" << endl
	  << "   255 mul cvi dup 255 gt {pop 255} if" << endl
	  << "   dup 0 lt {pop 0} if" << endl
	  << "  put" << endl
	  << "  str3 2" << endl
	  << "   bx x mul  by y mul add bc add" << endl
	  << "   255 mul cvi dup 255 gt {pop 255} if" << endl
	  << "   dup 0 lt {pop 0} if" << endl
	  << "  put" << endl
	  << "  /count count 1 add def" << endl
	  << "  /x x " << 1.0/res << " add def" << endl
	  << "  count maxcount ge {" << endl
	  << "   /count 0 def" << endl
	  << "   /x minx def" << endl
	  << "   /y y " << 1.0/res << " add def" << endl
	  << "  } if" << endl
	  << "  str3" << endl
	  << "} def" << endl << endl;
	*s << "%%EndProlog" << endl << endl
	  << "%%Page: 1 1" << endl << endl
	  << "gsave" << endl 
	  << "0 setlinewidth" << endl;

}

void psrender::drawtriangle (double x[3], double y[3],
	double r, double g, double b, bool e[3]) {

	*s << r << ' ' << g
	  << ' ' << b << endl
	  << ' ' << sc*x[0] << ' ' << sc*y[0]
	  << ' ' << sc*x[1] << ' ' << sc*y[1]
	  << ' ' << sc*x[2] << ' ' << sc*y[2] << " drawtri" << endl;
	if (e[0] || ol)
		*s << sc*x[0] << ' ' << sc*y[0] << ' '
		  << sc*x[1] << ' ' << sc*y[1] << ' '
		  << olr << ' ' << olg << ' '
		  << olb << " drawlinesegm" << endl;
	if (e[1] || ol)
		*s << sc*x[1] << ' ' << sc*y[1] << ' '
		  << sc*x[2] << ' ' << sc*y[2] << ' '
		  << olr << ' ' << olg << ' '
		  << olb << " drawlinesegm" << endl;
	if (e[2] || ol)
		*s << sc*x[2] << ' ' << sc*y[2] << ' '
		  << sc*x[0] << ' ' << sc*y[0] << ' '
		  << olr << ' ' << olg << ' '
		  << olb << " drawlinesegm" << endl;
}

void psrender::rendertriangle(double x[3], double y[3], bool e[3],
                double cx[3], double cy[3],
                double r[3], double g[3], double b[3]) {
	if (solid) {
		drawtriangle(x,y,(r[0]+r[1]+r[2])/3,
			(g[0]+g[1]+g[2])/3,(b[0]+b[1]+b[2])/3,e);
	} else {
		drawgradstriangle(cx,cy,r,g,b,x,y,e);
	}
}

void psrender::drawgradstriangle(double x[3], double y[3], double r[3],
	double g[3], double b[3], double px[3], double py[3], bool e[3]) {

	double rx,ry,rc,gx,gy,gc,bx,by,bc;
	double maxx,maxy,minx,miny;
	double den[2];

	den[0] = ((x[2]-x[0])*(y[1]-y[0]) + (x[1]-x[0])*(y[0]-y[2]));
	den[1] = ((y[2]-y[0])*(x[1]-x[0]) + (y[1]-y[0])*(x[0]-x[2]));
	maxx = sc*sc*0.001;
	if ((den[0]<maxx && den[0]>-maxx) || (den[1]<maxx && den[1]>-maxx)) {
		drawtriangle(px,py,(r[0]+r[1]+r[2])/3,
			(g[0]+g[1]+g[2])/3,(b[0]+b[1]+b[2])/3,e);
		return;
	}
	rx = ((r[2]-r[0])*(y[1]-y[0]) + (r[1]-r[0])*(y[0]-y[2])) / den[0];
	ry = ((r[2]-r[0])*(x[1]-x[0]) + (r[1]-r[0])*(x[0]-x[2])) / den[1];
	rc = r[0] - sc*rx*x[0] - sc*ry*y[0];
	gx = ((g[2]-g[0])*(y[1]-y[0]) + (g[1]-g[0])*(y[0]-y[2])) / den[0];
	gy = ((g[2]-g[0])*(x[1]-x[0]) + (g[1]-g[0])*(x[0]-x[2])) / den[1];
	gc = g[0] - sc*gx*x[0] - sc*gy*y[0];
	bx = ((b[2]-b[0])*(y[1]-y[0]) + (b[1]-b[0])*(y[0]-y[2])) / den[0];
	by = ((b[2]-b[0])*(x[1]-x[0]) + (b[1]-b[0])*(x[0]-x[2])) / den[1];
	bc = b[0] - sc*bx*x[0] - sc*by*y[0];
	if (rx<0.01 && rx>-0.01 && ry<0.01 && ry>-0.01 &&
	    gx<0.01 && gx>-0.01 && gy<0.01 && gy>-0.01 &&
	    bx<0.01 && bx>-0.01 && by<0.01 && by>-0.01) {
		drawtriangle(px,py,(r[0]+r[1]+r[2])/3,
			(g[0]+g[1]+g[2])/3,(b[0]+b[1]+b[2])/3,e);
		return;
	}
	if (px[1] < px[0]) { maxx = px[0]; minx = px[1]; }
	else { maxx = px[1]; minx = px[0]; }
	if (py[1] < py[0]) { maxy = py[0]; miny = py[1]; }
	else { maxy = py[1]; miny= py[0]; }
	if (maxx < px[2]) maxx = px[2];
	if (minx > px[2]) minx = px[2];
	if (maxy < py[2]) maxy = py[2];
	if (miny > py[2]) miny = py[2];
	*s << rx << ' ' << ry << ' ' << rc << ' ' << gx << ' ' << gy << ' ' 
	  << gc << ' ' << bx << ' ' << by << ' ' << bc << ' ' << endl
	  << "  " << sc*minx << ' ' << sc*miny << ' '
	  << sc*maxx << ' ' << sc*maxy << ' ' <<endl
	  << "  " << sc*px[0] << ' ' << sc*py[0] << ' '
	  << sc*px[1] << ' ' << sc*py[1] << ' '
	  << sc*px[2] << ' ' << sc*py[2] << " drawtrig" << endl;
	if (e[0] || ol)
		*s << sc*px[0] << ' ' << sc*py[0] << ' '
		  << sc*px[1] << ' ' << sc*py[1] << ' '
		  << olr << ' ' << olg << ' '
		  << olb << " drawlinesegm" << endl;
	if (e[1] || ol)
		*s << sc*px[1] << ' ' << sc*py[1] << ' '
		   << sc*px[2] << ' ' << sc*py[2] << ' '
		  << olr << ' ' << olg << ' '
		  << olb << " drawlinesegm" << endl;
	if (e[2] || ol)
		*s << sc*px[2] << ' ' << sc*py[2] << ' '
		  << sc*px[0] << ' ' << sc*py[0] << ' '
		  << olr << ' ' << olg << ' '
		  << olb << " drawlinesegm" << endl;
}

void psrender::end() {
	*s << "grestore" << endl << "grestore" << endl << "showpage" << endl;
}
