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
#include <loaders/vtkreader.h>
#include <containers/ilist.h>
#include <matrix/vec.h>

#include <stdio.h>
#include <fstream>
#include <string.h>

#define DEPTH 0.001

using namespace std;

vtkreader::vtkreader(const char *filename, bool hashon) {
	FILE *fp; // yes, we are doing this the old-fashioned way!
	char buffer[100];
	bool binary;

	fp = fopen(filename,"r");
	if (fp==NULL) { ret = NULL; return; }
	fscanf(fp,"%90[^\n0123456789]",buffer);
	if (strcmp(buffer,"# vtk DataFile Version ")) {
		ret = NULL;
		fclose(fp);
		return;
	}
	fscanf(fp,"%*d.%*d%*[^\n]");
	fscanf(fp,"\n%*[^\n]");
	fscanf(fp,"\n");
	fscanf(fp,"%[^\n]\n",buffer);
	if (strcmp(buffer,"BINARY")==0)
		binary = true;
	else if (strcmp(buffer,"ASCII")==0)
		binary = false;
	else {
		ret = NULL;
		fclose(fp);
		return;
	}
	fscanf(fp,"%[^\n]\n",buffer);
	if (strcmp(buffer,"DATASET POLYDATA")) {
		ret = NULL;
		fclose(fp);
		return;
	}
	
	// so I'm making a couple of bad assumptions:
	//  1. the file will always have the POINTS section first
	//  2. when adding vertices to the mesh, they will get sequential #'s
	//  3. the file gives the coordinates in floats (for BINARY)

	vec vmin((FP)HUGE_VAL);
	vec vmax(-(FP)HUGE_VAL);
	vec div(1);
	vec vpos;

	fscanf(fp,"%s",buffer);
	if (strcmp(buffer,"POINTS")) {
		ret = NULL;
		fclose(fp);
		return;
	}
	int npts;
	fscanf(fp,"%d%s%*[\n]",&npts,buffer);
	fscanf(fp,"\n");
	if (binary && strcmp(buffer,"float")) {
		ret = NULL;
		fclose(fp);
		return;
	}
	int i,j;
	float rv[3];
	ilist<vec> vlist;

	for(i=0;i<npts;i++) {
		if (binary) fread(rv,sizeof(float),3,fp);
		else fscanf(fp,"%f%f%f",rv,rv+1,rv+2);
		vpos[0] = rv[0];
		vpos[1] = rv[1];
		vpos[2] = rv[2];
		vpos[3] = vpos[4] = vpos[5] = 0.8;
		vmin = min(vmin,vpos);
		vmax = max(vmax,vpos);
		vlist += vpos;
		printf("%d/%d\r",i,npts);
		fflush(stdout);
	}
	printf ("done with vertices\n");
	vmin[3]=vmin[4]=vmin[5] = 0;
	vmax[3]=vmax[4]=vmax[5] = 1;
	ret = new trimesh(vmin,vmax,div,hashon);
	for(i=0;i<npts;i++)
		ret->addvertex(vlist[i]);
	
	int npoly,size,nvert;
	while(!feof(fp)) {
		fscanf(fp,"%s%d%d%*[^\n]",buffer,&npoly,&size);
		if (feof(fp)) break;
		fscanf(fp,"\n");
		if (strcmp(buffer,"POLYGONS")==0) {
			for(i=0;i<npoly;i++) {
				printf ("P:%d/%d\r",i,npoly);
				fflush(stdout);
				if (binary) fread(&nvert,sizeof(int),1,fp);
				else fscanf(fp,"%d",&nvert);
				int base,v1,v2;
				if (binary) fread(&base,sizeof(int),1,fp);
				else fscanf(fp,"%d",&base);
				if (binary) fread(&v1,sizeof(int),1,fp);
				else fscanf(fp,"%d",&v1);
				for(j=2;j<nvert;j++) {
					if (binary) fread(&v2,sizeof(int),1,fp);
					else fscanf(fp,"%d",&v2);
					ret->addface(base,v1,v2);
					v1 = v2;
				}
			}
			printf ("done P\n");
		} else if (strcmp(buffer,"TRIANGLE_STRIPS")==0) {
			for(i=0;i<npoly;i++) {
				printf ("T:%d/%d\r",i,npoly);
				fflush(stdout);
				if (binary) fread(&nvert,sizeof(int),1,fp);
				else fscanf(fp,"%d",&nvert);
				int v1,v2,v3;
				if (binary) fread(&v1,sizeof(int),1,fp);
				else fscanf(fp,"%d",&v1);
				if (binary) fread(&v2,sizeof(int),1,fp);
				else fscanf(fp,"%d",&v2);
				for(j=2;j<nvert;j++) {
					if (binary) fread(&v3,sizeof(int),1,fp);
					else fscanf(fp,"%d",&v3);
					ret->addface(v1,v2,v3);
					v1 = v2;
					v2 = v3;
				}
			}
			printf ("done P\n");
		} else {
			if (binary) fseek(fp,sizeof(float)*size,SEEK_CUR);
			else for(i=0;i<size;i++) fscanf(fp,"%*f");
		}
	}
}


vtkreader::~vtkreader() {
}

