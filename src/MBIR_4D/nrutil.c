#include <stdlib.h>
#include <stdio.h>
#include "XT_Constants.h"

char *mget_spc();
void exit();

void nrerror(char error_text[])
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}



Real_t *vector(int nl,int nh)
{
	Real_t *v;

	v=(Real_t *)mget_spc(nh-nl+1,sizeof(Real_t));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)mget_spc(nh-nl+1,sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

Real_t *dvector(int nl, int nh)
{
	Real_t *v;

	v=(Real_t *)mget_spc(nh-nl+1,sizeof(Real_t));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}



Real_t **matrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	Real_t **m;

	m=(Real_t **) mget_spc(nrh-nrl+1,sizeof(Real_t*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(Real_t *) mget_spc(nch-ncl+1,sizeof(Real_t));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

Real_t **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	Real_t **m;

	m=(Real_t **) mget_spc(nrh-nrl+1,sizeof(Real_t*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(Real_t *) mget_spc(nch-ncl+1,sizeof(Real_t));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,**m;

	m=(int **)mget_spc(nrh-nrl+1,sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)mget_spc(nch-ncl+1,sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



Real_t **submatrix(Real_t **a,int oldrl,int oldrh,int oldcl,
                  int oldch,int newrl,int newcl)
{
	int i,j;
	Real_t **m;

	m=(Real_t **) mget_spc(oldrh-oldrl+1,sizeof(Real_t*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(Real_t *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_ivector(int *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_dvector(Real_t *v,int nl,int nh)
{
	free((char*) (v+nl));
}



void free_matrix(Real_t **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(Real_t **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_submatrix(Real_t **b,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (b+nrl));
}



Real_t **convert_matrix(Real_t *a,int nrl,int nrh,int ncl,int nch)
{
	int i,j,nrow,ncol;
	Real_t **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (Real_t **) mget_spc(nrow,sizeof(Real_t*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix(Real_t **b,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (b+nrl));
}
