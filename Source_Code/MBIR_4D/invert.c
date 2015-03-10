#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "nrutil.h"
#include "allocate.h"
#include "XT_Constants.h"

 void ludcmp(Real_t **a,int n,int *indx,Real_t *d);
 void lubksb(Real_t **a,int n,int *indx,Real_t b[]);

void exit();


/* For invert 1 and 2 */
/*
main()
{
	int i,j;
	Real_t aa[3][3],*a[3];
        void invert1(),invert2();

        a[0]=&(aa[0][0]);
        a[1]=&(aa[1][0]);
        a[2]=&(aa[2][0]);

	a[0][0]=1; a[0][1]=2; a[0][2]=3;
	a[1][0]=1; a[1][1]=2; a[1][2]=2;
	a[2][0]=2; a[2][1]=1; a[2][2]=1;

	invert2((Real_t **)a,3);
	for(i=0; i<1000; i++) invert1((Real_t **)a,3);
	printf("finished\n");

	for(i=0; i<3; i++) {
	  for(j=0; j<3; j++) printf("%lf ",a[i][j]);
          printf("\n");
        }

        return(0);
}
*/


/* For invert 3 and 4 */
/*
main()
{
	int i,j;
	Real_t a[3][3];

	a[0][0]=1; a[0][1]=2; a[0][2]=3;
	a[1][0]=1; a[1][1]=2; a[1][2]=2;
	a[2][0]=2; a[2][1]=1; a[2][2]=1;

	invert((Real_t *)a,3);
	for(i=0; i<1000; i++) invert((Real_t *)a,3);
	printf("finished\n");

	for(i=0; i<3; i++) {
	  for(j=0; j<3; j++) printf("%lf ",a[i][j]);
          printf("\n");
        }

        return(0);
}
*/

/* answer = (
(0.0, -0.33333333333333337, 0.66666666666666663)
(-1.0, 1.66666666666666674, -0.333333333333333315)
(1.0, -1.0, 0.0)
)*/





/* inverts a matrix input as a 2D array. Requires n<=6 */
void invert1(
  Real_t **a, /* input/output matrix */
  int	n)    /* dimension */
{
	int  i,j,indx[7];
	Real_t  y[7][7],col[7],d;

	if(n>6) { printf("invert: n too large\n"); exit(-1); }

        /* Shift to index to a[1-n][1-n] */
	a--;
	for(i=1; i<=n; i++) a[i]--;

        /* Invert matrix */
	ludcmp(a,n,indx,&d);
	for(j=1; j<=n; j++) {
	  for(i=1; i<=n; i++) col[i]=0.0;
	  col[j]=1.0;
	  lubksb(a,n,indx,col);
	  for(i=1; i<=n; i++) y[i][j]=col[i];
	}

	for(i=1; i<=n; i++)
	for(j=1; j<=n; j++) a[i][j]=y[i][j];

        /* Shift to index back */
	for(i=1; i<=n; i++) a[i]++;
	a++;
}

/* inverts a matrix of arbitrary size input as a 2D array. */ 
void invert2(
  Real_t **a, /* input/output matrix */
  int	n)    /* dimension */
{
	int  i,j,*indx;
	Real_t  **y,*col,d,**dmatrix();

        /* Shift to index to a[1-n][1-n] */
	a--;
	for(i=1; i<=n; i++) a[i]--;

	indx = ivector(1,n);
	y = dmatrix(1,n,1,n); 
	col = (Real_t *)mget_spc(n,sizeof(Real_t)); col--;

	ludcmp(a,n,indx,&d);
	for(j=1; j<=n; j++) {
	  for(i=1; i<=n; i++) col[i]=0.0;
	  col[j]=1.0;
	  lubksb(a,n,indx,col);
	  for(i=1; i<=n; i++) y[i][j]=col[i];
	} 

	for(i=1; i<=n; i++)
	for(j=1; j<=n; j++) a[i][j]=y[i][j];

        /* Shift to index as a[1-n][1-n] */
	for(i=1; i<=n; i++) a[i]++;
	a++;

	free_ivector(indx,1,n);
	free_dmatrix(y,1,n,1,n);
	free((char *)(col+1));
}

/* inverts a matrix input as a 1D vector. Requires n<=6 */
void invert3(Real_t *a,int n)
{
	int  i,j,indx[7];
	Real_t  **aa,yy[7][7],col[7],d;

	if(n>6) { printf("invert: n too large\n"); exit(-1); }

	a--;
	aa = (Real_t **)mget_spc(n,sizeof(Real_t *)); aa--;
	for(i=1; i<=n; i++) aa[i]= a+(i-1)*n;

	ludcmp(aa,n,indx,&d);
	for(j=1; j<=n; j++) {
	  for(i=1; i<=n; i++) col[i]=0.0;
	  col[j]=1.0;
	  lubksb(aa,n,indx,col);
	  for(i=1; i<=n; i++) yy[i][j]=col[i];
	}

	for(i=1; i<=n; i++)
	for(j=1; j<=n; j++) aa[i][j]=yy[i][j];

	free((char *)(aa+1));
}


/* inverts a matrix of arbitrary size input as a 1D vector. */ 
void invert4(Real_t *a,int n)
{
	int  i,j,*indx;
	Real_t  **aa,**yy,*col,d,**dmatrix();

	a--;
	aa = (Real_t **)mget_spc(n,sizeof(Real_t *)); aa--;
	for(i=1; i<=n; i++) aa[i]= a+(i-1)*n;

	indx = ivector(1,n);
	yy = dmatrix(1,n,1,n);
	col = (Real_t *)mget_spc(n,sizeof(Real_t)); col--;

	ludcmp(aa,n,indx,&d);
	for(j=1; j<=n; j++) {
	  for(i=1; i<=n; i++) col[i]=0.0;
	  col[j]=1.0;
	  lubksb(aa,n,indx,col);
	  for(i=1; i<=n; i++) yy[i][j]=col[i];
	}

	for(i=1; i<=n; i++)
	for(j=1; j<=n; j++) aa[i][j]=yy[i][j];

	free_ivector(indx,1,n);
	free_dmatrix(yy,1,n,1,n);
	free((char *)(aa+1));
	free((char *)(col+1));
}


#define TINY 1.0e-20;

 void ludcmp(Real_t **a,int n,int *indx,Real_t *d)
{
	int i,imax,j,k;
	Real_t big,dum,sum,temp;
	Real_t *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}

#undef TINY

 void lubksb(Real_t **a,int n,int *indx,Real_t b[])
{
	int i,ii=0,ip,j;
	Real_t sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

