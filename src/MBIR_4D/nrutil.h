#include "XT_Constants.h"
void nrerror(char error_text[]);
Real_t *vector(int nl,int nh);
int *ivector(int nl,int nh);
Real_t *dvector(int nl, int nh);
Real_t **matrix(int nrl,int nrh,int ncl,int nch);
Real_t **dmatrix(int nrl,int nrh,int ncl,int nch);
int **imatrix(int nrl,int nrh,int ncl,int nch);
Real_t **submatrix(Real_t **a,int oldrl,int oldrh,int oldcl,
                  int oldch,int newrl,int newcl);
void free_vector(Real_t *v,int nl,int nh);
void free_ivector(int *v,int nl,int nh);
void free_dvector(Real_t *v,int nl,int nh);
void free_matrix(Real_t **m,int nrl,int nrh,int ncl,int nch);
void free_dmatrix(Real_t **m,int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
void free_submatrix(Real_t **b,int nrl,int nrh,int ncl,int nch);
Real_t **convert_matrix(Real_t *a,int nrl,int nrh,int ncl,int nch);
void free_convert_matrix(Real_t **b,int nrl,int nrh,int ncl,int nch);
