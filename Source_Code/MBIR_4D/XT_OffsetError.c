#include <stdio.h>
#include "nrutil.h" 
#include "XT_Constants.h"
#include "XT_Structures.h"
#include <mpi.h>
#include <math.h>
#include "XT_IOMisc.h"
#include "invert.h"
#include "allocate.h"

void gen_offset_constraint_windows (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr)
{
	int32_t r_size, t_size, num = 0, i, j, k, l, dim[4], N_t_node, N_r, N_t, node_rank, node_num, node_idx, k_idx, l_idx;
	char constraint_file[100] = "proj_constraint";

	node_rank = TomoInputsPtr->node_rank;
	node_num = TomoInputsPtr->node_num;
	N_r = SinogramPtr->N_r;
	N_t_node = SinogramPtr->N_t;
	N_t = N_t_node*node_num;

	r_size = 2*N_r/((int32_t)(sqrt(N_r) + 0.5));
	t_size = 2*N_t/((int32_t)(sqrt(N_t) + 0.5));
	
	for (i = 0; i <= N_r - r_size/2; i = i + r_size/2)	
	for (j = 0; j <= N_t - t_size/2; j = j + t_size/2) num++;	
	
	SinogramPtr->off_constraint = (Real_arr_t***)multialloc(sizeof(Real_arr_t), 3, num, N_r, N_t_node); 
	memset(&(SinogramPtr->off_constraint[0][0][0]), 0, num*N_r*N_t_node*sizeof(Real_arr_t));
	for (num = 0, i = 0; i <= N_r - r_size/2; i = i + r_size/2)
	for (j = 0; j <= N_t - t_size/2; j = j + t_size/2)	
	{	
		for (k = i; k < i + r_size; k++)
		for (l = j; l < j + t_size; l++)
		{
			node_idx = node_rank*N_t_node;
			k_idx = k % N_r;
			l_idx = l % N_t;
			if (l_idx >= node_idx && l_idx < node_idx + N_t_node)
			{
				SinogramPtr->off_constraint[num][k_idx][l_idx-node_idx] = (k-i) < r_size/2 ? (k-i+1): r_size-(k-i);
				SinogramPtr->off_constraint[num][k_idx][l_idx-node_idx] *= (l-j) < t_size/2 ? (l-j+1): t_size-(l-j);
			}
		}
		num++;
	}

	SinogramPtr->off_constraint_num = num;

	dim[0] = 1; dim[1] = num; dim[2] = SinogramPtr->N_r; dim[3] = SinogramPtr->N_t;
	sprintf(constraint_file, "%s_n%d", constraint_file, node_rank);
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (constraint_file, dim, 0, 1, 2, 3, &(SinogramPtr->off_constraint[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);

	fprintf(TomoInputsPtr->debug_file_ptr, "gen_offset_constraint_windows: r_size = %d, t_size = %d, number of constraints = %d\n", r_size, t_size, SinogramPtr->off_constraint_num);
/*	SinogramPtr->off_constraint_size = SinogramPtr->N_r;
	SinogramPtr->off_constraint = (Real_t**)multialloc(sizeof(Real_t), 2, 1, SinogramPtr->N_r); 
	for (j = 0; j < SinogramPtr->N_r; j++)
		SinogramPtr->off_constraint[0][j] = 1;
	SinogramPtr->off_constraint_num = 1;*/	
}


void constrained_quad_opt (Real_t** Lambda, Real_t** b, Real_arr_t*** A, Real_arr_t** x, int32_t Nr, int32_t Nt, int32_t M, TomoInputs* TomoInputsPtr)
{
	Real_t **D, **Dinv;
	Real_t *temp, *temp2;
	int32_t i, j, k, l;
	D = (Real_t**)multialloc(sizeof(Real_t), 2, M, M);
	Dinv = (Real_t**)multialloc(sizeof(Real_t), 2, M, M);
	temp = (Real_t*)get_spc(M, sizeof(Real_t));
	temp2 = (Real_t*)get_spc(M, sizeof(Real_t));
	memset(&(D[0][0]), 0, M*M*sizeof(Real_t));
	memset(&(Dinv[0][0]), 0, M*M*sizeof(Real_t));

    	#pragma omp parallel for collapse(2) private(k, l)
	for (i = 0; i < M; i++)	
	for (j = 0; j < M; j++)
	for (k = 0; k < Nr; k++)
	for (l = 0; l < Nt; l++)
	{
		D[i][j] += A[i][k][l]*A[j][k][l]/Lambda[k][l];
		/*sum += A[i][k]*A[j][k]/Lambda[k];*/
	}
	
/*  	TomoInputsPtr->t0_mpired1 = time(NULL); */
  	MPI_Allreduce(&(D[0][0]), &(Dinv[0][0]), M*M, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
/*  	TomoInputsPtr->time_mpired1 += difftime(time(NULL), TomoInputsPtr->t0_mpired1);*/
/*	printf("Checksum is %f\n", sum);*/
	invert2(Dinv, M);
	
    	#pragma omp parallel for private(j, k)
	for (i = 0; i < M; i++)
	{
		temp[i] = 0;
		for (j = 0; j < Nr; j++)
		for (k = 0; k < Nt; k++)
			temp[i] += A[i][j][k]*b[j][k]/Lambda[j][k];	
	}
 /* 	TomoInputsPtr->t0_mpired1 = time(NULL);*/ 
  	MPI_Allreduce(&(temp[0]), &(temp2[0]), M, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
/*  	TomoInputsPtr->time_mpired1 += difftime(time(NULL), TomoInputsPtr->t0_mpired1);*/
	
    	#pragma omp parallel for private(j)
	for (i = 0; i < M; i++)
	{
		temp[i] = 0;
		for (j = 0; j < M; j++)
			temp[i] += Dinv[i][j]*temp2[j];	
	}
/*  	TomoInputsPtr->t0_mpired1 = time(NULL); 
	MPI_Allreduce(&(temp[0]), &(temp2[0]), M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	TomoInputsPtr->time_mpired1 += difftime(time(NULL), TomoInputsPtr->t0_mpired1);*/

    	#pragma omp parallel for collapse(2) private(k)
	for (i = 0; i < Nr; i++)
	for (j = 0; j < Nt; j++)
	{
		x[i][j] = 0;
		for (k = 0; k < M; k++)
			x[i][j] += A[k][i][j]*temp[k];
	}	
	
    	#pragma omp parallel for collapse(2)
	for (i = 0; i < Nr; i++)
	for (j = 0; j < Nt; j++)
		x[i][j] = (b[i][j] - x[i][j])/Lambda[i][j];

	free(temp);
	free(temp2);
	multifree(D, 2);
	multifree(Dinv, 2);
}

void compute_d_constraint (Real_arr_t*** A, Real_arr_t **d, int32_t Nr, int32_t Nt, int32_t M, FILE* debug_file_ptr)
{
	int32_t i, j, k;
	Real_t *temp, *val;
 
	temp = (Real_t*)get_spc(M, sizeof(Real_t));
	val = (Real_t*)get_spc(M, sizeof(Real_t));
   	#pragma omp parallel for private(j, k)
	for (i = 0; i < M; i++)
	{
		temp[i] = 0;
		for (j = 0; j < Nr; j++)
		for (k = 0; k < Nt; k++)
			temp[i] += A[i][j][k]*d[j][k];	
	}
  	MPI_Allreduce(&(temp[0]), &(val[0]), M, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);

	for (i = 0; i < M; i++)
		fprintf(debug_file_ptr, "compute_d_constraint: The i th constraint on offset error is %f\n", val[i]);
	free(temp);
	free(val);
}
