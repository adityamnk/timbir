/* ============================================================================
 * Copyright (c) 2013 K. Aditya Mohan (Purdue University)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of K. Aditya Mohan, Purdue
 * University, nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "XT_Structures.h"
#include "XT_IOMisc.h"
#include <math.h>
/*Sends the x-y slices in from each z-block assigned to each node to the neighboring node's (in terms of rank) z-blocks.
This ensures that when each node does ICD on the assigned blocks they have information about the neighboring 
x-y slices which aren't updated in the same node.*/
/*send_reqs - contains information about send requests
recv_reqs - contains information on receive requests 
select - chooses whether to communicate the top x-y slice or the bottom in relation to odd and even time sllices*/
void MPI_Send_Recv_Z_Slices (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, MPI_Request* send_reqs, MPI_Request* recv_reqs, uint8_t select)
{
	int32_t i, num, N_z, off1, off2;
	N_z = ScannedObjectPtr->N_z;
	num = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;

	if (select == 0)
	{
		off1 = 0;
		off2 = 1;
	}
	else
	{
		off1 = 1;
		off2 = 0;
	}

	for (i = off1; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Isend(&(ScannedObjectPtr->Object[i][1][0][0]), num, MPI_REAL_ARR_DATATYPE, TomoInputsPtr->node_rank - 1, i, MPI_COMM_WORLD, &(send_reqs[i]));
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Irecv(&(ScannedObjectPtr->Object[i][N_z+1][0][0]), num, MPI_REAL_ARR_DATATYPE, TomoInputsPtr->node_rank + 1, i, MPI_COMM_WORLD, &(recv_reqs[i]));	
	}
		
	for (i = off2; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Irecv(&(ScannedObjectPtr->Object[i][0][0][0]), num, MPI_REAL_ARR_DATATYPE, TomoInputsPtr->node_rank - 1, i, MPI_COMM_WORLD, &(recv_reqs[i]));	
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Isend(&(ScannedObjectPtr->Object[i][N_z][0][0]), num, MPI_REAL_ARR_DATATYPE, TomoInputsPtr->node_rank + 1, i, MPI_COMM_WORLD, &(send_reqs[i]));
	}
}
			
/*Once the x-y slices are initiated for communication, the offset error and variance parameter are updated since they
don't require the neighboring slices from neighboring nodes. Then, the function below is used to wait for all communication to end.
send_reqs - contains information about send requests
recv_reqs - contains information on receive requests 
select - chooses whether to communicate the top x-y slice or the bottom in relation to odd and even time sllices*/
void MPI_Wait_Z_Slices (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, MPI_Request* send_reqs, MPI_Request* recv_reqs, uint8_t select)
{
	int32_t i, off1, off2;
	/*N_z = ScannedObjectPtr->N_z;
	num = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;*/

	if (select == 0)
	{
		off1 = 0;
		off2 = 1;
	}
	else
	{
		off1 = 1;
		off2 = 0;
	}

	for (i = off1; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Wait(&(send_reqs[i]), MPI_STATUSES_IGNORE);
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Wait(&(recv_reqs[i]), MPI_STATUSES_IGNORE);
	}
		
	for (i = off2; i < ScannedObjectPtr->N_time; i = i + 2)
	{
		if (TomoInputsPtr->node_rank > 0)
			MPI_Wait(&(recv_reqs[i]), MPI_STATUSES_IGNORE);
		
		if (TomoInputsPtr->node_rank < TomoInputsPtr->node_num - 1)	
			MPI_Wait(&(send_reqs[i]), MPI_STATUSES_IGNORE);
	}
}
			

/*void compute_RMSE_Converged_Object(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, int32_t Iter)
{
	int32_t i, j, k, l;
	Real_t RMSE = 0, RMSE_ALL = 0, avg_conv = 0, avg_obj = 0;
	char rmse_filename[100] = "rmse_converged";

	for (i = 0; i < ScannedObjectPtr->N_time; i++)
		for (j = 0; j < ScannedObjectPtr->N_z; j++)
			for (k = 0; k < ScannedObjectPtr->N_y; k++)
				for (l = 0; l < ScannedObjectPtr->N_x; l++)
				{
					RMSE += (ScannedObjectPtr->Object[i][j+1][k][l] - ScannedObjectPtr->Conv_Object[i][j][k][l])*(ScannedObjectPtr->Object[i][j+1][k][l] - ScannedObjectPtr->Conv_Object[i][j][k][l]);
					avg_conv += ScannedObjectPtr->Conv_Object[i][j][k][l];
					avg_obj += ScannedObjectPtr->Object[i][j+1][k][l];
				}
		
  	MPI_Reduce(&RMSE, &RMSE_ALL, 1, MPI_REAL_DATATYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  	if (TomoInputsPtr->node_rank == 0)
	{
		RMSE_ALL = sqrt(RMSE_ALL/(TomoInputsPtr->node_num*ScannedObjectPtr->N_time*ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x));
		avg_conv = avg_conv/(TomoInputsPtr->node_num*ScannedObjectPtr->N_time*ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x);
		avg_obj = avg_obj/(TomoInputsPtr->node_num*ScannedObjectPtr->N_time*ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x);
		fprintf(TomoInputsPtr->debug_file_ptr, "compute_RMSE_Converged_Object: The RMSE between the reconstruction and converged result is %fHU, average of converged object is %f, average of reconstruction is %f\n", convert_um2HU(RMSE_ALL), avg_conv, avg_obj);
		if (Iter == 0)
    			Write2Bin(rmse_filename, 1, 1, 1, 1, &RMSE_ALL, TomoInputsPtr->debug_file_ptr);	
		else
    			Append2Bin(rmse_filename, 1, 1, 1, 1, &RMSE_ALL, TomoInputsPtr->debug_file_ptr);	
	}
}
*/ 
