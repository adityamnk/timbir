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


#include "randlib.h"
#include "XT_Structures.h"
#include "allocate.h"
#include "XT_NHICD.h"
#include "XT_IOMisc.h"

void randomize_list_in_place (int32_t* x_array, int32_t* y_array, int32_t num)    
{
        int32_t temp, i, j;    
    
        for (i = 0; i < num; i++)    
        {    
                j = ceil(random2()*(num - i));
                j = i + j - 1;
    
                temp = x_array[j];
                x_array[j] = x_array[i];
                x_array[i] = temp;    
    
                temp = y_array[j];
                y_array[j] = y_array[i];
                y_array[i] = temp;    
    
 /*               printf ("i = %d, j = %d\n", i, j);   */ 
        }   
}


Real_t RandomizedSelect(Real_arr_t* A, uint32_t p, uint32_t r, uint32_t i)
{
	uint32_t q, k;
	if (p == r)
    	{
		return A[p];
    	}
	
	q = RandomizedPartition(A, p, r);
	k = q - p + 1;
	
	if (i == k)
    	{
		return A[q];
    	}
	else if  (i < k)
    	{
		return RandomizedSelect(A, p, q - 1, i) ;
    	}
	else return RandomizedSelect(A, q + 1, r, i - k);
}

uint32_t Partition (Real_arr_t* A, uint32_t p, uint32_t r)
{
	Real_t x, temp;
	uint32_t i, j;

	x = A[r];
	i = p - 1;
	for (j = p; j < r; j++)
    	{
		if (A[j] <= x)
        	{	
			i++;
			temp = A[i];
			A[i] = A[j];
			A[j] = temp;
        	}
    	}
	
	temp = A[i + 1];
	A[i + 1] = A[r];
	A[r] = temp;
	return (i + 1);
}

uint32_t RandomizedPartition(Real_arr_t* A, uint32_t p, uint32_t r)
{
	Real_t temp;
	uint32_t j, rnum;
	
	rnum = floor(random2()*(r - p + 1));
	rnum = (rnum == r - p + 1)? r - p: rnum;

	/*if (rnum != 0)
		printf ("rnum is non-zero and equal to %d\n", rnum); */

	j = p + rnum;
	temp = A[r];
	A[r] = A[j];
	A[j] = temp;
	
	return Partition(A, p, r);
}

void ComputeVSC (ScannedObject* ScannedObjectPtr, Real_arr_t** magUpdateMap, Real_arr_t** filtMagUpdateMap, Real_arr_t** filtMagUpdateMap_copy)
{
	int32_t i, j, p, q;
    	Real_t filter_op = 0;
	Real_t HammingWindow[5][5] = {{0.0064, 0.0432, 0.0800, 0.0432, 0.0064},
				  {0.0432, 0.2916, 0.5400, 0.2916, 0.0432},
				  {0.0800, 0.5400, 1.0000, 0.5400, 0.0800},
				  {0.0432, 0.2916, 0.5400, 0.2916, 0.0432},
				  {0.0064, 0.0432, 0.0800, 0.0432, 0.0064}};    

           /*for (p = -2; p <= 2; p++)
           for (q = -2; q <= 2; q++)
		printf("p = %d, q = %d, val = %f\n", p, q, HammingWindow[p+2][q+2]);*/
    	for (i = 0; i < ScannedObjectPtr->N_y; i++)
    	{
		for (j = 0; j < ScannedObjectPtr->N_x; j++)
        	{
            		filter_op = 0;
            		for (p = -2; p <= 2; p++)
        /*    		for (p = 0; p <= 0; p++)*/
            		{
                		for (q = -2; q <= 2; q++)
          /*      		for (q = 0; q <= 0; q++)*/
                		{
                    			if(i + p >= 0 && i + p < ScannedObjectPtr->N_y && j + q >= 0 && j + q < ScannedObjectPtr->N_x)
                    			{
                				filter_op += HammingWindow[p + 2][q + 2] * magUpdateMap[i + p][j + q];
                    			}
                		}
            		}
        		filtMagUpdateMap[i][j] = filter_op;
			filtMagUpdateMap_copy[i][j] = filter_op;
		}
    	}
}

int32_t VSC_based_Voxel_Line_Select (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t**** MagUpdateMap)
{
	Real_t thresh;
	int32_t i, j, k, idx = 0, block;
	/*Real_t avg = 0;*/
/*	char filename[] = "updatemap";
	char magname[] = "magupmap";
	char filtmagname[] = "filtmagupmap";
	int dim[4];
	Real_t*** updatemap = (Real_t***)multialloc(sizeof(Real_t), 3, ScannedObjectPtr->N_time, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);*/
	
	Real_arr_t** filtMagUpdateMap = (Real_arr_t**)multialloc(sizeof(Real_arr_t), 2, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);
	Real_arr_t** filtMagUpdateMap_copy = (Real_arr_t**)multialloc(sizeof(Real_arr_t), 2, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);

	for (i = 0; i < ScannedObjectPtr->N_time; i++)
	for (block = 0; block < TomoInputsPtr->num_z_blocks; block++)
	{
		idx = 0;
		TomoInputsPtr->NHICDSelectNum[i][block] = floor(TomoInputsPtr->UpdateSelectNum[i][block]/ScannedObjectPtr->NHICD_Iterations);	
		/*avg = 0;*/
		ComputeVSC (ScannedObjectPtr, MagUpdateMap[i][block], filtMagUpdateMap, filtMagUpdateMap_copy);
		thresh = RandomizedSelect(&(filtMagUpdateMap[0][0]), 0, ScannedObjectPtr->N_x*ScannedObjectPtr->N_y - 1, ScannedObjectPtr->N_x*ScannedObjectPtr->N_y - TomoInputsPtr->NHICDSelectNum[i][block]);
		for (j = 0; j < ScannedObjectPtr->N_y; j++)
			for (k = 0; k < ScannedObjectPtr->N_x; k++)
			{
				/*updatemap[i][j][k]=0;*/
				if (filtMagUpdateMap_copy[j][k] > thresh)
				{
					TomoInputsPtr->x_NHICD_select[i][block][idx] = k;
					TomoInputsPtr->y_NHICD_select[i][block][idx] = j;
					/*updatemap[i][j][k] = 1.0;*/
					/*avg += filtMagUpdateMap[j][k];*/
					idx++;
				}
			}
		
/*		if (idx != TomoInputsPtr->NHICDSelectNum[i][block])
			fprintf(stderr, "WARNING: VSC_based_Voxel_Line_Select: Number of elements = %d above threshold does not match required = %d\n", idx, TomoInputsPtr->NHICDSelectNum[i][block]);*/
		TomoInputsPtr->NHICDSelectNum[i][block] = idx;
/*		randomize_list_in_place (TomoInputsPtr->x_rand_select[i], TomoInputsPtr->y_rand_select[i], idx);*/
/*		printf ("VSC_based_Voxel_Line_Select: The number of lines in time %d is %d, average is %f, threshold is %f\n", i, idx, avg/idx, thresh);*/	
	}	
	
/*	dim[0] = 1; dim[1] = ScannedObjectPtr->N_time; dim[2] = ScannedObjectPtr->N_y; dim[3] = ScannedObjectPtr->N_x; 
	WriteMultiDimArray2Tiff (filename, dim, 0, 1, 2, 3, &(updatemap[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	WriteMultiDimArray2Tiff (filtmagname, dim, 0, 1, 2, 3, &(filtMagUpdateMap_copy[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	WriteMultiDimArray2Tiff (magname, dim, 0, 1, 2, 3, &(MagUpdateMap[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	multifree(updatemap, 3);*/
	
	multifree(filtMagUpdateMap, 2);
	multifree(filtMagUpdateMap_copy, 2);
	return (0);
} 

