/* ============================================================================
 * Copyright (c) 2015 K. Aditya Mohan (Purdue University)
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
#include <math.h>
#include <stdlib.h>
#include "hdf5.h"
#include <mpi.h>
#include "allocate.h"
#include "XT_Constants.h"

/*Instead of reading projection and weight data from binary files, use the function below to read those directly from
HDF files.*/
int32_t read_ProjWeightData (char data_filename[], float *projections, float *weights, int32_t datafile_row0, int32_t proj_rows, int32_t proj_cols, int32_t proj_start, int32_t proj_num, FILE* debug_file_ptr)
{
	hid_t projs_file_id, weights_file_id, weights_dataset, projs_dataset, weights_dataspace, projs_dataspace, weights_memspace, projs_memspace;
	hsize_t weights_dims[3], projs_dims[3], data_offset[3], data_count[3], mem_offset[3]; 
   	herr_t status;
	int32_t i, j, k, m, n, weights_rank, projs_rank, extras_r, true_length_r, ratio_r, ratio_t, idx, nodes_rank, nodes_num, slice_num;	
	float ***projs_img, temp, ***weights_img;

	MPI_Comm_size(MPI_COMM_WORLD, &nodes_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &nodes_rank);
	/*HDF file pointers*/
	projs_file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	weights_file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*dataset pointers*/
	weights_dataset = H5Dopen(weights_file_id, STD_WEIGHTS_DATASET_NAME, H5P_DEFAULT);
	projs_dataset = H5Dopen(projs_file_id, STD_PROJECTIONS_DATASET_NAME, H5P_DEFAULT);

	weights_dataspace = H5Dget_space (weights_dataset);    /* dataspace handle */
	projs_dataspace = H5Dget_space (projs_dataset);    /* dataspace handle */
	/*Gives the number of dimensions in a dataset*/
    	weights_rank = H5Sget_simple_extent_ndims (weights_dataspace);
    	projs_rank = H5Sget_simple_extent_ndims (projs_dataspace);
	
	if (weights_rank != 3)
		fprintf(debug_file_ptr, "ERROR: The rank of the dataset %s is not 3\n", STD_WEIGHTS_DATASET_NAME);
	if (projs_rank != 3)
		fprintf(debug_file_ptr, "ERROR: The rank of the dataset %s is not 3\n", STD_PROJECTIONS_DATASET_NAME);

	/*finds the dimension of the dataset and stores them in dims_wd and dims_proj*/	
	status = H5Sget_simple_extent_dims (weights_dataspace, weights_dims, NULL);
	status = H5Sget_simple_extent_dims (projs_dataspace, projs_dims, NULL);
	fprintf(debug_file_ptr, "Size of weights (%s) dataset is %dx%dx%d\n", STD_WEIGHTS_DATASET_NAME, (int32_t)weights_dims[0], (int32_t)weights_dims[1], (int32_t)weights_dims[2]);
	fprintf(debug_file_ptr, "Size of projections (%s) dataset is %dx%dx%d\n", STD_PROJECTIONS_DATASET_NAME, (int32_t)projs_dims[0], (int32_t)projs_dims[1], (int32_t)projs_dims[2]);

	if (weights_dims[0] != projs_dims[0] || weights_dims[1] != projs_dims[1] || weights_dims[2] != projs_dims[2])
	{
		fprintf(debug_file_ptr, "ERROR: Dimensions of weights dataset and projection datasets don't match\n");
		return(-1);
	}
	
        extras_r = projs_dims[2] % proj_cols;
        true_length_r = projs_dims[2] - extras_r;
	ratio_r = true_length_r/proj_cols;
	ratio_t = ratio_r;

	proj_rows = proj_rows/nodes_num;
	slice_num = proj_rows*ratio_t;
 
	data_offset[0] = proj_start;
	data_offset[1] = datafile_row0 + nodes_rank*slice_num;
    	data_offset[2] = extras_r/2;

	data_count[0] = proj_num;
	data_count[1] = slice_num;
	data_count[2] = true_length_r;
	
	if (data_offset[0] + data_count[0] > projs_dims[0] || data_offset[1] + data_count[1] > projs_dims[1] || data_offset[2] + data_count[2] > projs_dims[2])
	{
		fprintf(debug_file_ptr, "ERROR: Dataset size is inconsistent with inputs\n");
		return(-1);
	} 

	fprintf(debug_file_ptr, "Data sub-sampling factor along x-axis = %d, data sub-sampling factor along z-axis = %d\n", ratio_r, ratio_t);
	weights_img = (float***)multialloc(sizeof(float), 3, data_count[0], data_count[1], data_count[2]);
	projs_img = (float***)multialloc(sizeof(float), 3, data_count[0], data_count[1], data_count[2]);

	if (weights_img == NULL || projs_img == NULL)
                fprintf(debug_file_ptr, "ERROR: calloc() returned NULL!");

	/*Selects ROI in the dataset which should be read into arrays*/
    	status = H5Sselect_hyperslab (weights_dataspace, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
    	status = H5Sselect_hyperslab (projs_dataspace, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
    	weights_memspace = H5Screate_simple (3, data_count, NULL);  
    	projs_memspace = H5Screate_simple (3, data_count, NULL);   
	mem_offset[0] = 0; mem_offset[1] = 0; mem_offset[2] = 0;
    	status = H5Sselect_hyperslab (weights_memspace, H5S_SELECT_SET, mem_offset, NULL, data_count, NULL);
    	status = H5Sselect_hyperslab (projs_memspace, H5S_SELECT_SET, mem_offset, NULL, data_count, NULL);

	fprintf(debug_file_ptr,"Reading HDF5 dataset ...\n");
    	status = H5Dread(weights_dataset, H5T_NATIVE_FLOAT, weights_memspace, weights_dataspace, H5P_DEFAULT, &(weights_img[0][0][0]));
	fprintf(debug_file_ptr,"Read dataset %s\n", STD_WEIGHTS_DATASET_NAME);
    	status = H5Dread(projs_dataset, H5T_NATIVE_FLOAT, projs_memspace, projs_dataspace, H5P_DEFAULT, &(projs_img[0][0][0]));
	fprintf(debug_file_ptr,"Read dataset %s\n", STD_PROJECTIONS_DATASET_NAME);

  	#pragma omp parallel for private(j, k, m, n, temp, idx)
	for (i = 0; i < proj_num; i++)
	{
		for (j = 0; j < proj_rows; j++)
		for (k = 0; k < proj_cols; k++)
		{
			temp = 0;
			idx = i*proj_rows*proj_cols + j*proj_cols + k;
			weights[idx] = 0;
			for (m = 0; m < ratio_t; m++)
			{
				for (n = 0; n < ratio_r; n++)
				{
					temp += projs_img[i][j*ratio_t + m][k*ratio_r + n];
					weights[idx] += weights_img[i][j*ratio_t + m][k*ratio_r + n];
				}
			}
			temp = temp/(ratio_r*ratio_t);
			projections[idx] = BH_QUAD_COEF*temp*temp + temp;
		}
	}

	fprintf(debug_file_ptr,"Generated projections and weight data with beamhardening coefficient of %f\n", (float)BH_QUAD_COEF);
	
/*	if (TomoInputsPtr->Write2Tiff == 1)
	{
		dim[0] = 1; dim[1] = SinogramPtr->N_p; dim[2] = SinogramPtr->N_r; dim[3] = total_t_slices;
		WriteMultiDimArray2Tiff (projs_filename, dim, 0, 3, 1, 2, &(Projection[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		WriteMultiDimArray2Tiff (weights_filename, dim, 0, 3, 1, 2, &(Weight[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}
	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Wrote projections and weight data to binary files\n");
*/	
	multifree(weights_img, 3);
	multifree(projs_img, 3);

	H5Sclose(weights_memspace);
	H5Sclose(projs_memspace);
	H5Sclose(weights_dataspace);
	H5Sclose(projs_dataspace);
	H5Dclose(weights_dataset);
	H5Dclose(projs_dataset);
	H5Fclose(weights_file_id);
	H5Fclose(projs_file_id);
	return(0);
} 

int32_t read_AngleTimeReconList (char data_filename[], float *proj_angles, float *proj_times, /*float *recon_times,*/ int32_t proj_start, int32_t proj_num, /*int32_t recon_num,*/ FILE *debug_file_ptr)
{
	hid_t angles_file_id, times_file_id, /*recon_file_id,*/ angles_dataset, times_dataset /*,recon_dataset*/;
	hid_t angles_dataspace, times_dataspace, /*recon_dataspace,*/ angles_memspace, times_memspace /*,recon_memspace*/;
	hsize_t angles_dims[1], times_dims[1], /*recon_dims[1],*/ data_offset[1], data_count[1], mem_offset[1]; 
   	herr_t status;
	
	int32_t angles_rank, times_rank/*, recon_rank*/;	

	angles_file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	times_file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*recon_file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);*/

	angles_dataset = H5Dopen(angles_file_id, STD_PROJ_ANGLES_DATASET_NAME, H5P_DEFAULT);
	times_dataset = H5Dopen(times_file_id, STD_PROJ_TIMES_DATASET_NAME, H5P_DEFAULT);
	/*recon_dataset = H5Dopen(recon_file_id, STD_RECON_TIMES_DATASET_NAME, H5P_DEFAULT);*/

	angles_dataspace = H5Dget_space (angles_dataset);    /* dataspace handle */
	times_dataspace = H5Dget_space (times_dataset);    /* dataspace handle */
	/*recon_dataspace = H5Dget_space (recon_dataset);*/    /* dataspace handle */
    	
	angles_rank = H5Sget_simple_extent_ndims (angles_dataspace);
    	times_rank = H5Sget_simple_extent_ndims (times_dataspace);
    	/*recon_rank = H5Sget_simple_extent_ndims (recon_dataspace);*/
	
	if (angles_rank != 1)
		fprintf(debug_file_ptr, "ERROR: The rank of the dataset %s is not 1\n", STD_PROJ_ANGLES_DATASET_NAME);
	if (times_rank != 1)
		fprintf(debug_file_ptr, "ERROR: The rank of the dataset %s is not 1\n", STD_PROJ_TIMES_DATASET_NAME);
	/*if (recon_rank != 1)
		fprintf(debug_file_ptr, "ERROR: The rank of the dataset %s is not 1\n", STD_RECON_TIMES_DATASET_NAME);*/
	
	status = H5Sget_simple_extent_dims (angles_dataspace, angles_dims, NULL);
	status = H5Sget_simple_extent_dims (times_dataspace, times_dims, NULL);
	/*status = H5Sget_simple_extent_dims (recon_dataspace, recon_dims, NULL);*/
	
	fprintf(debug_file_ptr, "Size of projection angles (%s) dataset is %d\n", STD_PROJ_ANGLES_DATASET_NAME, (int32_t)angles_dims[0]);
	fprintf(debug_file_ptr, "Size of projection times (%s) dataset is %d\n", STD_PROJ_TIMES_DATASET_NAME, (int32_t)times_dims[0]);
	/*fprintf(debug_file_ptr, "Size of reconstruction times (%s) dataset is %d\n", STD_RECON_TIMES_DATASET_NAME, (int32_t)recon_dims[0]);*/
	
	if (times_dims[0] != angles_dims[0])
	{
		fprintf(debug_file_ptr, "ERROR: Size of angles and times list does not match\n");
		return(-1);
	}

	data_offset[0] = proj_start;
	data_count[0] = proj_num;
	if (data_offset[0] + data_count[0] > angles_dims[0])
	{
		fprintf(debug_file_ptr, "ERROR: Dataset size is inconsistent with inputs\n");
		return(-1);
	} 
	
	status = H5Sselect_hyperslab (angles_dataspace, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
    	status = H5Sselect_hyperslab (times_dataspace, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
    	angles_memspace = H5Screate_simple (1, data_count, NULL);   
    	times_memspace = H5Screate_simple (1, data_count, NULL);  
	mem_offset[0] = 0; 
    	status = H5Sselect_hyperslab (angles_memspace, H5S_SELECT_SET, mem_offset, NULL, data_count, NULL);
    	status = H5Sselect_hyperslab (times_memspace, H5S_SELECT_SET, mem_offset, NULL, data_count, NULL);

    	status = H5Dread(angles_dataset, H5T_NATIVE_FLOAT, angles_memspace, angles_dataspace, H5P_DEFAULT, &(proj_angles[0]));
	fprintf(debug_file_ptr, "Read dataset %s\n", STD_PROJ_ANGLES_DATASET_NAME);
    	status = H5Dread(times_dataset, H5T_NATIVE_FLOAT, times_memspace, times_dataspace, H5P_DEFAULT, &(proj_times[0]));
	fprintf(debug_file_ptr, "Read dataset %s\n", STD_PROJ_TIMES_DATASET_NAME);
  /*	
	data_offset[0] = 0;
	data_count[0] = recon_num + 1;
	if (data_offset[0] + data_count[0] > recon_dims[0])
	{
		fprintf(debug_file_ptr, "ERROR: Dataset size is inconsistent with inputs\n");
		return(-1);
	} 
	status = H5Sselect_hyperslab (recon_dataspace, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
    	recon_memspace = H5Screate_simple (1, data_count, NULL);   
	mem_offset[0] = 0; 
    	status = H5Sselect_hyperslab (recon_memspace, H5S_SELECT_SET, mem_offset, NULL, data_count, NULL);
    	status = H5Dread(recon_dataset, H5T_NATIVE_FLOAT, recon_memspace, recon_dataspace, H5P_DEFAULT, &(recon_times[0]));
	fprintf(debug_file_ptr, "Read dataset %s\n", STD_RECON_TIMES_DATASET_NAME);
*/
	H5Sclose(angles_memspace);
	H5Sclose(times_memspace);
/*	H5Sclose(recon_memspace);*/
	H5Sclose(angles_dataspace);
	H5Sclose(times_dataspace);
/*	H5Sclose(recon_dataspace);*/
	H5Dclose(angles_dataset);
	H5Dclose(times_dataset);
/*	H5Dclose(recon_dataset);*/
	H5Fclose(angles_file_id);
	H5Fclose(times_file_id);
/*	H5Fclose(recon_file_id);*/
	return (0);
}



