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
int32_t read_data (char data_filename[], char whites_filename[], char darks_filename[], float *projections, float *weights, int32_t datafile_row0, int32_t proj_rows, int32_t proj_cols, int32_t proj_start, int32_t proj_num, FILE* debug_file_ptr)
{
	hid_t data_file_id, white_file_id, dark_file_id, white_dataset, dark_dataset, proj_dataset;
	hid_t white_dataspace, dark_dataspace, proj_dataspace, white_memspace, dark_memspace, proj_memspace;
	hsize_t dims_white[3], dims_dark[3], dims_proj[3], white_offset[3], dark_offset[3], proj_offset[3], white_count[3], dark_count[3], proj_count[3], mem_offset[3]; 
   	herr_t status;
	int32_t i, j, k, m, n, idx, white_rank, dark_rank, proj_rank, extras_r, true_length_r, ratio_r, ratio_t, nodes_num, nodes_rank, slice_num;	
	uint16_t ***white_img, ***dark_img, ***proj_img;
	float temp, **white_2D_img, **dark_2D_img, **wd_dwnsmpl_img;
	
	MPI_Comm_size(MPI_COMM_WORLD, &nodes_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &nodes_rank);
	
	/*HDF file pointers*/
	data_file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	white_file_id = H5Fopen(whites_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dark_file_id = H5Fopen(darks_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	/*dataset pointers*/
	white_dataset = H5Dopen(white_file_id, "/exchange/data_white", H5P_DEFAULT);
	dark_dataset = H5Dopen(dark_file_id, "/exchange/data_dark", H5P_DEFAULT);
	proj_dataset = H5Dopen(data_file_id, "/exchange/data", H5P_DEFAULT);

	white_dataspace = H5Dget_space (white_dataset);    /* dataspace handle */
	dark_dataspace = H5Dget_space (dark_dataset);    /* dataspace handle */
	proj_dataspace = H5Dget_space (proj_dataset);    /* dataspace handle */
	/*Gives the number of dimensions in a dataset*/
    	white_rank = H5Sget_simple_extent_ndims (white_dataspace);
    	dark_rank = H5Sget_simple_extent_ndims (dark_dataspace);
    	proj_rank = H5Sget_simple_extent_ndims (proj_dataspace);
	
	if (white_rank != 3 || dark_rank != 3 || proj_rank != 3)
	{
		fprintf(debug_file_ptr, "ERROR: The rank of one of the datasets in one of the HDF files is not 3.\n");
		return(-1);
	}

	/*finds the dimension of the dataset and stores them in dims_wd and dims_proj*/	
	status = H5Sget_simple_extent_dims (white_dataspace, dims_white, NULL);
	status = H5Sget_simple_extent_dims (dark_dataspace, dims_dark, NULL);
	status = H5Sget_simple_extent_dims (proj_dataspace, dims_proj, NULL);
	fprintf(debug_file_ptr, "Size of white (/exchange/data_white) dataset is %dx%dx%d\n", (int32_t)dims_white[0], (int32_t)dims_white[1], (int32_t)dims_white[2]);
	fprintf(debug_file_ptr, "Size of dark (/exchange/data_dark) dataset is %dx%dx%d\n", (int32_t)dims_dark[0], (int32_t)dims_dark[1], (int32_t)dims_dark[2]);
	fprintf(debug_file_ptr, "Size of count (/exchange/data) dataset is %dx%dx%d\n", (int32_t)dims_proj[0], (int32_t)dims_proj[1], (int32_t)dims_proj[2]);

	if (dims_white[2] != dims_proj[2] || dims_white[2] < proj_cols)
	{
		fprintf(debug_file_ptr, "ERROR: dims_white[2] = %d, dims_proj[2] = %d, proj_cols = %d\n", (int32_t)dims_white[2], (int32_t)dims_proj[2], proj_cols);
		return(-1);
	}
	
	if (dims_dark[2] != dims_proj[2] || dims_dark[2] < proj_cols)
	{
		fprintf(debug_file_ptr, "ERROR: dims_dark[2] = %d, dims_proj[2] = %d, proj_cols = %d\n", (int32_t)dims_dark[2], (int32_t)dims_proj[2], proj_cols);
		return(-1);
	}
	
        extras_r = dims_proj[2] % proj_cols;
        true_length_r = dims_proj[2] - extras_r;
	ratio_r = true_length_r/proj_cols;
	ratio_t = ratio_r; 
	
	proj_rows = proj_rows/nodes_num;
	slice_num = proj_rows*ratio_t;
	
	if (dims_white[1] != dims_proj[1] || dims_white[1] < datafile_row0 + slice_num*nodes_num)
	{
		fprintf(debug_file_ptr, "ERROR: dims_white[1] = %d and dims_proj[1] = %d should be at least %d\n", (int32_t)dims_white[1], (int32_t)dims_proj[1], datafile_row0 + slice_num*nodes_num);
		return(-1);
	}
	
	if (dims_dark[1] != dims_proj[1] || dims_dark[1] < datafile_row0 + slice_num*nodes_num)
	{
		fprintf(debug_file_ptr, "ERROR: dims_dark[1] = %d and dims_proj[1] = %d should be at least %d\n", (int32_t)dims_dark[1], (int32_t)dims_proj[1], datafile_row0 + slice_num*nodes_num);
		return(-1);
	}

	white_offset[0] = 1;
	white_offset[1] = datafile_row0 + nodes_rank*slice_num;
    	white_offset[2] = extras_r/2;
	
	dark_offset[0] = 1;
	dark_offset[1] = datafile_row0 + nodes_rank*slice_num;
    	dark_offset[2] = extras_r/2;
	
	proj_offset[0] = proj_start;
	proj_offset[1] = datafile_row0 + nodes_rank*slice_num;
    	proj_offset[2] = extras_r/2;

	white_count[0] = dims_white[0] - 2;
	white_count[1] = slice_num;
	white_count[2] = true_length_r;
	
	dark_count[0] = dims_dark[0] - 2;
	dark_count[1] = slice_num;
	dark_count[2] = true_length_r;

	proj_count[0] = proj_num;
	proj_count[1] = slice_num;
	proj_count[2] = true_length_r;

	white_img = (uint16_t***)multialloc(sizeof(uint16_t), 3, white_count[0], white_count[1], white_count[2]);
	dark_img = (uint16_t***)multialloc(sizeof(uint16_t), 3, dark_count[0], dark_count[1], dark_count[2]);
	proj_img = (uint16_t***)multialloc(sizeof(uint16_t), 3, proj_count[0], proj_count[1], proj_count[2]);
	white_2D_img = (float**)multialloc(sizeof(float), 2, white_count[1], white_count[2]);
	dark_2D_img = (float**)multialloc(sizeof(float), 2, dark_count[1], dark_count[2]);
	wd_dwnsmpl_img = (float**)multialloc(sizeof(float), 2, proj_rows, proj_cols);

	/*Selects ROI in the dataset which should be read into arrays*/
    	status = H5Sselect_hyperslab (white_dataspace, H5S_SELECT_SET, white_offset, NULL, white_count, NULL);
    	status = H5Sselect_hyperslab (dark_dataspace, H5S_SELECT_SET, dark_offset, NULL, dark_count, NULL);
    	status = H5Sselect_hyperslab (proj_dataspace, H5S_SELECT_SET, proj_offset, NULL, proj_count, NULL);
    	white_memspace = H5Screate_simple (3, white_count, NULL);  
    	dark_memspace = H5Screate_simple (3, dark_count, NULL);  
    	proj_memspace = H5Screate_simple (3, proj_count, NULL);   
	mem_offset[0] = 0; mem_offset[1] = 0; mem_offset[2] = 0;
    	status = H5Sselect_hyperslab (white_memspace, H5S_SELECT_SET, mem_offset, NULL, white_count, NULL);
    	status = H5Sselect_hyperslab (dark_memspace, H5S_SELECT_SET, mem_offset, NULL, dark_count, NULL);
    	status = H5Sselect_hyperslab (proj_memspace, H5S_SELECT_SET, mem_offset, NULL, proj_count, NULL);

	fprintf(debug_file_ptr,"Reading HDF5 dataset ...\n");
    	status = H5Dread(white_dataset, H5T_NATIVE_UINT16, white_memspace, white_dataspace, H5P_DEFAULT, &(white_img[0][0][0]));
    	status = H5Dread(dark_dataset, H5T_NATIVE_UINT16, dark_memspace, dark_dataspace, H5P_DEFAULT, &(dark_img[0][0][0]));
    	status = H5Dread(proj_dataset, H5T_NATIVE_UINT16, proj_memspace, proj_dataspace, H5P_DEFAULT, &(proj_img[0][0][0]));

	fprintf(debug_file_ptr, "Read the white, dark and the count data.\n");
	fprintf(debug_file_ptr, "Data sub-sampling factor along x-axis = %d, data sub-sampling factor along z-axis = %d\n", ratio_r, ratio_t);

	for (j = 0; j < white_count[1]; j++)
	for (k = 0; k < white_count[2]; k++)
	{
		white_2D_img[j][k] = 0;
		for (i = 0; i < white_count[0]; i++)
		{
			white_2D_img[j][k] += white_img[i][j][k];
		}
		white_2D_img[j][k] /= white_count[0];
	}

	for (j = 0; j < dark_count[1]; j++)
        for (k = 0; k < dark_count[2]; k++)
        {
                dark_2D_img[j][k] = 0;
                for (i = 0; i < dark_count[0]; i++)
                {
                        dark_2D_img[j][k] += dark_img[i][j][k];
                }
                dark_2D_img[j][k] /= dark_count[0];
        }

	fprintf(debug_file_ptr,"Generated the downsampled versions of white and dark images\n");
		
	for (i = 0; i < proj_num; i++)
	for (j = 0; j < proj_rows; j++)
	for (k = 0; k < proj_cols; k++)
	{
		idx = i*proj_rows*proj_cols + j*proj_cols + k;
		projections[idx] = 0;
		weights[idx] = 0;
		wd_dwnsmpl_img[j][k] = 0;
		for (m = 0; m < ratio_t; m++)
			for (n = 0; n < ratio_r; n++)
			{
				weights[idx] += fabs(proj_img[i][j*ratio_t + m][k*ratio_r + n] - dark_2D_img[j*ratio_t + m][k*ratio_r + n]);
				wd_dwnsmpl_img[j][k] += fabs(white_2D_img[j*ratio_t + m][k*ratio_r + n] - dark_2D_img[j*ratio_t + m][k*ratio_r + n]);
			}
		temp = log((wd_dwnsmpl_img[j][k])/(weights[idx]));
		projections[idx] = BH_QUAD_COEF*temp*temp + temp;
	}

	fprintf(debug_file_ptr,"Generated projections and weight data with beamhardening coefficient of %f\n", (float)BH_QUAD_COEF);
/*	
	if (TomoInputsPtr->Write2Tiff == 1)
	{
		dim[0] = 1; dim[1] = SinogramPtr->N_p; dim[2] = SinogramPtr->N_r; dim[3] = total_t_slices;
		WriteMultiDimArray2Tiff (proj_filename, dim, 0, 3, 1, 2, &(Projection[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		WriteMultiDimArray2Tiff (weight_filename, dim, 0, 3, 1, 2, &(Weight[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		dim[0] = 1; dim[1] = 1; dim[2] = SinogramPtr->N_r; dim[3] = total_t_slices;
		WriteMultiDimArray2Tiff (wd_filename, dim, 0, 1, 2, 3, &(wd_dwnsmpl_img[0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}
	fprintf(TomoInputsPtr->debug_file_ptr,"gen_projection_4m_HDF: Wrote projections, weight and bright field data to binary files\n");
	*/

	multifree(white_img, 3);
	multifree(dark_img, 3);
	multifree(proj_img, 3);
	multifree(wd_dwnsmpl_img, 2);
	multifree(white_2D_img, 2);
	multifree(dark_2D_img, 2);

	H5Sclose(white_memspace);
	H5Sclose(dark_memspace);
	H5Sclose(proj_memspace);
	H5Sclose(white_dataspace);
	H5Sclose(dark_dataspace);
	H5Sclose(proj_dataspace);
	H5Dclose(white_dataset);
	H5Dclose(dark_dataset);
	H5Dclose(proj_dataset);
	H5Fclose(white_file_id);
	H5Fclose(dark_file_id);
	H5Fclose(data_file_id);

	return (0);
} 



