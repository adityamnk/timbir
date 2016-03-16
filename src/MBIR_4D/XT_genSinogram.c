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





/*#include <iostream>*/
/*#include "TiffUtilities.h"*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "XT_Structures.h"
#include "XT_Constants.h"
#include "allocate.h"
#include <math.h>
#include "XT_IOMisc.h"
#include "XT_AMatrix.h"
#include "XT_Profile.h"
#include "randlib.h"
#include "XT_Filter.h"
#include "XT_Init.h"

/*generates projection data from phantom*/
int32_t ForwardProject (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, float *projections, float *weights, float *offset_sim)
{
	FILE *fp;
	AMatrixCol *VoxelLineResponse;
	long int stream_offset, size, result;
	int32_t i, j, k, m, n, idx, t, slice; 
	Real_t proj_avg, weight_avg, pixel, val, **H_r, *H_t;
  	uint8_t AvgNumXElements, AvgNumZElements, zinger_added = 0;
	char projection_file[100] = PROJECTION_FILENAME;
	char offgndtruth_file[100] = "proj_offset_truth";
	char weight_file[100] = WEIGHT_MATRIX_FILENAME;
	char phantom_file[1000];
	char detect_file[] = "detector_forwardproj";
	char offset_sim_file[100] = "proj_offset_sim";
	int dimTiff[4];

	AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	/*AvgNumXElements over estimates the total number of entries in a single column of A matrix when indexed by both voxel and angle*/
  	AvgNumXElements = (uint8_t)ceil(3*ScannedObjectPtr->delta_xy/(SinogramPtr->delta_r));
  	AMatrixPtr->values = (Real_t*)get_spc((int32_t)AvgNumXElements,sizeof(Real_t));
  	AMatrixPtr->index  = (int32_t*)get_spc((int32_t)AvgNumXElements,sizeof(int32_t));

	object = (float***)multialloc(sizeof(float), 3, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);
	memset(&(projections[0]), 0, SinogramPtr->N_p*SinogramPtr->N_t*SinogramPtr->N_r*sizeof(Real_t));	

	H_r = (Real_t **)multialloc(sizeof(Real_t), 2, SinogramPtr->N_p, DETECTOR_RESPONSE_BINS+1);
	H_t = (Real_t *)get_spc(DETECTOR_RESPONSE_BINS + 1, sizeof(Real_t));
	DetectorResponseProfile (H_r, H_t, SinogramPtr, ScanndObjectPtr, TomoInputsPtr);
	
  	AvgNumZElements = (uint8_t)((ScannedObjectPtr->delta_z/SinogramPtr->delta_t) + 2);
	
	VoxelLineResponse = (AMatrixCol*)get_spc(ScannedObjectPtr->N_z,sizeof(AMatrixCol));
	for (t = 0; t < ScannedObjectPtr->N_z; t++){
    		VoxelLineResponse[t].values = (Real_t*)get_spc(AvgNumZElements, sizeof(Real_t));
    		VoxelLineResponse[t].index = (int32_t*)get_spc(AvgNumZElements, sizeof(int32_t));
	}
	storeVoxelLineResponse(H_t, VoxelLineResponse, ScannedObjectPtr, SinogramPtr);

	sprintf(phantom_file, "%s", PATH_TO_PHANTOM);
	fp = fopen (phantom_file, "rb");
	
	check_error(fp==NULL, TomoInputsPtr->debug_file_ptr, "Error in reading file %s\n", phantom_file);		
	size = (long int)ScannedObjectPtr->N_z*(long int)ScannedObjectPtr->N_y*(long int)ScannedObjectPtr->N_x;

	check_info(TomoInputsPtr->debug_file_ptr, "Forward projecting phantom ...\n");	
	for (i=0; i<SinogramPtr->N_p; i++){
		stream_offset = (long int)i*(long int)PHANTOM_Z_SIZE*(long int)ScannedObjectPtr->N_y*(long int)ScannedObjectPtr->N_x;
		stream_offset += (long int)ScannedObjectPtr->N_z*(long int)ScannedObjectPtr->N_y*(long int)ScannedObjectPtr->N_x*(long int)TomoInputsPtr->node_rank;
		/*stream_offset += (long int)SinogramPtr->slice_begin*(long int)ScannedObjectPtr->N_y*(long int)ScannedObjectPtr->N_x;*/
		result = fseek (fp, stream_offset*sizeof(float), SEEK_SET);

  		check_error(result != 0, TomoInputsPtr->debug_file_ptr, "ERROR: Error in seeking file %s, i = %d, stream_offset = %ld\n",phantom_file,i,stream_offset);

		result = fread (&(object[0][0][0]), sizeof(float), size, fp);
  		check_error(result != size, TomoInputsPtr->debug_file_ptr, "ERROR: Reading file %s, Number of elements read does not match required, number of elements read=%ld, i=%d, stream_offset=%ld, size=%ld\n",phantom_file,result,i,stream_offset,size);

		for (j=0; j<ScannedObjectPtr->N_y; j++)
		for (k=0; k<ScannedObjectPtr->N_x; k++){	
	   	    	calcAMatrixColumnforAngle(SinogramPtr, ScannedObjectPtr, H_r, AMatrixPtr, j, k, i); 
                	for (slice=0; slice<ScannedObjectPtr->N_z; slice++){
			    	pixel = (Real_t)(object[slice][j][k]);
	     	          	for (m=0; m<AMatrixPtr->count; m++){
                            		idx=AMatrixPtr->index[m];
                            		val=AMatrixPtr->values[m];
                            		for (n=0; n<VoxelLineResponse[slice].count; n++)
                                    		projections[i*SinogramPtr->N_t*SinogramPtr->N_r + VoxelLineResponse[slice].index[n]*SinogramPtr->N_r + idx] += pixel*val*VoxelLineResponse[slice].values[n];
	     			}
				/*proj_avg += SinogramPtr->Projection[i][idx][VoxelLineResponse[slice].index[n]];
				if (proj_avg != proj_avg)
				{
					printf("genSinogramFromPhantom: proj_avg = %f, j = %d, k = %d, slice = %d, m = %d, n = %d, pixel = %f, val = %f, line resp = %f, proj = %f\n",proj_avg,j,k,slice,m,n,pixel,val,VoxelLineResponse[slice].values[n],SinogramPtr->Projection[i][idx][VoxelLineResponse[slice].index[n]]);
					exit(1);
				}*/
			}
	  	 }
	}

	dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_p; dimTiff[3] = DETECTOR_RESPONSE_BINS+1;
	if (TomoInputsPtr->Write2Tiff == 1)
		if (WriteMultiDimArray2Tiff (detect_file, dimTiff, 0, 1, 2, 3, &(H_r[0][0]), 0, TomoInputsPtr->debug_file_ptr)) {goto error;}
       

/*	offset_sim = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
	Read4mBin(offset_sim_file, 1, 1, SinogramPtr->N_r, SinogramPtr->N_t, &(offset_sim[0][0]), TomoInputsPtr->debug_file_ptr);*/
	dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_t; dimTiff[3] = SinogramPtr->N_r;
	if (TomoInputsPtr->Write2Tiff == 1)
		if (WriteMultiDimArray2Tiff (offgndtruth_file, dimTiff, 0, 1, 2, 3, &(offset_sim[0]), 0, TomoInputsPtr->debug_file_ptr)) {goto error;}

	proj_avg = 0; weight_avg = 0;
	if (TomoInputsPtr->No_Projection_Noise == 0)
		check_info(TomoInputsPtr->debug_file_ptr, "The variance parameter, sigma^2, used to add noise is %f\n", TomoInputsPtr->sim_var_param);	
	check_info(TomoInputsPtr->debug_file_ptr, "The expected count is %f\n", TomoInputsPtr->sim_exp_counts);	
	for (i=0; i < SinogramPtr->N_p; i++)
	for (slice=0; slice < SinogramPtr->N_t; slice++)
	for (j=0; j < SinogramPtr->N_r; j++)
	{
		idx = i*SinogramPtr->N_t*SinogramPtr->N_r + slice*SinogramPtr->N_r + j;
		val = projections[idx];
		val = TomoInputsPtr->sim_exp_counts*exp(-val);
		check_error(val == 0, "ERROR: val = %f, projection = %f, i = %d, j = %d, slice = %d\n", val, projections[idx], i, j, slice);
		/*TomoInputsPtr->Weight[i][j] = val + sqrt(val)*random2();*/
		if (TomoInputsPtr->No_Projection_Noise == 1)
			weights[idx] = fabs(val);
		else
			weights[idx] = fabs(val + sqrt(TomoInputsPtr->sim_var_param*val)*normal());
		/*Scaling noise S.D. this way, is equivalent to using a incident photon count of EXPECTED_COUNTS_FOR_PHANTOM_DATA/VAR_PARAM_4_SIM_DATA,
 		adding Poisson noise (mean = standard deviation) and then scaling the result by a multiplicative factor of VAR_PARAM_4_SIM_DATA*/
		
		projections[idx] = log(TomoInputsPtr->sim_exp_counts/weights[idx]) - offset_sim[slice*SinogramPtr->N_r + j];
		if (random2() > 0.999)
		{
			zinger_added = 1;
			projections[idx] = 0;
		}	
		weight_avg += weights[idx];	
		proj_avg += projections[idx];	
	}

	proj_avg /= size;
	weight_avg /= size;
	printf("genSinogramFromPhantom: The average of all projection data after weight computation with/without noise is %f\n", proj_avg);
	printf("genSinogramFromPhantom: The average of all weight data is %f\n", weight_avg);
	
	check_info (zinger_added == 1, TomoInputsPtr->debug_file_ptr, "Simulated the occurence of zingers\n");
	
	free(AMatrixPtr->values);
	free(AMatrixPtr->index);
        free(VoxelLineResponse->values);
        free(VoxelLineResponse->index);
	multifree(H_r,2);
	free(H_t);
        free(AMatrixPtr);
        free(VoxelLineResponse);
	fclose(fp);	
	multifree(object,3);
	return (0);
error:
	free(AMatrixPtr->values);
	free(AMatrixPtr->index);
        free(VoxelLineResponse->values);
        free(VoxelLineResponse->index);
	multifree(H_r,2);
	free(H_t);
        free(AMatrixPtr);
        free(VoxelLineResponse);
	fclose(fp);	
	multifree(object,3);
	return (-1);
}


