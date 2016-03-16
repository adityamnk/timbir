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


/*initializes structures used when generating projection data from phantom*/
void initPhantomStructures(Sinogram* Sino, ScannedObject* ScanObj, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	*Sino = *SinogramPtr;
	*ScanObj = *ScannedObjectPtr;

    	ScanObj->N_x = TomoInputsPtr->phantom_N_xy;
        ScanObj->N_y = TomoInputsPtr->phantom_N_xy;
    	ScanObj->N_z = SinogramPtr->slice_num/TomoInputsPtr->node_num;
    
    	ScanObj->delta_xy = Sino->Length_R/ScanObj->N_x;
    	ScanObj->delta_z = Sino->Length_T/ScanObj->N_z;
 
	Sino->OffsetR = (ScanObj->delta_xy/sqrt(2.0)+Sino->delta_r/2.0)/DETECTOR_RESPONSE_BINS;
	Sino->OffsetT = ((ScanObj->delta_z/2) + Sino->delta_t/2)/DETECTOR_RESPONSE_BINS;
}

/*generates projection data from phantom*/
void genSinogramFromPhantom (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
	FILE *fp;
	Sinogram Phantom_Sino;
	ScannedObject Phantom_ScanObj;
	AMatrixCol *VoxelLineResponse;
	long int stream_offset, size, result;
	int32_t i, j, k, m, n, idx, t, slice; 
        float ***object;
	Real_t proj_avg, weight_avg, pixel, val, **H_r, *H_t, **offset_sim;
  	uint8_t AvgNumXElements, AvgNumZElements, zinger_added = 0;
	char projection_file[100] = PROJECTION_FILENAME;
	char offgndtruth_file[100] = "proj_offset_truth";
	char weight_file[100] = WEIGHT_MATRIX_FILENAME;
	char phantom_file[1000];
	char detect_file[] = "detector_forwardproj";
	char offset_sim_file[100] = "proj_offset_sim";
	int dimTiff[4];

	sprintf(projection_file, "%s_n%d", projection_file, TomoInputsPtr->node_rank);
	sprintf(weight_file, "%s_n%d", weight_file, TomoInputsPtr->node_rank);
	if (TomoInputsPtr->phantom_N_xy % SinogramPtr->N_r != 0 || TomoInputsPtr->phantom_N_z < SinogramPtr->slice_begin + SinogramPtr->slice_num){
		printf("ERROR: genSinogramFromPhantom: N_r = %d does not divide phantom_N_xy = %d or number of slices = %d is more than phantom_N_z = %d\n", SinogramPtr->N_r, TomoInputsPtr->phantom_N_xy, SinogramPtr->slice_num, TomoInputsPtr->phantom_N_z);
		exit (1);
	}
	
	/*printf("\nsizes = %d, %d, %d, %d\n", sizeof(int), sizeof(int64_t), sizeof(long int), sizeof(long long));*/
	
	initPhantomStructures(&Phantom_Sino, &Phantom_ScanObj, SinogramPtr, ScannedObjectPtr, TomoInputsPtr);	

	AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(1,sizeof(AMatrixCol));
	/*AvgNumXElements over estimates the total number of entries in a single column of A matrix when indexed by both voxel and angle*/
  	AvgNumXElements = (uint8_t)ceil(3*Phantom_ScanObj.delta_xy/(Phantom_Sino.delta_r));
  	AMatrixPtr->values = (Real_t*)get_spc((int32_t)AvgNumXElements,sizeof(Real_t));
  	AMatrixPtr->index  = (int32_t*)get_spc((int32_t)AvgNumXElements,sizeof(int32_t));

	object = (float***)multialloc(sizeof(float), 3, Phantom_ScanObj.N_z, Phantom_ScanObj.N_y, Phantom_ScanObj.N_x);
	memset(&(Phantom_Sino.Projection[0][0][0]), 0, Phantom_Sino.N_p*Phantom_Sino.N_t*Phantom_Sino.N_r*sizeof(Real_t));	

	H_r = (Real_t **)multialloc(sizeof(Real_t), 2, Phantom_Sino.N_p, DETECTOR_RESPONSE_BINS+1);
	H_t = (Real_t *)get_spc(DETECTOR_RESPONSE_BINS + 1, sizeof(Real_t));
	DetectorResponseProfile (H_r, H_t, &Phantom_Sino, &Phantom_ScanObj, TomoInputsPtr);
	
  	AvgNumZElements = (uint8_t)((Phantom_ScanObj.delta_z/Phantom_Sino.delta_t) + 2);
	
	VoxelLineResponse = (AMatrixCol*)get_spc(Phantom_ScanObj.N_z,sizeof(AMatrixCol));
	for (t = 0; t < Phantom_ScanObj.N_z; t++){
    		VoxelLineResponse[t].values = (Real_t*)get_spc(AvgNumZElements, sizeof(Real_t));
    		VoxelLineResponse[t].index = (int32_t*)get_spc(AvgNumZElements, sizeof(int32_t));
	}
	storeVoxelLineResponse(H_t, VoxelLineResponse, &Phantom_ScanObj, &Phantom_Sino);

	sprintf(phantom_file, "%s", PATH_TO_PHANTOM);
	fp = fopen (phantom_file, "rb");
	
	if (fp==NULL) {fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: genSinogramFromPhantom: error in reading file %s\n",phantom_file); exit (1);}		
	size = (long int)Phantom_ScanObj.N_z*(long int)Phantom_ScanObj.N_y*(long int)Phantom_ScanObj.N_x;

	fprintf(TomoInputsPtr->debug_file_ptr, "DebugMsg: genSinogramFromPhantom: Forward projecting phantoms - \n");	
	for (i=0; i<Phantom_Sino.N_p; i++){
		stream_offset = (long int)i*(long int)TomoInputsPtr->phantom_N_z*(long int)Phantom_ScanObj.N_y*(long int)Phantom_ScanObj.N_x;
		stream_offset += (long int)Phantom_ScanObj.N_z*(long int)Phantom_ScanObj.N_y*(long int)Phantom_ScanObj.N_x*(long int)TomoInputsPtr->node_rank;
		stream_offset += (long int)SinogramPtr->slice_begin*(long int)Phantom_ScanObj.N_y*(long int)Phantom_ScanObj.N_x;
		result = fseek (fp, stream_offset*sizeof(float), SEEK_SET);

  		if (result != 0) 
		{fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Error in seeking file %s, i = %d, stream_offset = %ld\n",phantom_file,i,stream_offset);}

		result = fread (&(object[0][0][0]), sizeof(float), size, fp);
  		if (result != size) 
		{fprintf(TomoInputsPtr->debug_file_ptr, "ERROR: Reading file %s, Number of elements read does not match required, number of elements read=%ld, i=%d, stream_offset=%ld, size=%ld\n",phantom_file,result,i,stream_offset,size);}

		for (j=0; j<Phantom_ScanObj.N_y; j++)
		for (k=0; k<Phantom_ScanObj.N_x; k++){	
	   	    	calcAMatrixColumnforAngle(&Phantom_Sino, &Phantom_ScanObj, H_r, AMatrixPtr, j, k, i); 
                	for (slice=0; slice<Phantom_ScanObj.N_z; slice++){
#ifdef PHANTOM_IN_HU
			    	pixel = convert_HU2um((Real_t)object[slice][j][k]);
#else
			    	pixel = (Real_t)(object[slice][j][k]);
#endif
	     	          	for (m=0; m<AMatrixPtr->count; m++){
                            		idx=AMatrixPtr->index[m];
                            		val=AMatrixPtr->values[m];
                            		for (n=0; n<VoxelLineResponse[slice].count; n++)
                                    		Phantom_Sino.Projection[i][idx][VoxelLineResponse[slice].index[n]] += pixel*val*VoxelLineResponse[slice].values[n];
	     			}
				/*proj_avg += Phantom_Sino.Projection[i][idx][VoxelLineResponse[slice].index[n]];
				if (proj_avg != proj_avg)
				{
					printf("genSinogramFromPhantom: proj_avg = %f, j = %d, k = %d, slice = %d, m = %d, n = %d, pixel = %f, val = %f, line resp = %f, proj = %f\n",proj_avg,j,k,slice,m,n,pixel,val,VoxelLineResponse[slice].values[n],Phantom_Sino.Projection[i][idx][VoxelLineResponse[slice].index[n]]);
					exit(1);
				}*/
			}
	  	 }
	}

	dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_p; dimTiff[3] = DETECTOR_RESPONSE_BINS+1;
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (detect_file, dimTiff, 0, 1, 2, 3, &(H_r[0][0]), 0, TomoInputsPtr->debug_file_ptr);
       
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

	offset_sim = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
	Read4mBin(offset_sim_file, 1, 1, SinogramPtr->N_r, SinogramPtr->N_t, &(offset_sim[0][0]), TomoInputsPtr->debug_file_ptr);
	dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
	if (TomoInputsPtr->Write2Tiff == 1)
		WriteMultiDimArray2Tiff (offgndtruth_file, dimTiff, 0, 1, 2, 3, &(offset_sim[0][0]), 0, TomoInputsPtr->debug_file_ptr);

	proj_avg = 0; weight_avg = 0;
	if (TomoInputsPtr->No_Projection_Noise == 0)
		fprintf(TomoInputsPtr->debug_file_ptr, "The variance parameter, sigma^2, used to add noise is %f\n", VAR_PARAM_4_SIM_DATA);	
	fprintf(TomoInputsPtr->debug_file_ptr,"The expected count is %f\n", EXPECTED_COUNTS_FOR_PHANTOM_DATA);	
	for (i=0; i < SinogramPtr->N_p; i++)
	for (j=0; j < SinogramPtr->N_r; j++)
	for (slice=0; slice < SinogramPtr->N_t; slice++)
	{
		val = SinogramPtr->Projection[i][j][slice];
		val = EXPECTED_COUNTS_FOR_PHANTOM_DATA*exp(-val);
		if (val == 0)
		{
			printf("ERROR: val = %f, projection = %f, i = %d, j = %d, slice = %d\n",val,SinogramPtr->Projection[i][j][slice], i, j, slice);
		}
		/*TomoInputsPtr->Weight[i][j] = val + sqrt(val)*random2();*/
		if (TomoInputsPtr->No_Projection_Noise == 1)
			TomoInputsPtr->Weight[i][j][slice] = fabs(val);
		else
			TomoInputsPtr->Weight[i][j][slice] = fabs(val + sqrt(VAR_PARAM_4_SIM_DATA*val)*normal());
		/*Scaling noise S.D. this way, is equivalent to using a incident photon count of EXPECTED_COUNTS_FOR_PHANTOM_DATA/VAR_PARAM_4_SIM_DATA,
 		adding Poisson noise (mean = standard deviation) and then scaling the result by a multiplicative factor of VAR_PARAM_4_SIM_DATA*/
		
		SinogramPtr->Projection[i][j][slice] = log(EXPECTED_COUNTS_FOR_PHANTOM_DATA/TomoInputsPtr->Weight[i][j][slice]) - offset_sim[j][slice];
		if (random2() > 0.999)
		{
			zinger_added = 1;
			SinogramPtr->Projection[i][j][slice] = 0;
		}	
		weight_avg += TomoInputsPtr->Weight[i][j][slice];	
		proj_avg += SinogramPtr->Projection[i][j][slice];	
	}

	proj_avg /= size;
	weight_avg /= size;
	printf("genSinogramFromPhantom: The average of all projection data after weight computation with/without noise is %f\n", proj_avg);
	printf("genSinogramFromPhantom: The average of all weight data is %f\n", weight_avg);
	Write2Bin (projection_file, 1, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t, &(SinogramPtr->Projection[0][0][0]), TomoInputsPtr->debug_file_ptr);
	
	dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
	if (TomoInputsPtr->Write2Tiff == 1)
	{	
		WriteMultiDimArray2Tiff (projection_file, dimTiff, 0, 3, 1, 2, &(SinogramPtr->Projection[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
		WriteMultiDimArray2Tiff (weight_file, dimTiff, 0, 3, 1, 2, &(TomoInputsPtr->Weight[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
	}
	
	if (zinger_added == 1)
		fprintf(TomoInputsPtr->debug_file_ptr, "genSinogramFromPhantom: Zingers were added\n");
	Write2Bin (weight_file, 1, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t, &(TomoInputsPtr->Weight[0][0][0]), TomoInputsPtr->debug_file_ptr);
	multifree(offset_sim, 2); 
}


