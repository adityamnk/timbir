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
#include <math.h>
#include "XT_Structures.h"
#include "XT_AMatrix.h"
#include "allocate.h"
#include "XT_ICD_update.h"
#include "randlib.h"
#include "XT_Init.h"
#include "XT_Constants.h"
#include <time.h>
#include "XT_IOMisc.h"
#include <ctype.h>
#include <mpi.h>
#include "XT_Debug.h"

/*
	- Function Name : reconstruct
	- Inputs (in order)
		- float** object : Address of the pointer to the reconstructed object.
			*object is a pointer to a 1D array in raster order of size recon_num x proj_rows x proj_cols x proj_cols where 'recon_num' is the number of time samples in the reconstruction, 'proj_rows' is the number of projection slices (along the axis of rotation), and 'proj_cols' is the number columns in the projection (or number of pixels along a row of detector bins)
		- float* projections : Pointer to the projection data. 
			'projections' is a pointer to a 1D array in raster order of size proj_num x proj_cols x proj_rows. A projection is typically computed as the logarithm of the ratio of light intensity incident on the object to the measured intensity.
		- float* weights : Pointer to the weight data.
			'weights' is a pointer to a 1D array in raster order of size proj_num x proj_cols x proj_rows. Every entry of 'weights' is used to appropriately weigh each term of the 'projections' array in the likelihood term of MBIR.
		- float* proj_angles : Pointer to the list of angles at which the projections are acquired.
		- float* proj_times : Pointer to the list of times at which the projections are acquired.
		- float* recon_times : Pointer to a array of reconstruction times. The reconstruction is assumed to be peicewise constant with fixed/varying step-sizes. Thus, the array 'recon_times' contains the times at which the steps occur. For example, to do a single 3D reconstruction, we use a array of two elements, first element being the time of the first projection and the second element being the time of the last projection.  
		- int32_t proj_rows : Number of rows in the projection data i.e., the number of projection slices (along the axis of rotation)
		- int32_t proj_cols : Number of columns in the projection data i.e., the number of pixels of the detector perpendicular to the axis of rotation.
		- int32_t proj_num : Total number of projections used for reconstruction
		- int32_t recon_num : Total number of reconstruction time steps. Note that the number of elements in the array pointed by 'recon_times' is 'recon_num + 1'.
		- float vox_wid : Side length of each cubic voxel. The unit of reconstruction is the inverse unit of 'vox_wid'.
		- float rot_center : Center of rotation. For example, if the center of rotation coincides with the center of the detector then rot_center = proj_rows/2. 
		- float sig_s : Spatial regularization parameter to be varied to achieve the optimum reconstruction quality. Reducing 'sig_s' will make the reconstruction smoother and increasing it will make it sharper but also noisier.
		- float sig_t : Temporal regularization parameter to be varied to achieve best quality. Reducing it will increase temporal smoothness which can improve quality. However, excessive smoothing along time might introduce artifacts.
		- float c_s : Parameter of the spatial qGGMRF prior. It should be chosen such that c_s < 0.01*D/sig_s where 'D' is a rough estimate for the maximum change in value of the reconstruction along an edge in space.
		- float c_t : Parameter of the temporal qGGMRF prior. It should be chosen such that c_t < 0.01*D/sig_t where 'D' is a rough estimate for the maximum change in value of the reconstruction along a temporal edge.
		- float convg_thresh : Convergence threshold expressed as a percentage (chosen in the range of 0 to 100).
		- float remove_rings : Legal values are '0', '1', '2' and '3'. If '0', then the algorithm does not do any ring correction. If '1', it does ring correction by estimating an offset error in the projection data. If it is equal to or greater than '1', it does ring correction. '1' does a uncontrained optimization. '2' enforces a zero mean constraint on the offset errors. '3' enforces a zero constraint on the weighted average of the offset errors over overlapping rectangular patches. 
		- int32_t quad_convex : Legal values are '0' and '1'. If '1', then the algorithm uses a convex quadratic forward model. This model does not account for the zinger measurements which causes streak artifacts in the reconstruction. If '0', then the algorithm uses a generalized Huber function which models the effect of zingers. This reduces streak artifacts in the reconstruction. Also, using '1' disables estimation of variance parameter 'sigma' and '0' enables it.
		- float huber_delta : The parameter \delta of the generalized Huber function which models the effect of zingers. Legal values are in the range 0 to 1.
		- float huber_T : The threshold parameter T of the generalized Huber function. All positive values are legal values.
		- uint8_t restart : Legal values are '0' and '1'. It is typically used if the reconstruction gets killed for any reason. If '1', the reconstruction starts from the previously run multi-resolution stage.
		- FILE *debug_msg_ptr : Pointer to the file to which the debug messages should be directed. Use 'stdout' if you do not want to direct messages to a file on disk.
	- Outputs (Return value) : '0' implies a safe return.  
*/
int reconstruct (float **object, float *projections, float *weights, float *proj_angles, float *proj_times, float *recon_times, int32_t proj_rows, int32_t proj_cols, int32_t proj_num, int32_t recon_num, float vox_wid, float rot_center, float sig_s, float sig_t, float c_s, float c_t, float convg_thresh, int32_t remove_rings, int32_t quad_convex, float huber_delta, float huber_T, uint8_t restart, FILE *debug_msg_ptr)
{
	time_t start;	
	int32_t flag, multres_num, i, mult_xy[MAX_MULTRES_NUM], mult_z[MAX_MULTRES_NUM], num_nodes, rank, last_multres;
	Sinogram *SinogramPtr = (Sinogram*)get_spc(1,sizeof(Sinogram));
	ScannedObject *ScannedObjectPtr = (ScannedObject*)get_spc(1,sizeof(ScannedObject));
	TomoInputs *TomoInputsPtr = (TomoInputs*)get_spc(1,sizeof(TomoInputs));

	MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	check_info(rank==0, debug_msg_ptr, "Reconstructing the data ....\n");

	start = time(NULL);
	srandom2(761521);
	
	check_error(proj_rows % num_nodes != 0, rank==0, debug_msg_ptr, "The total number of nodes requested should divide the number of rows in the projection data.\n");
	check_error(proj_rows < MIN_PROJECTION_ROWS, rank==0, debug_msg_ptr, "The minimum number of projection rows should be %d.\n", MIN_PROJECTION_ROWS);
	check_error(proj_cols % 2 != 0, rank==0, debug_msg_ptr, "The number of projection columns should be even.\n");
	check_error(proj_rows/num_nodes % 2 != 0 || proj_rows/num_nodes < MIN_ROWS_PER_NODE, rank==0, debug_msg_ptr, "The number of projection rows divided by the number of nodes should be an even number greater than or equal to %d.\n", MIN_ROWS_PER_NODE);

	check_mem(SinogramPtr,rank==0,debug_msg_ptr);
	check_mem(ScannedObjectPtr,rank==0,debug_msg_ptr);
	check_mem(TomoInputsPtr,rank==0,debug_msg_ptr);	
	TomoInputsPtr->debug_file_ptr = debug_msg_ptr;

	/***************Logic to auto initialize number of resolutions ****************/

	multres_num = (int32_t)(log(((float)proj_cols)/MIN_XY_RECON_RES)/log(2.0) + 1);
	if (multres_num > MAX_MULTRES_NUM)
		multres_num = MAX_MULTRES_NUM;
	if (multres_num < 2)
		multres_num = 2;

	mult_xy[0] = 1; 
	mult_z[0] = 1;
	for (i = 1; i < multres_num; i++)
	{
		if (proj_cols % (mult_xy[i-1]*2) == 0)
		{	
			mult_xy[i] = mult_xy[i-1]*2;
			if (proj_rows/num_nodes % (mult_z[i-1]*2) == 0 && proj_rows/num_nodes/(mult_z[i-1]*2) >= MIN_ROWS_PER_NODE)
				mult_z[i] = mult_z[i-1]*2;
			else
				mult_z[i] = mult_z[i-1];
		}
		else
		{
			multres_num = i;
			break;
		}
	}

	int32_t multres_xy[MAX_MULTRES_NUM], multres_z[MAX_MULTRES_NUM];
	for (i = 0; i < multres_num; i++)
	{
		multres_xy[i] = mult_xy[multres_num-1-i];
		multres_z[i] = mult_z[multres_num-1-i];
	}

	/********************End of auto multi-res init code******************/

	if (restart == 1) /*Used to restart from an older run */
	{
		if (Read4mBin (RUN_STATUS_FILENAME, 1, 1, 1, 1, sizeof(int32_t), &last_multres, TomoInputsPtr->debug_file_ptr)) {goto error;}
	}
	else 	
	{
		last_multres = 0;
	}

	check_info(rank==0, TomoInputsPtr->debug_file_ptr, "Number of multi-resolution stages is %d.\n", multres_num);


	for (i = last_multres; i < multres_num; i++) /*Main loop performing the optimization*/
	{
		check_info(rank==0, TomoInputsPtr->debug_file_ptr, "Running multi-resolution stage %d with x-y voxel scale = %d and z voxel scale = %d.\n", i, multres_xy[i], multres_z[i]);
		check_info(rank==0, TomoInputsPtr->debug_file_ptr, "Initializing object with flag %d \n",TomoInputsPtr->initICD);		 
		if (initStructures (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, i, multres_xy, multres_z, projections, weights, proj_angles, proj_times, recon_times, proj_rows, proj_cols, proj_num, recon_num, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, quad_convex, huber_delta, huber_T)) {goto error;}
#ifdef EXTRA_DEBUG_MESSAGES
		check_debug(rank==0, TomoInputsPtr->debug_file_ptr, "SinogramPtr numerical variable values are N_r = %d, N_t = %d, N_p = %d, total_t_slices = %d, delta_r = %f, delta_t = %f, R0 = %f, RMax = %f, T0 = %f, TMax = %f, Length_R = %f, Length_T = %f, OffsetR = %f, OffsetT = %f, z_overlap_num = %d, off_constraint_size = %d, off_constraint_num = %d\n", SinogramPtr->N_r, SinogramPtr->N_t, SinogramPtr->N_p, SinogramPtr->total_t_slices, SinogramPtr->delta_r, SinogramPtr->delta_t, SinogramPtr->R0, SinogramPtr->RMax, SinogramPtr->T0, SinogramPtr->TMax, SinogramPtr->Length_R, SinogramPtr->Length_T, SinogramPtr->OffsetR, SinogramPtr->OffsetT, SinogramPtr->z_overlap_num, SinogramPtr->off_constraint_size, SinogramPtr->off_constraint_num);	
		check_debug(rank==0, TomoInputsPtr->debug_file_ptr, "ScannedObjectPtr numerical variable values are Length_X = %f, Length_Y = %f, Length_Z = %f, N_x = %d, N_y = %d, N_z = %d, N_time = %d, x0 = %f, y0 = %f, z0 = %f, delta_xy = %f, delta_z = %f, mult_xy = %f, mult_z = %f, BeamWidth = %f, Sigma_S = %f, Sigma_t = %f, C_S = %f, C_T = %f, NHICD_Iterations = %d, delta_recon = %f.\n", ScannedObjectPtr->Length_X, ScannedObjectPtr->Length_Y, ScannedObjectPtr->Length_Z, ScannedObjectPtr->N_x, ScannedObjectPtr->N_y, ScannedObjectPtr->N_z, ScannedObjectPtr->N_time, ScannedObjectPtr->x0, ScannedObjectPtr->y0, ScannedObjectPtr->z0, ScannedObjectPtr->delta_xy, ScannedObjectPtr->delta_z, ScannedObjectPtr->mult_xy, ScannedObjectPtr->mult_z, ScannedObjectPtr->BeamWidth, ScannedObjectPtr->Sigma_S, ScannedObjectPtr->Sigma_T, ScannedObjectPtr->C_S, ScannedObjectPtr->C_T, ScannedObjectPtr->NHICD_Iterations, ScannedObjectPtr->delta_recon);
		check_debug(rank==0, TomoInputsPtr->debug_file_ptr, "TomoInputsPtr numerical variable values are NumIter = %d, StopThreshold = %f, RotCenter = %f, radius_obj = %f, Sigma_S_Q = %f, Sigma_T_Q = %f, Sigma_S_Q_P = %f, Sigma_T_Q_P = %f, var_est = %f, alpha = %f, cost_thresh = %f, initICD = %d, Write2Tiff = %d, updateProjOffset = %d, no_NHICD = %d, WritePerIter = %d, num_z_blocks = %d, prevnum_z_blocks = %d, ErrorSinoThresh = %f, ErrorSinoDelta = %f, node_num = %d, node_rank = %d, updateVar = %d, initMagUpMap = %d, ErrorSinoCost = %f, Forward_Cost = %f, Prior_Cost = %f, num_threads = %d\n", TomoInputsPtr->NumIter, TomoInputsPtr->StopThreshold, TomoInputsPtr->RotCenter, TomoInputsPtr->radius_obj, TomoInputsPtr->Sigma_S_Q, TomoInputsPtr->Sigma_T_Q, TomoInputsPtr->Sigma_S_Q_P, TomoInputsPtr->Sigma_T_Q_P, TomoInputsPtr->var_est, TomoInputsPtr->alpha, TomoInputsPtr->cost_thresh, TomoInputsPtr->initICD, TomoInputsPtr->Write2Tiff, TomoInputsPtr->updateProjOffset, TomoInputsPtr->no_NHICD, TomoInputsPtr->WritePerIter, TomoInputsPtr->num_z_blocks, TomoInputsPtr->prevnum_z_blocks, TomoInputsPtr->ErrorSinoThresh, TomoInputsPtr->ErrorSinoDelta, TomoInputsPtr->node_num, TomoInputsPtr->node_rank, TomoInputsPtr->updateVar, TomoInputsPtr->initMagUpMap, TomoInputsPtr->ErrorSino_Cost, TomoInputsPtr->Forward_Cost, TomoInputsPtr->Prior_Cost, TomoInputsPtr->num_threads);
#endif /*end of extra debug messages*/
		flag = ICD_BackProject(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
		check_info(rank == 0, TomoInputsPtr->debug_file_ptr, "Time elapsed is %f minutes.\n", difftime(time(NULL), start)/60.0);
		check_error(flag != 0, rank == 0, TomoInputsPtr->debug_file_ptr, "Reconstruction failed!\n");
		if (TomoInputsPtr->WritePerIter == 0)
			if (write_ObjectProjOff2TiffBinPerIter (SinogramPtr, ScannedObjectPtr, TomoInputsPtr)) {goto error;}
		if (Write2Bin (RUN_STATUS_FILENAME, 1, 1, 1, 1, sizeof(int32_t), &last_multres, TomoInputsPtr->debug_file_ptr)) {goto error;}
	
		*object = Arr4DToArr1D(ScannedObjectPtr->Object);	
		if (i < multres_num - 1) {free(*object);}
		freeMemory(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
		check_info(rank==0, TomoInputsPtr->debug_file_ptr, "Completed multi-resolution stage %d.\n", i);
	}
	
	free(SinogramPtr);
	free(ScannedObjectPtr);
	free(TomoInputsPtr);
	/*free(projections);
	free(weights);*/
	check_info(rank==0, TomoInputsPtr->debug_file_ptr, "Exiting MBIR 4D\n");
	
	return (0);

error:
	freeMemory(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	if (SinogramPtr) free(SinogramPtr);
	if (ScannedObjectPtr) free(ScannedObjectPtr);
	if (TomoInputsPtr) free(TomoInputsPtr);
/*	if (projections) free(projections);
	if (weights) free(weights);*/
	return (-1);
	
}

/*
int project (char obj_filename[], float **projections, float **weights, float *proj_angles, float *offset_sim, int32_t proj_rows, int32_t proj_cols, int32_t proj_num, float vox_wid, float rot_center, int32_t add_rings, int32_t add_streaks, FILE *debug_msg_ptr)
{
	time_t start;	
	int32_t flag, i, num_nodes, rank;
	Sinogram *SinogramPtr = (Sinogram*)get_spc(1,sizeof(Sinogram));
	ScannedObject *ScannedObjectPtr = (ScannedObject*)get_spc(1,sizeof(ScannedObject));
	TomoInputs *TomoInputsPtr = (TomoInputs*)get_spc(1,sizeof(TomoInputs));

	MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	check_info(rank==0, debug_msg_ptr, "Forward projecting the object ....\n");

	start = time(NULL);
	srandom2(761521);

	check_error(proj_rows % num_nodes != 0, rank==0, debug_msg_ptr, "The total number of nodes requested should divide the number of rows in the projection data.\n");
	check_error(proj_rows < MIN_PROJECTION_ROWS, rank==0, debug_msg_ptr, "The minimum number of projection rows should be %d.\n", MIN_PROJECTION_ROWS);
	check_error(proj_cols % 2 != 0, rank==0, debug_msg_ptr, "The number of projection columns should be even.\n");
	check_error(proj_rows/num_nodes % 2 != 0 || proj_rows/num_nodes < MIN_ROWS_PER_NODE, rank==0, debug_msg_ptr, "The number of projection rows divided by the number of nodes should be an even number greater than or equal to %d.\n", MIN_ROWS_PER_NODE);

	TomoInputsPtr->debug_file_ptr = debug_msg_ptr;
	if (initPhantomStructures (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, offset_sim, proj_angles, proj_rows, proj_cols, proj_num, vox_wid, rot_center)) {goto error;}

#ifdef EXTRA_DEBUG_MESSAGES
		check_debug(rank==0, TomoInputsPtr->debug_file_ptr, "SinogramPtr numerical variable values are N_r = %d, N_t = %d, N_p = %d, total_t_slices = %d, delta_r = %f, delta_t = %f, R0 = %f, RMax = %f, T0 = %f, TMax = %f, Length_R = %f, Length_T = %f, OffsetR = %f, OffsetT = %f, z_overlap_num = %d\n", SinogramPtr->N_r, SinogramPtr->N_t, SinogramPtr->N_p, SinogramPtr->total_t_slices, SinogramPtr->delta_r, SinogramPtr->delta_t, SinogramPtr->R0, SinogramPtr->RMax, SinogramPtr->T0, SinogramPtr->TMax, SinogramPtr->Length_R, SinogramPtr->Length_T, SinogramPtr->OffsetR, SinogramPtr->OffsetT, SinogramPtr->z_overlap_num);	
		check_debug(rank==0, TomoInputsPtr->debug_file_ptr, "ScannedObjectPtr numerical variable values are Length_X = %f, Length_Y = %f, Length_Z = %f, N_x = %d, N_y = %d, N_z = %d, N_time = %d, x0 = %f, y0 = %f, z0 = %f, delta_xy = %f, delta_z = %f, mult_xy = %f, mult_z = %f, BeamWidth = %f, Sigma_S = %f, Sigma_t = %f, C_S = %f, C_T = %f, NHICD_Iterations = %d, delta_recon = %f.\n", ScannedObjectPtr->Length_X, ScannedObjectPtr->Length_Y, ScannedObjectPtr->Length_Z, ScannedObjectPtr->N_x, ScannedObjectPtr->N_y, ScannedObjectPtr->N_z, ScannedObjectPtr->N_time, ScannedObjectPtr->x0, ScannedObjectPtr->y0, ScannedObjectPtr->z0, ScannedObjectPtr->delta_xy, ScannedObjectPtr->delta_z, ScannedObjectPtr->mult_xy, ScannedObjectPtr->mult_z, ScannedObjectPtr->BeamWidth, ScannedObjectPtr->Sigma_S, ScannedObjectPtr->Sigma_T, ScannedObjectPtr->C_S, ScannedObjectPtr->C_T, ScannedObjectPtr->NHICD_Iterations, ScannedObjectPtr->delta_recon);
		check_debug(rank==0, TomoInputsPtr->debug_file_ptr, "TomoInputsPtr numerical variable values are NumIter = %d, StopThreshold = %f, RotCenter = %f, radius_obj = %f, Sigma_S_Q = %f, Sigma_T_Q = %f, Sigma_S_Q_P = %f, Sigma_T_Q_P = %f, var_est = %f, alpha = %f, cost_thresh = %f, initICD = %d, Write2Tiff = %d, updateProjOffset = %d, EnforceZeroMeanOffset = %d, no_NHICD = %d, WritePerIter = %d, num_z_blocks = %d, prevnum_z_blocks = %d, ErrorSinoThresh = %f, ErrorSinoDelta = %f, node_num = %d, node_rank = %d, updateVar = %d, initMagUpMap = %d, ErrorSinoCost = %f, Forward_Cost = %f, Prior_Cost = %f, num_threads = %d\n", TomoInputsPtr->NumIter, TomoInputsPtr->StopThreshold, TomoInputsPtr->RotCenter, TomoInputsPtr->radius_obj, TomoInputsPtr->Sigma_S_Q, TomoInputsPtr->Sigma_T_Q, TomoInputsPtr->Sigma_S_Q_P, TomoInputsPtr->Sigma_T_Q_P, TomoInputsPtr->var_est, TomoInputsPtr->alpha, TomoInputsPtr->cost_thresh, TomoInputsPtr->initICD, TomoInputsPtr->Write2Tiff, TomoInputsPtr->updateProjOffset, TomoInputsPtr->EnforceZeroMeanOffset, TomoInputsPtr->no_NHICD, TomoInputsPtr->WritePerIter, TomoInputsPtr->num_z_blocks, TomoInputsPtr->prevnum_z_blocks, TomoInputsPtr->ErrorSinoThresh, TomoInputsPtr->ErrorSinoDelta, TomoInputsPtr->node_num, TomoInputsPtr->node_rank, TomoInputsPtr->updateVar, TomoInputsPtr->initMagUpMap, TomoInputsPtr->ErrorSino_Cost, TomoInputsPtr->Forward_Cost, TomoInputsPtr->Prior_Cost, TomoInputsPtr->num_threads);
#endif
	flag = ForwardProject(SinogramPtr, ScannedObjectPtr, TomoInputsPtr, obj_filename);
	check_info(rank == 0, TomoInputsPtr->debug_file_ptr, "Time elapsed is %f minutes.\n", difftime(time(NULL), start)/60.0);
	check_error(flag != 0, rank == 0, TomoInputsPtr->debug_file_ptr, "Forward projection failed!\n");
	if (TomoInputsPtr->WritePerIter == 0)
		if (write_ObjectProjOff2TiffBinPerIter (SinogramPtr, ScannedObjectPtr, TomoInputsPtr)) {goto error;}
	if (Write2Bin (RUN_STATUS_FILENAME, 1, 1, 1, 1, sizeof(int32_t), &last_multres, TomoInputsPtr->debug_file_ptr)) {goto error;}
	freeMemory(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	
	free(SinogramPtr);
	free(ScannedObjectPtr);
	free(TomoInputsPtr);
	return (0);

error:
	freeMemory(SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
	if (SinogramPtr)
		free(SinogramPtr);
	if (ScannedObjectPtr)
		free(ScannedObjectPtr);
	if (TomoInputsPtr)
		free(TomoInputsPtr);
	return (-1);
	
}
*/
