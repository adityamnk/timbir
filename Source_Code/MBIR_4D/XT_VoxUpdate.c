

#include "XT_Constants.h"
#include <stdio.h>
#include "XT_Structures.h"
#include "XT_ICD_update.h"
#include "XT_AMatrix.h"
#include <math.h>
#include "allocate.h"

void compute_voxel_update_Atten (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino, AMatrixCol* AMatrixPtr, /*AMatrixCol* VoxelLineResponse,*/ Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1], int32_t i_new, int32_t slice, int32_t j_new, int32_t k_new)
{
  	int32_t p, q, r, sino_view, z_overlap_num;
	Real_t V,THETA1,THETA2,THETASelTemp;
	Real_t UpdatedVoxelValue, ProjectionEntry;
  	int32_t i_r, i_t;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; /*Store the present value of the voxel*/
	z_overlap_num = SinogramPtr->z_overlap_num;

	THETA1 = 0.0;
	THETA2 = 0.0;
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
	sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
	for (q = 0; q < AMatrixPtr[p].count; q++)
	{
      	    	i_r = (AMatrixPtr[p].index[q]);
       	    	ProjectionEntry = (AMatrixPtr[p].values[q]*SinogramPtr->delta_t);
/*       	ProjectionEntry = (AMatrixPtr[p].values[q]);
		for (r = 0; r < VoxelLineResponse[slice].count; r++)*/
		for (r = 0; r < z_overlap_num; r++)
		{ 
			/*i_t = VoxelLineResponse[slice].index[r];*/
			i_t = slice*z_overlap_num + r;
			if (SinogramPtr->ProjSelect[sino_view][i_r][i_t] == true)
			{
	           		/*THETA2 += (VoxelLineResponse[slice].values[r]*VoxelLineResponse[slice].values[r]*ProjectionEntry*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
               			THETA1 += -(ErrorSino[sino_view][i_r][i_t]*VoxelLineResponse[slice].values[r]*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);*/
	           		THETA2 += (ProjectionEntry*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
               			THETA1 += -(ErrorSino[sino_view][i_r][i_t]*ProjectionEntry*TomoInputsPtr->Weight[sino_view][i_r][i_t]);
            		}
			else
			{
				THETASelTemp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])/fabs(ErrorSino[sino_view][i_r][i_t]);
	            		/*THETA2 += (VoxelLineResponse[slice].values[r]*VoxelLineResponse[slice].values[r]*ProjectionEntry*ProjectionEntry*THETASelTemp);
            			THETA1 += -(ErrorSino[sino_view][i_r][i_t]*VoxelLineResponse[slice].values[r]*ProjectionEntry*THETASelTemp);*/
	            		THETA2 += (ProjectionEntry*ProjectionEntry*THETASelTemp);
            			THETA1 += -(ErrorSino[sino_view][i_r][i_t]*ProjectionEntry*THETASelTemp);
			}
            	}
	}
        }


            /*Solve the 1-D optimization problem
            TODO : What if theta1 = 0 ? Then this will give error*/


        UpdatedVoxelValue = CE_FunctionalSubstitution(V, THETA1, THETA2, ScannedObjectPtr, TomoInputsPtr, Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag);
              
        ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] = UpdatedVoxelValue;
	
	for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++){
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		for (q = 0; q < AMatrixPtr[p].count; q++)
        	{
               	    	i_r = (AMatrixPtr[p].index[q]);
        	    	ProjectionEntry = (AMatrixPtr[p].values[q]*SinogramPtr->delta_t);
        	    	/*ProjectionEntry = (AMatrixPtr[p].values[q]);
			for (r = 0; r < VoxelLineResponse[slice].count; r++)*/
			for (r = 0; r < z_overlap_num; r++)
			{ 
				/*i_t = VoxelLineResponse[slice].index[r];
	        		ErrorSino[sino_view][i_r][i_t] -= (ProjectionEntry*VoxelLineResponse[slice].values[r]*(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] - V));*/
				i_t = slice*z_overlap_num + r;
	        		ErrorSino[sino_view][i_r][i_t] -= (ProjectionEntry*(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] - V));
	   			if (fabs(ErrorSino[sino_view][i_r][i_t]*sqrt(TomoInputsPtr->Weight[sino_view][i_r][i_t])) < TomoInputsPtr->ErrorSinoThresh)
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = true;
				else
					SinogramPtr->ProjSelect[sino_view][i_r][i_t] = false;
	   		}
		}
	}
}

Real_t updateVoxels_Atten (int32_t time_begin, int32_t time_end, int32_t slice_begin, int32_t slice_end, int32_t xy_begin, int32_t xy_end, int32_t* x_rand_select, int32_t* y_rand_select, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino, Real_arr_t** DetectorResponse_XY, /*AMatrixCol* VoxelLineResponse,*/ int32_t Iter, long int *zero_count, Real_arr_t** MagUpdateMap, uint8_t** Mask)
{
  int32_t p,q,r,slice,i_new,j_new,k_new,idxr,idxq,idxp,index_xy;
  Real_t V;
  bool ZSFlag;
  int32_t sino_view;
  int32_t z_min, z_max;
  Real_t total_vox_mag = 0.0;

  z_min = 0;
  z_max = ScannedObjectPtr->N_z + 1;
  if (TomoInputsPtr->node_rank == 0)
	z_min = 1;
  if (TomoInputsPtr->node_rank == TomoInputsPtr->node_num - 1)
	z_max = ScannedObjectPtr->N_z;

    Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM]; 
    Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1]; 
    bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM];
    bool Time_BDFlag[NHOOD_TIME_MAXDIM-1];

  int32_t maxview = find_max(ScannedObjectPtr->ProjNum, ScannedObjectPtr->N_time);
   /*printf ("maxview = %d, size of AMatrixCol = %d\n",maxview,sizeof(AMatrixCol));*/
  AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(maxview, sizeof(AMatrixCol));
  uint8_t AvgNumXElements = (uint8_t)ceil(3*ScannedObjectPtr->delta_xy/SinogramPtr->delta_r);
  
  for (p = 0; p < maxview; p++)
  {
  	AMatrixPtr[p].values = (Real_t*)get_spc(AvgNumXElements,sizeof(Real_t));
  	AMatrixPtr[p].index  = (int32_t*)get_spc(AvgNumXElements,sizeof(int32_t));
  }  

   for (i_new = time_begin; i_new <= time_end; i_new++) 
   {
      for (index_xy = xy_begin; index_xy <= xy_end; index_xy++) 
      {
    /*    printf ("Entering index\n");*/ 
	k_new = x_rand_select[index_xy];
        j_new = y_rand_select[index_xy];
    	MagUpdateMap[j_new][k_new] = 0;  
           /*printf ("Entering mask\n"); */
	  for (p = 0; p < ScannedObjectPtr->ProjNum[i_new]; p++)
    	  {
		sino_view = ScannedObjectPtr->ProjIdxPtr[i_new][p];
		calcAMatrixColumnforAngle(SinogramPtr, ScannedObjectPtr, DetectorResponse_XY, &(AMatrixPtr[p]), j_new, k_new, sino_view);
    	  }
          for (slice = slice_begin; slice <= slice_end; slice++) {
        /*  	printf ("Entering slice\n");*/ 
            /*For a given (i,j,k) store its 26 point neighborhood*/           
	    if (Mask[j_new][k_new] == 1)
	    {   
	 	if (i_new - 1 >= 0){
			Time_Nhood[0] = ScannedObjectPtr->Object[i_new-1][slice+1][j_new][k_new];
			Time_BDFlag[0] = true;
		}
		else 
		{
			Time_Nhood[0] = 0.0;
			Time_BDFlag[0] = false;
		}

	    	if (i_new + 1 < ScannedObjectPtr->N_time){
			Time_Nhood[1] = ScannedObjectPtr->Object[i_new+1][slice+1][j_new][k_new];
			Time_BDFlag[1] = true;
		}
		else
		{
			Time_Nhood[1] = 0.0;
			Time_BDFlag[1] = false;
		}
	
	
	 for (p = 0; p < NHOOD_Z_MAXDIM; p++)
	 {
		idxp = slice + p;
		if (idxp >= z_min && idxp <= z_max)
		{
	 		for (q = 0; q < NHOOD_Y_MAXDIM; q++)
         		{
	 			idxq = j_new + q - 1;
                		if(idxq >= 0 && idxq < ScannedObjectPtr->N_y)
         			{
					for (r = 0; r < NHOOD_X_MAXDIM; r++)
					{
		    				idxr = k_new + r - 1;
                    				if(idxr >= 0 && idxr < ScannedObjectPtr->N_x){
	                				Spatial_Nhood[p][q][r] = ScannedObjectPtr->Object[i_new][idxp][idxq][idxr];
        	        				Spatial_BDFlag[p][q][r] = true;
                    				}
						else
						{
	                				Spatial_Nhood[p][q][r] = 0.0;
                    					Spatial_BDFlag[p][q][r] = false;
						}
					}
				}
		 		else
				{
         				for (r = 0; r < NHOOD_X_MAXDIM; r++){
	                			Spatial_Nhood[p][q][r] = 0.0;
                    				Spatial_BDFlag[p][q][r] = false;
					}
				}
                	}
		}
		else
        	{ 
			for (q = 0; q < NHOOD_Y_MAXDIM; q++){
				for (r = 0; r < NHOOD_X_MAXDIM; r++){
	              			Spatial_Nhood[p][q][r] = 0.0;
                   			Spatial_BDFlag[p][q][r] = false;
				}
			}
               }
	}

        Spatial_Nhood[(NHOOD_Y_MAXDIM-1)/2][(NHOOD_X_MAXDIM-1)/2][(NHOOD_Z_MAXDIM-1)/2] = 0.0;
        V = ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]; /*Store the present value of the voxel*/

#ifdef ZERO_SKIPPING
			  /*Zero Skipping Algorithm*/
			 ZSFlag = true;
			 if(V == 0.0 && Iter > 1) /*Iteration starts from 1. Iteration 0 corresponds to initial cost before ICD*/
			  {
					if (Time_Nhood[0] > 0.0 || Time_Nhood[1] > 0.0)
						ZSFlag = false;
			
					for(p = 0; p < NHOOD_Y_MAXDIM; p++)
						for(q = 0; q < NHOOD_X_MAXDIM; q++)
					  		for(r = 0; r < NHOOD_Z_MAXDIM; r++)
							  	if(Spatial_Nhood[p][q][r] > 0.0)
							  	{
									  ZSFlag = false;
								 	  break;
							  	}
			  }
			  else
			  {
				  ZSFlag = false;
			  }
#else
			  ZSFlag = false; /*do ICD on all voxels*/
#endif /*#ifdef ZERO_SKIPPING*/
	if(ZSFlag == false)
	{
		compute_voxel_update_Atten (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, AMatrixPtr, /*VoxelLineResponse,*/ Spatial_Nhood, Time_Nhood, Spatial_BDFlag, Time_BDFlag, i_new, slice, j_new, k_new);
	    	MagUpdateMap[j_new][k_new] += fabs(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new] - V);
	    	total_vox_mag += fabs(ScannedObjectPtr->Object[i_new][slice+1][j_new][k_new]);
 	}
		else
		    (*zero_count)++;
       }
       }
       }
}

    
     for (p=0; p<maxview; p++)
     {
     	free(AMatrixPtr[p].values);
     	free(AMatrixPtr[p].index);
     }
     free(AMatrixPtr);
      return (total_vox_mag);
}


Real_t updateVoxels (int32_t time_begin, int32_t time_end, int32_t slice_begin, int32_t slice_end, int32_t xy_begin, int32_t xy_end, int32_t* x_rand_select, int32_t* y_rand_select, Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino, Real_arr_t** DetectorResponse_XY, /*AMatrixCol* VoxelLineResponse,*/ int32_t Iter, long int *zero_count, Real_arr_t** MagUpdateMap, uint8_t** Mask)
{
	Real_t total_vox_mag;
	total_vox_mag = updateVoxels_Atten  (time_begin, time_end, slice_begin, slice_end, xy_begin, xy_end, x_rand_select, y_rand_select, SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse_XY, /*VoxelLineResponse,*/ Iter, zero_count, MagUpdateMap, Mask);
	return (total_vox_mag);
}
