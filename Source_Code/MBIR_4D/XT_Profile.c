
/* ============================================================================
 * Copyright (c) 2013 K. Aditya Mohan (Purdue University)
 * Copyright (c) 2013 Singanallur Venkatakrishnan (Purdue University)
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
 * Neither the name of K. Aditya Mohan, Singanallur Venkatakrishnan, Purdue
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
#include "XT_Profile.h"
#include "allocate.h"
#include "XT_Constants.h"
#include <math.h>
#include "XT_IOMisc.h"
#include "XT_Debug.h"

/* 'calculateVoxelProfile' computes the voxel profile as a function of angle and detector index. Here we assume that we have PROFILE_RESOULUTION number of detector bins. Note that these bins are virtual and is not the actual resolution of the detector. All distance computations are normalized.
'VoxProfile' will contain the voxel profile*/
void calculateVoxelProfile(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t** VoxProfile)
{
	Real_t angle,MaxValLineIntegral=0;
	Real_t temp,dist1,dist2,LeftCorner,LeftNear,RightNear,RightCorner,t;

/*'MaxValLineIntegral' gives the maximum value of the line integral of a single voxel at a particular angle 'theta'. Also, note that in between the distances of 'LeftNear' and 'RightNear' the line integral will equal 'MaxValLineIntegral'. Beyond these points, to the left till 'LeftCorner' and to the right till 'RightCorner' the line integral is given by linear interpolation. To the left beyond 'LeftCorner' and to the right beyond 'RightCorner' the line integral is zero*/

	Real_t checksum=0;
	int32_t i, j;

	for (i=0;i<(int32_t)SinogramPtr->N_p;i++)
	{
		angle = SinogramPtr->ViewPtr[i];
		while(angle > M_PI_2)
			angle -= M_PI_2;

		while(angle < 0)
			angle +=M_PI_2;

		if(angle <= M_PI_4)
		{
			/*MaxValLineIntegral = ScannedObjectPtr->delta_xy/cos(angle)/SinogramPtr->delta_t;*/
			MaxValLineIntegral = ScannedObjectPtr->delta_xy/cos(angle);
		}
		else
		{
			/*MaxValLineIntegral = ScannedObjectPtr->delta_xy/cos(M_PI_2-angle)/SinogramPtr->delta_t;*/
			MaxValLineIntegral = ScannedObjectPtr->delta_xy/cos(M_PI_2-angle);
		}
/*MaxValLineIntegral gives the line integral through a voxel at a given angle*/
/*Divide by delta_t because it will be appropriately scaled by VoxelLineResponse*/
		temp=cos(M_PI_4);
		dist1 = temp * cos((M_PI_4 - angle));
		dist2 = temp * fabs((cos((M_PI_4 + angle))));
		LeftCorner = 1-dist1;
		LeftNear = 1-dist2;
		RightNear = 1+dist2;
		RightCorner = 1+dist1;

		for(j = 0;j<PROFILE_RESOLUTION;j++)
		{
			t = 2.0*j / PROFILE_RESOLUTION; /* 2 is the normalized length of the profile (basically equal to 2*delta_xy)*/
			if(t <= LeftCorner || t >= RightCorner)
				VoxProfile[i][j] = 0;
			else if(t > RightNear)
				VoxProfile[i][j] = MaxValLineIntegral*(RightCorner-t)/(RightCorner-RightNear);
			else if(t >= LeftNear)
				VoxProfile[i][j] = MaxValLineIntegral;
			else
				VoxProfile[i][j] = MaxValLineIntegral*(t-LeftCorner)/(LeftNear-LeftCorner);

			checksum+=VoxProfile[i][j];
		}
	}

	check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "MaxValLineIntegral = %f\n", MaxValLineIntegral);
	check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Calculated Voxel Profile with Check sum = %f\n",checksum);
}

/*The response of a single detector bin is not uniform all over its area. In general, the response is maximum at its center and goes down as we move away from the center. In this code, we assume that the response follows the Hamming Window type response.
BeamProfile contains the beam profile*/
void initializeBeamProfile(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t *BeamProfile)
{
  int32_t i;
 Real_t sum=0;
  for (i=0; i < BEAM_RESOLUTION ;i++)
  {
    BeamProfile[i] = 0.54 - 0.46*cos((2*M_PI/BEAM_RESOLUTION)*i);
    sum=sum+BeamProfile[i];
  }

  /*Normalize the beam to have an area of 1*/
  for (i=0; i < BEAM_RESOLUTION ;i++)
  {
    BeamProfile[i]/=sum;
  }
  check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Generated Beam Profile.\n");
}

/*'DetectorResponseProfile' computes the detector response as a function of the distance between the center of the voxel and the center of the detector bin. This response is computed for all the angles of rotation. Note that the response is the same irrespective of whether the voxel is to the right or the left of the center of the detector bin as long as it is at the same distance. Also, the distance is quantized to DETECTOR_RESPONSE_BINS number of bins.
H_r - Detector response along r - axis
H_t - Detector response along t - axis*/
void DetectorResponseProfile (Real_arr_t** H_r, Real_arr_t* H_t, Sinogram* SinogramPtr, ScannedObject *ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  Real_t r,sum=0,rmin,ProfileCenterR,ProfileCenterT,TempConst;
  Real_t r0 = -(ScannedObjectPtr->BeamWidth)/2;
  Real_t StepSize = (ScannedObjectPtr->BeamWidth)/BEAM_RESOLUTION;
  int32_t i,k,p,ProfileIndex;
  Real_arr_t** VoxProfile;
  Real_t* BeamProfile;
  char filename[100]="VoxelProfile";
  int dimTiff[4];  

  VoxProfile = (Real_arr_t**)multialloc(sizeof(Real_arr_t),2,SinogramPtr->N_p,PROFILE_RESOLUTION);
  calculateVoxelProfile(SinogramPtr, ScannedObjectPtr, TomoInputsPtr, VoxProfile);
  
  dimTiff[0] = 1;
  dimTiff[1] = 1;
  dimTiff[2] = SinogramPtr->N_p;
  dimTiff[3] = PROFILE_RESOLUTION;
  sprintf(filename, "%s_n%d", filename, TomoInputsPtr->node_rank);
  if (TomoInputsPtr->Write2Tiff == 1)
  WriteMultiDimArray2Tiff (filename, dimTiff, 0, 1, 2, 3, &(VoxProfile[0][0]), 0, TomoInputsPtr->debug_file_ptr);

  BeamProfile=(Real_t*)get_spc(BEAM_RESOLUTION,sizeof(Real_t));
  initializeBeamProfile(ScannedObjectPtr, TomoInputsPtr, BeamProfile);

  TempConst=(PROFILE_RESOLUTION)/(2*ScannedObjectPtr->delta_xy);

  for(k = 0 ; k < (int32_t) SinogramPtr->N_p; k++)
  {
    for (i = 0; i < DETECTOR_RESPONSE_BINS; i++) 
    {
      ProfileCenterR = i*SinogramPtr->OffsetR;
      rmin = ProfileCenterR - ScannedObjectPtr->delta_xy;
        
	sum = 0;
        for (p=0; p < BEAM_RESOLUTION; p++)
        {
          r = r0 + p*StepSize;

          ProfileIndex = (int32_t)floor((r - rmin) * TempConst);
          if(ProfileIndex >= 0 && ProfileIndex <= PROFILE_RESOLUTION-1)
         	 sum += (VoxProfile[k][ProfileIndex] * BeamProfile[p]);
        }
        H_r[k][i] = sum;
      }
      H_r[k][i] = 0.0;
    }
	free(BeamProfile);
	multifree(VoxProfile,2);

    for (i = 0; i < DETECTOR_RESPONSE_BINS; i++)
    {
      ProfileCenterT = i * SinogramPtr->OffsetT;
      if(ScannedObjectPtr->delta_z >= SinogramPtr->delta_t)
      {
        if(ProfileCenterT <= ((ScannedObjectPtr->delta_z/2) - (SinogramPtr->delta_t/2)))
        {
          H_t[i] = SinogramPtr->delta_t;
        }
        else
        {
          H_t[i] = -ProfileCenterT + (ScannedObjectPtr->delta_z / 2) + SinogramPtr->delta_t / 2;
        }
        if(H_t[i] < 0)
        {
          H_t[i] = 0;
        }
      }
      else
      {
        if(ProfileCenterT <= SinogramPtr->delta_t/2 - ScannedObjectPtr->delta_z/2)
        {
          H_t[i] = ScannedObjectPtr->delta_z;
        }
        else
        {
          H_t[i] = -ProfileCenterT + (ScannedObjectPtr->delta_z/2) + SinogramPtr->delta_t/2;
        }

        if(H_t[i] < 0)
        {
          H_t[i] = 0;
        }

      }
      H_t[i] = H_t[i]/SinogramPtr->delta_t;/*Normalization*/
    }
	H_t[DETECTOR_RESPONSE_BINS] = 0;

  }

/*Generates the voxel line response from H_t*/
void storeVoxelLineResponse(Real_t* H_t,  AMatrixCol* VoxelLineResponse, ScannedObject* ScannedObjectPtr, Sinogram* SinogramPtr)
{
  
  Real_t ProfileThickness = 0.0;
  Real_t y = 0.0;
  Real_t t = 0.0;
  Real_t tmin;
  Real_t tmax;
  int32_t i_t, i, slice_index_min, slice_index_max;
  Real_t center_t,delta_t;
  int32_t index_delta_t;
  Real_t w3,w4;

  /*Storing the response along t-direction for each voxel line*/
  for (i = 0; i < (int32_t)ScannedObjectPtr->N_z; i++)
  {
    VoxelLineResponse[i].count = 0;
    
    y = ((Real_t)i + 0.5) * ScannedObjectPtr->delta_z + ScannedObjectPtr->z0;
    t = y;
    tmin = (t - ScannedObjectPtr->delta_z / 2) > SinogramPtr->T0 ? t - ScannedObjectPtr->delta_z / 2 : SinogramPtr->T0;
    tmax = (t + ScannedObjectPtr->delta_z / 2) <= SinogramPtr->TMax ? t + ScannedObjectPtr->delta_z / 2 : SinogramPtr->TMax;

    slice_index_min = (int32_t)(floor((tmin - SinogramPtr->T0) / SinogramPtr->delta_t));
    slice_index_max = (int32_t)(floor((tmax - SinogramPtr->T0) / SinogramPtr->delta_t));

    if(slice_index_min < 0)
    {
      slice_index_min = 0;
    }
    if(slice_index_max >= SinogramPtr->N_t)
    {
      slice_index_max = SinogramPtr->N_t - 1;
    }

    /*printf("%d %d\n",slice_index_min,slice_index_max);*/

    for (i_t = slice_index_min; i_t <= slice_index_max; i_t++)
    {
      center_t = ((Real_t)i_t + 0.5) * SinogramPtr->delta_t + SinogramPtr->T0;
      delta_t = fabs(center_t - t);
      index_delta_t = (int32_t)(floor(delta_t / SinogramPtr->OffsetT));
      if(index_delta_t < DETECTOR_RESPONSE_BINS)
      {
        w3 = delta_t - (Real_t)(index_delta_t) * SinogramPtr->OffsetT;
        w4 = ((Real_t)index_delta_t + 1) * SinogramPtr->OffsetT - delta_t;
        ProfileThickness = (w4 / SinogramPtr->OffsetT) * H_t[index_delta_t]+ (w3 / SinogramPtr->OffsetT) * H_t[index_delta_t+1];
    /*  ProfileThickness = (w4 / OffsetT) * detectorResponse->d[0][uint16_t(floor(SinogramPtr->N_theta/2))][index_delta_t]
        + (w3 / OffsetT) * detectorResponse->d[0][uint16_t(floor(SinogramPtr->N_theta/2))][index_delta_t + 1 < DETECTOR_RESPONSE_BINS ? index_delta_t + 1 : DETECTOR_RESPONSE_BINS - 1];*/
    }
      else
      {
        ProfileThickness = 0;
      }

      if(ProfileThickness != 0) /*Store the response of this slice */
      {
        VoxelLineResponse[i].values[VoxelLineResponse[i].count] = ProfileThickness;/*VoxelLineResponse is normalized value*/
        VoxelLineResponse[i].index[VoxelLineResponse[i].count++] = i_t;
/*	printf ("i = %d, i_t = %d, value = %f\n", i, i_t, ProfileThickness);*/
      }
    }
  }
}



