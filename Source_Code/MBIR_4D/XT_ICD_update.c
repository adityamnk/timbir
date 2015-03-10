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

#include "XT_Constants.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "allocate.h"
#include "randlib.h"
#include <time.h>
#include "XT_AMatrix.h"
#include "XT_Profile.h"
#include "XT_Structures.h"
#include "XT_IOMisc.h"
#include "XT_NHICD.h"
#include "omp.h"
#include "XT_MPI.h"
#include <mpi.h>
#include "XT_VoxUpdate.h"
#include "XT_ForwardProject.h"
#include "XT_MPIIO.h"
#include "XT_Debug.h"
#include "XT_OffsetError.h"

/*computes the location of (i,j,k) th element in a 1D array*/
int32_t array_loc_1D (int32_t i, int32_t j, int32_t k, int32_t N_j, int32_t N_k)
{
  return (i*N_j*N_k + j*N_k + k);
}
/*finds the maximum in a array 'array_in' with number of elements being 'num'*/
int32_t find_max(int32_t* array_in, int32_t num)
{
  int32_t i, maxnum;
  maxnum = array_in[0];
  for (i=1; i<num; i++)
  if (array_in[i] > maxnum)
  maxnum = array_in[i];
  
  return(maxnum);
}
/*converts the value 'val' to hounsfield units and returns it*/
Real_t convert2Hounsfield (Real_t val)
{
  Real_t slope, c;
  
  slope=(HOUNSFIELD_WATER_MAP-HOUNSFIELD_AIR_MAP)/(WATER_MASS_ATT_COEFF*WATER_DENSITY-AIR_MASS_ATT_COEFF*AIR_DENSITY)/HFIELD_UNIT_CONV_CONST;
  c=-slope*(AIR_MASS_ATT_COEFF*AIR_DENSITY*HFIELD_UNIT_CONV_CONST);
  
  return (slope*val + c);
}
/*Computes the qGGMRF spatial prior cost value at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Spatial_Value(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return ((pow(fabs(delta),MRF_Q)/TomoInputsPtr->Sigma_S_Q)/(ScannedObjectPtr->C_S + pow(fabs(delta),MRF_Q - MRF_P)/TomoInputsPtr->Sigma_S_Q_P));
}
/*Computes the qGGMRF temporal prior cost value at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Temporal_Value(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return ((pow(fabs(delta),MRF_Q)/TomoInputsPtr->Sigma_T_Q)/(ScannedObjectPtr->C_T + pow(fabs(delta),MRF_Q - MRF_P)/TomoInputsPtr->Sigma_T_Q_P));
}
/*Computes the qGGMRF spatial prior derivative at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Spatial_Derivative(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  Real_t temp1,temp2,temp3;
  temp1=pow(fabs(delta),MRF_Q - MRF_P)/(TomoInputsPtr->Sigma_S_Q_P);
  temp2=pow(fabs(delta),MRF_Q - 1);
  temp3 = ScannedObjectPtr->C_S + temp1;
  if(delta < 0)
  return ((-1*temp2/(temp3*TomoInputsPtr->Sigma_S_Q))*(MRF_Q - ((MRF_Q-MRF_P)*temp1)/(temp3)));
  else
  {
    return ((temp2/(temp3*TomoInputsPtr->Sigma_S_Q))*(MRF_Q - ((MRF_Q-MRF_P)*temp1)/(temp3)));
  }
}
/*Computes the qGGMRF temporal prior derivative at delta = x_i - x_j. i & j being the voxel and its neighbor*/
Real_t CE_QGGMRF_Temporal_Derivative(Real_t delta, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  Real_t temp1,temp2,temp3;
  temp1 = pow(fabs(delta),MRF_Q - MRF_P)/(TomoInputsPtr->Sigma_T_Q_P);
  temp2 = pow(fabs(delta),MRF_Q - 1);
  temp3 = ScannedObjectPtr->C_T + temp1;
  if(delta < 0)
  return ((-1*temp2/(temp3*TomoInputsPtr->Sigma_T_Q))*(MRF_Q - ((MRF_Q-MRF_P)*temp1)/(temp3)));
  else
  {
    return ((temp2/(temp3*TomoInputsPtr->Sigma_T_Q))*(MRF_Q - ((MRF_Q-MRF_P)*temp1)/(temp3)));
  }
}
/*Computes the qGGMRF spatial prior second derivative at delta = 0*/
Real_t CE_QGGMRF_Spatial_SecondDerivative(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return MRF_Q/(TomoInputsPtr->Sigma_S_Q*ScannedObjectPtr->C_S);
}
/*Computes the qGGMRF spatial prior second derivative at delta = 0*/
Real_t CE_QGGMRF_Temporal_SecondDerivative(ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
{
  return MRF_Q/(TomoInputsPtr->Sigma_T_Q*ScannedObjectPtr->C_T);
}
/*Computes the voxel update and returns it. V is the present value of voxel.
THETA1 and THETA2 are the values used in voxel update. Spatial_Nhood and Time_Nhood gives the
values of voxels in the neighborhood of V. Time_BDFlag and Spatial_BDFlag are masks which determine
whether a neighbor should be included in the neighorhood or not.*/
Real_t CE_FunctionalSubstitution(Real_t V, Real_t THETA1, Real_t THETA2, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_t Spatial_Nhood[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], Real_t Time_Nhood[NHOOD_TIME_MAXDIM-1], bool Spatial_BDFlag[NHOOD_Y_MAXDIM][NHOOD_X_MAXDIM][NHOOD_Z_MAXDIM], bool Time_BDFlag[NHOOD_TIME_MAXDIM-1])
{
  Real_t u,temp1=0,temp2=0,temp_const,RefValue=0,Delta0;
  Real_t QGGMRF_Params;
  int32_t i,j,k;
  RefValue = V;
  /*Need to Loop this for multiple iterations of substitute function*/
  for (i=0; i < NHOOD_Y_MAXDIM; i++)
  for (j=0; j < NHOOD_X_MAXDIM; j++)
  for (k=0; k < NHOOD_Z_MAXDIM; k++)
  {
    if(Spatial_BDFlag[i][j][k] == true && (i != (NHOOD_Y_MAXDIM-1)/2 || j != (NHOOD_X_MAXDIM-1)/2 || k != (NHOOD_Z_MAXDIM-1)/2))
    {
      Delta0 = (RefValue - Spatial_Nhood[i][j][k]);
      if(Delta0 != 0)
      QGGMRF_Params = CE_QGGMRF_Spatial_Derivative(Delta0,ScannedObjectPtr,TomoInputsPtr)/(Delta0);
      else {
        QGGMRF_Params = CE_QGGMRF_Spatial_SecondDerivative(ScannedObjectPtr,TomoInputsPtr);
      }
      temp_const = TomoInputsPtr->Spatial_Filter[i][j][k]*QGGMRF_Params;
      temp1 += temp_const*Spatial_Nhood[i][j][k];
      temp2 += temp_const;
    }
  }
  for (i=0; i < NHOOD_TIME_MAXDIM - 1; i++)
  {
    if(Time_BDFlag[i] == true)
    {
      Delta0 = (RefValue - Time_Nhood[i]);
      if(Delta0 != 0)
      QGGMRF_Params = CE_QGGMRF_Temporal_Derivative(Delta0,ScannedObjectPtr,TomoInputsPtr)/(Delta0);
      else {
        QGGMRF_Params = CE_QGGMRF_Temporal_SecondDerivative(ScannedObjectPtr,TomoInputsPtr);
      }
      
      temp_const = TomoInputsPtr->Time_Filter[0]*QGGMRF_Params;
      temp1 += temp_const*Time_Nhood[i];
      temp2 += temp_const;
    }
  }
  
  u=(temp1+ (THETA2*V) - THETA1)/(temp2 + THETA2);
  
  RefValue = RefValue + TomoInputsPtr->alpha*(u-RefValue);
  #ifdef POSITIVITY_CONSTRAINT
  if (RefValue <= 0)
  RefValue = 0;
  #endif
  return RefValue;
}

/*computes the value of cost function. 'ErrorSino' is the error sinogram*/
Real_t computeCost(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino)
{
  Real_t cost=0,temp=0, forward=0, prior=0;
  Real_t delta;
  int32_t i,j,k,p,N_z;
  bool j_minus, k_minus, i_plus, j_plus, k_plus, p_plus;
  
  #pragma omp parallel for private(j, k, temp) reduction(+:cost)
  for (i = 0; i < SinogramPtr->N_p; i++)
  for (j = 0; j < SinogramPtr->N_r; j++)
  for (k = 0; k < SinogramPtr->N_t; k++)
  {
    temp = ErrorSino[i][j][k] * sqrt(TomoInputsPtr->Weight[i][j][k]);
    if (SinogramPtr->ProjSelect[i][j][k] == true)
    temp = temp*temp;
    else
    temp = 2.0*TomoInputsPtr->ErrorSinoDelta*TomoInputsPtr->ErrorSinoThresh*fabs(temp) + TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoThresh*(1.0-2.0*TomoInputsPtr->ErrorSinoDelta);
    cost += temp;
  }
  cost /= 2.0;
  /*When computing the cost of the prior term it is important to make sure that you don't include the cost of any pair of neighbors more than once. In this code, a certain sense of causality is used to compute the cost. We also assume that the weghting kernel given by 'Filter' is symmetric. Let i, j and k correspond to the three dimensions. If we go forward to i+1, then all neighbors at j-1, j, j+1, k+1, k, k-1 are to be considered. However, if for the same i, if we go forward to j+1, then all k-1, k, and k+1 should be considered. For same i and j, only the neighbor at k+1 is considred.*/
  temp = 0;
  N_z = ScannedObjectPtr->N_z + 2;
  if (TomoInputsPtr->node_rank == TomoInputsPtr->node_num-1)
  N_z = ScannedObjectPtr->N_z + 1;
  #pragma omp parallel for private(delta, p, j, k, j_minus, k_minus, p_plus, i_plus, j_plus, k_plus) reduction(+:temp)
  for (i = 0; i < ScannedObjectPtr->N_time; i++)
  for (p = 1; p < ScannedObjectPtr->N_z + 1; p++)
  for (j = 0; j < ScannedObjectPtr->N_y; j++)
  {
    for (k = 0; k < ScannedObjectPtr->N_x; k++)
    {
      j_minus = (j - 1 >= 0)? true : false;
      k_minus = (k - 1 >= 0)? true : false;
      
      p_plus = (p + 1 < N_z)? true : false;
      i_plus = (i + 1 < ScannedObjectPtr->N_time)? true : false;
      j_plus = (j + 1 < ScannedObjectPtr->N_y)? true : false;
      k_plus = (k + 1 < ScannedObjectPtr->N_x)? true : false;
      
      if(k_plus == true) {
        delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j][k + 1]);
        temp += TomoInputsPtr->Spatial_Filter[1][1][2] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
      }
      if(j_plus == true) {
        if(k_minus == true) {
          delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j + 1][k - 1]);
          temp += TomoInputsPtr->Spatial_Filter[1][2][0] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
        }
        delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j + 1][k]);
        temp += TomoInputsPtr->Spatial_Filter[1][2][1] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
        if(k_plus == true) {
          delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p][j + 1][k + 1]);
          temp += TomoInputsPtr->Spatial_Filter[1][2][2] * CE_QGGMRF_Spatial_Value(delta,ScannedObjectPtr,TomoInputsPtr);
        }
      }
      if (p_plus == true)
      {
        if(j_minus == true)
        {
          delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j - 1][k];
          temp += TomoInputsPtr->Spatial_Filter[2][0][1] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
        }
        
        delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p+1][j][k];
        temp += TomoInputsPtr->Spatial_Filter[2][1][1] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
        if(j_plus == true)
        {
          delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p+1][j + 1][k];
          temp += TomoInputsPtr->Spatial_Filter[2][2][1] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
        }
        if(j_minus == true)
        {
          if(k_minus == true)
          {
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j - 1][k - 1];
            temp += TomoInputsPtr->Spatial_Filter[2][0][0] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }
          if(k_plus == true)
          {
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j - 1][k + 1];
            temp += TomoInputsPtr->Spatial_Filter[2][0][2] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }
        }
        if(k_minus == true)
        {
          delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j][k - 1];
          temp += TomoInputsPtr->Spatial_Filter[2][1][0] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
        }
        if(j_plus == true)
        {
          if(k_minus == true)
          {
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j + 1][k - 1];
            temp += TomoInputsPtr->Spatial_Filter[2][2][0] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }
          if(k_plus == true)
          {
            delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j + 1][k + 1];
            temp += TomoInputsPtr->Spatial_Filter[2][2][2] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
          }
        }
        if(k_plus == true)
        {
          delta = ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i][p + 1][j][k + 1];
          temp += TomoInputsPtr->Spatial_Filter[2][1][2] * CE_QGGMRF_Spatial_Value(delta, ScannedObjectPtr, TomoInputsPtr);
        }
      }
      if(i_plus == true) {
        delta = (ScannedObjectPtr->Object[i][p][j][k] - ScannedObjectPtr->Object[i+1][p][j][k]);
        temp += TomoInputsPtr->Time_Filter[0] * CE_QGGMRF_Temporal_Value(delta,ScannedObjectPtr,TomoInputsPtr);
      }
    }
  }
  /*Use MPI reduction operation to add the forward and prior costs from all nodes*/
  MPI_Reduce(&cost, &forward, 1, MPI_REAL_DATATYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&temp, &prior, 1, MPI_REAL_DATATYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (TomoInputsPtr->node_rank == 0)
  {
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Scaled error sino cost = %f\n",forward);
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Decrease in scaled error sino cost = %f\n",TomoInputsPtr->ErrorSino_Cost-forward);
    TomoInputsPtr->ErrorSino_Cost = forward;
    forward += (Real_t)TomoInputsPtr->node_num*(Real_t)SinogramPtr->N_p*(Real_t)SinogramPtr->N_r*(Real_t)SinogramPtr->N_t*log(TomoInputsPtr->var_est)/2;
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Forward cost = %f\n",forward);
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Prior cost = %f\n",prior);
    TomoInputsPtr->Forward_Cost = forward;
    TomoInputsPtr->Prior_Cost = prior;
    cost = forward + prior;
  }
  
  /*Broadcase the value of cost to all nodes*/
  MPI_Bcast(&cost, 1, MPI_REAL_DATATYPE, 0, MPI_COMM_WORLD);
  return cost;
}


/*Upsamples the (N_time x N_z x N_y x N_x) size 'Init' by a factor of 2 along the x-y plane and stores it in 'Object'*/
void upsample_bilinear_2D (Real_arr_t**** Object, Real_arr_t**** Init, int32_t N_time, int32_t N_z, int32_t N_y, int32_t N_x)
{
  int32_t i, j, k, m;
  Real_arr_t **buffer;
  
  #pragma omp parallel for private(buffer, m, j, k)
  for (i=0; i < N_time; i++)
  for (m=0; m < N_z; m++)
  {
    buffer = (Real_arr_t**)multialloc(sizeof(Real_arr_t), 2, N_y, 2*N_x);
    for (j=0; j < N_y; j++){
      buffer[j][0] = Init[i][m][j][0];
      buffer[j][1] = (3.0*Init[i][m][j][0] + Init[i][m][j][1])/4.0;
      buffer[j][2*N_x - 1] = Init[i][m][j][N_x - 1];
      buffer[j][2*N_x - 2] = (Init[i][m][j][N_x - 2] + 3.0*Init[i][m][j][N_x - 1])/4.0;
      for (k=1; k < N_x - 1; k++){
        buffer[j][2*k] = (Init[i][m][j][k-1] + 3.0*Init[i][m][j][k])/4.0;
        buffer[j][2*k + 1] = (3.0*Init[i][m][j][k] + Init[i][m][j][k+1])/4.0;
      }
    }
    for (k=0; k < 2*N_x; k++){
      Object[i][m][0][k] = buffer[0][k];
      Object[i][m][1][k] = (3.0*buffer[0][k] + buffer[1][k])/4.0;
      Object[i][m][2*N_y-1][k] = buffer[N_y-1][k];
      Object[i][m][2*N_y-2][k] = (buffer[N_y-2][k] + 3.0*buffer[N_y-1][k])/4.0;
    }
    for (j=1; j<N_y-1; j++){
      for (k=0; k<2*N_x; k++){
        Object[i][m][2*j][k] = (buffer[j-1][k] + 3.0*buffer[j][k])/4.0;
        Object[i][m][2*j + 1][k] = (3*buffer[j][k] + buffer[j+1][k])/4.0;
      }
    }
    multifree(buffer,2);
  }
}
/*Upsamples the (N_z x N_y x N_x) size 'Init' by a factor of 2 along the x-y plane and stores it in 'Object'*/
void upsample_object_bilinear_2D (Real_arr_t*** Object, Real_arr_t*** Init, int32_t N_z, int32_t N_y, int32_t N_x)
{
  int32_t j, k, slice;
  Real_arr_t **buffer;
  
  
  buffer = (Real_arr_t**)multialloc(sizeof(Real_arr_t), 2, N_y, 2*N_x);
  for (slice=0; slice < N_z; slice++){
    for (j=0; j < N_y; j++){
      buffer[j][0] = Init[slice][j][0];
      buffer[j][1] = (3.0*Init[slice][j][0] + Init[slice][j][1])/4.0;
      buffer[j][2*N_x - 1] = Init[slice][j][N_x - 1];
      buffer[j][2*N_x - 2] = (Init[slice][j][N_x - 2] + 3.0*Init[slice][j][N_x - 1])/4.0;
      for (k=1; k < N_x - 1; k++){
        buffer[j][2*k] = (Init[slice][j][k-1] + 3.0*Init[slice][j][k])/4.0;
        buffer[j][2*k + 1] = (3.0*Init[slice][j][k] + Init[slice][j][k+1])/4.0;
      }
    }
    for (k=0; k < 2*N_x; k++){
      Object[slice+1][0][k] = buffer[0][k];
      Object[slice+1][1][k] = (3.0*buffer[0][k] + buffer[1][k])/4.0;
      Object[slice+1][2*N_y-1][k] = buffer[N_y-1][k];
      Object[slice+1][2*N_y-2][k] = (buffer[N_y-2][k] + 3.0*buffer[N_y-1][k])/4.0;
    }
    for (j=1; j<N_y-1; j++){
      for (k=0; k<2*N_x; k++){
        Object[slice+1][2*j][k] = (buffer[j-1][k] + 3.0*buffer[j][k])/4.0;
        Object[slice+1][2*j + 1][k] = (3*buffer[j][k] + buffer[j+1][k])/4.0;
      }
    }
  }
  multifree(buffer,2);
}

void upsample_bilinear_3D (Real_arr_t**** Object, Real_arr_t**** Init, int32_t N_time, int32_t N_z, int32_t N_y, int32_t N_x)
{
  int32_t i, j, k, slice;
  Real_t ***buffer2D, ***buffer3D;
  
  #pragma omp parallel for private(buffer2D, buffer3D, slice, j, k)
  for (i=0; i < N_time; i++)
  {
	  buffer2D = (Real_t***)multialloc(sizeof(Real_t), 3, N_z, N_y, 2*N_x);
	  buffer3D = (Real_t***)multialloc(sizeof(Real_t), 3, N_z, 2*N_y, 2*N_x);
	  for (slice=0; slice < N_z; slice++){
	    for (j=0; j < N_y; j++){
	      buffer2D[slice][j][0] = Init[i][slice][j][0];
	      buffer2D[slice][j][1] = (3.0*Init[i][slice][j][0] + Init[i][slice][j][1])/4.0;
	      buffer2D[slice][j][2*N_x - 1] = Init[i][slice][j][N_x - 1];
	      buffer2D[slice][j][2*N_x - 2] = (Init[i][slice][j][N_x - 2] + 3.0*Init[i][slice][j][N_x - 1])/4.0;
	      for (k=1; k < N_x - 1; k++){
        	buffer2D[slice][j][2*k] = (Init[i][slice][j][k-1] + 3.0*Init[i][slice][j][k])/4.0;
        	buffer2D[slice][j][2*k + 1] = (3.0*Init[i][slice][j][k] + Init[i][slice][j][k+1])/4.0;
      	     }
    	    }
    	    for (k=0; k < 2*N_x; k++){
      		buffer3D[slice][0][k] = buffer2D[slice][0][k];
      		buffer3D[slice][1][k] = (3.0*buffer2D[slice][0][k] + buffer2D[slice][1][k])/4.0;
      		buffer3D[slice][2*N_y-1][k] = buffer2D[slice][N_y-1][k];
      		buffer3D[slice][2*N_y-2][k] = (buffer2D[slice][N_y-2][k] + 3.0*buffer2D[slice][N_y-1][k])/4.0;
    	    }
    		for (j=1; j<N_y-1; j++)
    		for (k=0; k<2*N_x; k++){
      			buffer3D[slice][2*j][k] = (buffer2D[slice][j-1][k] + 3.0*buffer2D[slice][j][k])/4.0;
      			buffer3D[slice][2*j + 1][k] = (3*buffer2D[slice][j][k] + buffer2D[slice][j+1][k])/4.0;
    		}
  	   }
  
  	for (j=0; j<2*N_y; j++)
  	for (k=0; k<2*N_x; k++){
    		Object[i][0][j][k] = buffer3D[0][j][k];
    		Object[i][1][j][k] = (3.0*buffer3D[0][j][k] + buffer3D[1][j][k])/4.0;
    		Object[i][2*N_z-1][j][k] = buffer3D[N_z-1][j][k];
    		Object[i][2*N_z-2][j][k] = (3.0*buffer3D[N_z-1][j][k] + buffer3D[N_z-2][j][k])/4.0;
  	}
  
  	for (slice=1; slice < N_z-1; slice++)
  	for (j=0; j<2*N_y; j++)
  	for (k=0; k<2*N_x; k++){
    		Object[i][2*slice][j][k] = (buffer3D[slice-1][j][k] + 3.0*buffer3D[slice][j][k])/4.0;
    		Object[i][2*slice+1][j][k] = (3.0*buffer3D[slice][j][k] + buffer3D[slice+1][j][k])/4.0;
  	}
  
  	multifree(buffer2D,3);
  	multifree(buffer3D,3);
  }
}

/*'InitObject' intializes the Object to be reconstructed to either 0 or an interpolated version of the previous reconstruction. It is used in multi resolution reconstruction in which after every coarse resolution reconstruction the object should be intialized with an interpolated version of the reconstruction following which the object will be reconstructed at a finer resolution.*/
/*Upsamples the (N_time x N_z x N_y x N_x) size 'Init' by a factor of 2 along the in 3D x-y-z coordinates and stores it in 'Object'*/
void upsample_object_bilinear_3D (Real_arr_t*** Object, Real_arr_t*** Init, int32_t N_z, int32_t N_y, int32_t N_x)
{
  int32_t j, k, slice;
  Real_t ***buffer2D, ***buffer3D;
  
  buffer2D = (Real_t***)multialloc(sizeof(Real_t), 3, N_z, N_y, 2*N_x);
  buffer3D = (Real_t***)multialloc(sizeof(Real_t), 3, N_z, 2*N_y, 2*N_x);
  for (slice=0; slice < N_z; slice++){
    for (j=0; j < N_y; j++){
      buffer2D[slice][j][0] = Init[slice][j][0];
      buffer2D[slice][j][1] = (3.0*Init[slice][j][0] + Init[slice][j][1])/4.0;
      buffer2D[slice][j][2*N_x - 1] = Init[slice][j][N_x - 1];
      buffer2D[slice][j][2*N_x - 2] = (Init[slice][j][N_x - 2] + 3.0*Init[slice][j][N_x - 1])/4.0;
      for (k=1; k < N_x - 1; k++){
        buffer2D[slice][j][2*k] = (Init[slice][j][k-1] + 3.0*Init[slice][j][k])/4.0;
        buffer2D[slice][j][2*k + 1] = (3.0*Init[slice][j][k] + Init[slice][j][k+1])/4.0;
      }
    }
    for (k=0; k < 2*N_x; k++){
      buffer3D[slice][0][k] = buffer2D[slice][0][k];
      buffer3D[slice][1][k] = (3.0*buffer2D[slice][0][k] + buffer2D[slice][1][k])/4.0;
      buffer3D[slice][2*N_y-1][k] = buffer2D[slice][N_y-1][k];
      buffer3D[slice][2*N_y-2][k] = (buffer2D[slice][N_y-2][k] + 3.0*buffer2D[slice][N_y-1][k])/4.0;
    }
    for (j=1; j<N_y-1; j++)
    for (k=0; k<2*N_x; k++){
      buffer3D[slice][2*j][k] = (buffer2D[slice][j-1][k] + 3.0*buffer2D[slice][j][k])/4.0;
      buffer3D[slice][2*j + 1][k] = (3*buffer2D[slice][j][k] + buffer2D[slice][j+1][k])/4.0;
    }
  }
  
  for (j=0; j<2*N_y; j++)
  for (k=0; k<2*N_x; k++){
    Object[1][j][k] = buffer3D[0][j][k];
    Object[2][j][k] = (3.0*buffer3D[0][j][k] + buffer3D[1][j][k])/4.0;
    Object[2*N_z][j][k] = buffer3D[N_z-1][j][k];
    Object[2*N_z-1][j][k] = (3.0*buffer3D[N_z-1][j][k] + buffer3D[N_z-2][j][k])/4.0;
  }
  
  for (slice=1; slice < N_z-1; slice++)
  for (j=0; j<2*N_y; j++)
  for (k=0; k<2*N_x; k++){
    Object[2*slice+1][j][k] = (buffer3D[slice-1][j][k] + 3.0*buffer3D[slice][j][k])/4.0;
    Object[2*slice+2][j][k] = (3.0*buffer3D[slice][j][k] + buffer3D[slice+1][j][k])/4.0;
  }
  
  multifree(buffer2D,3);
  multifree(buffer3D,3);
}

/*randomly select the voxels lines which need to be updated along the x-y plane for each z-block and time slice*/
void randomly_select_x_y (ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, uint8_t*** Mask)
{
  int32_t i, j, num,n, Index, col, row, *Counter, ArraySize, block;
  ArraySize = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
  Counter = (int32_t*)get_spc(ArraySize, sizeof(int32_t));
  for (i=0; i<ScannedObjectPtr->N_time; i++)
  for (block=0; block<TomoInputsPtr->num_z_blocks; block++)
  {
    ArraySize = ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
    for (Index = 0; Index < ArraySize; Index++)
    Counter[Index] = Index;
    
    TomoInputsPtr->UpdateSelectNum[i][block] = 0;
    for (j=0; j<ScannedObjectPtr->N_x*ScannedObjectPtr->N_y; j++){
      Index = floor(random2() * ArraySize);
      Index = (Index == ArraySize)?ArraySize-1:Index;
      col = Counter[Index] % ScannedObjectPtr->N_x;
      row = Counter[Index] / ScannedObjectPtr->N_x;
      for (n = block*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks); n < (block+1)*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks); n++)
      if (Mask[i][row][col] == 1)
      {
        num = TomoInputsPtr->UpdateSelectNum[i][block];
        TomoInputsPtr->x_rand_select[i][block][num] = col;
        TomoInputsPtr->y_rand_select[i][block][num] = row;
        (TomoInputsPtr->UpdateSelectNum[i][block])++;
        break;
      }
      Counter[Index] = Counter[ArraySize - 1];
      ArraySize--;
    }
  }
  free(Counter);
}


/*'InitObject' intializes the Object to be reconstructed to either 0 or an interpolated version of the previous reconstruction. It is used in multi resolution reconstruction in which after every coarse resolution reconstruction the object should be intialized with an interpolated version of the reconstruction following which the object will be reconstructed at a finer resolution.
--initICD--
If 1, initializes the object to 0
If 2, the code uses bilinear interpolation to initialize the object if the previous reconstruction was at a lower resolution
The function also initializes the magnitude update map 'MagUpdateMap' from the previous coarser resolution
reconstruction. */
int32_t initObject (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t**** MagUpdateMap)
{
  char object_file[100];
  int dimTiff[4];
  int32_t i, j, k, l, size, flag = 0;
  Real_arr_t ***Init, ****UpMapInit;
  
  for (i = 0; i < ScannedObjectPtr->N_time; i++)
  for (j = 0; j < ScannedObjectPtr->N_z; j++)
  for (k = 0; k < ScannedObjectPtr->N_y; k++)
  for (l = 0; l < ScannedObjectPtr->N_x; l++)
  	ScannedObjectPtr->Object[i][j+1][k][l] = OBJECT_INIT_VAL;
  
  if (TomoInputsPtr->initICD > 3 || TomoInputsPtr->initICD < 0){
	sentinel(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "ERROR: initICD value not recognized.\n");
  }
  else if (TomoInputsPtr->initICD == 1)
  {
	size = ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
      	for (i = 0; i < ScannedObjectPtr->N_time; i++)
      	{
        	sprintf(object_file, "%s_time_%d", OBJECT_FILENAME,i);
		if (read_SharedBinFile_At (object_file, &(ScannedObjectPtr->Object[i][1][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) flag = -1;
      	}
	if (TomoInputsPtr->initMagUpMap == 1)
      	{
		size = ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
		if (read_SharedBinFile_At (MAG_UPDATE_FILENAME, &(MagUpdateMap[0][0][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) flag = -1;
      	}
  }
  else if (TomoInputsPtr->initICD == 2 || TomoInputsPtr->initICD == 3)
  {
      	if (TomoInputsPtr->initICD == 3)
      	{
        	Init = (Real_arr_t***)multialloc(sizeof(Real_arr_t), 3, ScannedObjectPtr->N_z/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
	        check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Interpolating object using 3D bilinear interpolation.\n");
	        for (i = 0; i < ScannedObjectPtr->N_time; i++)
        	{
         		 sprintf(object_file, "%s_time_%d", OBJECT_FILENAME, i);
			 size = ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x/8;
			 if (read_SharedBinFile_At (object_file, &(Init[0][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) flag = -1;
          		 upsample_object_bilinear_3D (ScannedObjectPtr->Object[i], Init, ScannedObjectPtr->N_z/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
        	}
       		multifree(Init,3);
        	check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Done with interpolating object using 3D bilinear interpolation.\n");
      	}
	else
      	{
        	Init = (Real_arr_t***)multialloc(sizeof(Real_arr_t), 3, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
	        check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Interpolating object using 2D bilinear interpolation.\n");
        	for (i = 0; i < ScannedObjectPtr->N_time; i++)
        	{
         		sprintf(object_file, "%s_time_%d", OBJECT_FILENAME,i);
	  		size = ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x/4;
	  		if (read_SharedBinFile_At (object_file, &(Init[0][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) flag = -1;
          		upsample_object_bilinear_2D (ScannedObjectPtr->Object[i], Init, ScannedObjectPtr->N_z, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
        	}
        	multifree(Init,3);
        	check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Done with interpolating object using 2D bilinear interpolation.\n");
      	}
        if (TomoInputsPtr->initMagUpMap == 1)
        {
          	if (TomoInputsPtr->prevnum_z_blocks == TomoInputsPtr->num_z_blocks)
          	{	
			UpMapInit = (Real_arr_t****)multialloc(sizeof(Real_arr_t), 4, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
			size = ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x/4;
			if (read_SharedBinFile_At (MAG_UPDATE_FILENAME, &(UpMapInit[0][0][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) flag = -1;
          		check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Interpolating magnitude update map using 2D bilinear interpolation.\n");
          		upsample_bilinear_2D (MagUpdateMap, UpMapInit, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
          		multifree(UpMapInit,4);
	  	}
		else if (TomoInputsPtr->prevnum_z_blocks == TomoInputsPtr->num_z_blocks/2)
	  	{
			UpMapInit = (Real_arr_t****)multialloc(sizeof(Real_arr_t), 4, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
			size = ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x/8;
			if (read_SharedBinFile_At (MAG_UPDATE_FILENAME, &(UpMapInit[0][0][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) flag = -1; 
          		check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Interpolating magnitude update map using 3D bilinear interpolation.\n");
			upsample_bilinear_3D (MagUpdateMap, UpMapInit, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks/2, ScannedObjectPtr->N_y/2, ScannedObjectPtr->N_x/2);
          		multifree(UpMapInit,4);
	  	}
	  	else
	  	{
			check_warn(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Number of axial blocks is incompatible with previous stage of multi-resolution.\n");
			check_warn(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Initializing the multi-resolution map to zeros.\n");
	  	}	
          }
      }
  
  	dimTiff[0] = ScannedObjectPtr->N_time; dimTiff[1] = TomoInputsPtr->num_z_blocks; dimTiff[2] = ScannedObjectPtr->N_y; dimTiff[3] = ScannedObjectPtr->N_x;
  	sprintf(object_file, "%s_n%d", MAG_UPDATE_FILENAME, TomoInputsPtr->node_rank);
  	if (TomoInputsPtr->Write2Tiff == 1)
  		if (WriteMultiDimArray2Tiff (object_file, dimTiff, 0, 1, 2, 3, &(MagUpdateMap[0][0][0][0]), 0, TomoInputsPtr->debug_file_ptr))
			flag = -1;
  
  	for (i = 0; i < ScannedObjectPtr->N_time; i++)
  	{
		sprintf (object_file, "%s_n%d", INIT_OBJECT_FILENAME, TomoInputsPtr->node_rank);
	    	sprintf (object_file, "%s_time_%d", object_file, i);
    		dimTiff[0] = 1; dimTiff[1] = ScannedObjectPtr->N_z; dimTiff[2] = ScannedObjectPtr->N_y; dimTiff[3] = ScannedObjectPtr->N_x;
    		if (TomoInputsPtr->Write2Tiff == 1)
    			if (WriteMultiDimArray2Tiff (object_file, dimTiff, 0, 1, 2, 3, &(ScannedObjectPtr->Object[i][1][0][0]), 0, TomoInputsPtr->debug_file_ptr))
				flag = -1;
  	}
	
	return (flag);
error:
	return (-1);
}


/*'initErrorSinogram' is used to initialize the error sinogram before start of ICD. It computes e = y - Ax - d. Ax is computed by forward projecting the obkject x.*/
int32_t initErrorSinogam (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t** DetectorResponse, Real_arr_t*** ErrorSino/*, AMatrixCol* VoxelLineResponse*/)
{
  Real_t pixel, avg=0;
  int32_t dimTiff[4], i, j, k, p, sino_idx, slice, flag = 0;
  AMatrixCol* AMatrixPtr = (AMatrixCol*)get_spc(ScannedObjectPtr->N_time, sizeof(AMatrixCol));
  uint8_t AvgNumXElements = (uint8_t)ceil(3*ScannedObjectPtr->delta_xy/SinogramPtr->delta_r);
  char error_file[100] = "error_sinogram";
  sprintf(error_file, "%s_n%d", error_file, TomoInputsPtr->node_rank);
  for (i = 0; i < ScannedObjectPtr->N_time; i++)
  {
    AMatrixPtr[i].values = (Real_t*)get_spc(AvgNumXElements, sizeof(Real_t));
    AMatrixPtr[i].index = (int32_t*)get_spc(AvgNumXElements, sizeof(int32_t));
  }
  memset(&(ErrorSino[0][0][0]), 0, SinogramPtr->N_p*SinogramPtr->N_t*SinogramPtr->N_r*sizeof(Real_arr_t));
  #pragma omp parallel for private(j, k, p, sino_idx, slice, pixel)
  for (i=0; i<ScannedObjectPtr->N_time; i++)
  {
    for (j=0; j<ScannedObjectPtr->N_y; j++)
    {
      for (k=0; k<ScannedObjectPtr->N_x; k++){
        for (p=0; p<ScannedObjectPtr->ProjNum[i]; p++){
          sino_idx = ScannedObjectPtr->ProjIdxPtr[i][p];
          calcAMatrixColumnforAngle(SinogramPtr, ScannedObjectPtr, DetectorResponse, &(AMatrixPtr[i]), j, k, sino_idx);
          for (slice=0; slice<ScannedObjectPtr->N_z; slice++){
            /*	printf("count = %d, idx = %d, val = %f\n", VoxelLineResponse[slice].count, VoxelLineResponse[slice].index[0], VoxelLineResponse[slice].values[0]);*/
            pixel = ScannedObjectPtr->Object[i][slice+1][j][k]; /*slice+1 to account for extra z slices required for MPI*/
            forward_project_voxel (SinogramPtr, pixel, ErrorSino, &(AMatrixPtr[i])/*, &(VoxelLineResponse[slice])*/, sino_idx, slice);
          }
        }
      }
    }
  }
  
  #pragma omp parallel for private(j, k) reduction(+:avg)
  for(i=0; i < SinogramPtr->N_p; i++)
  for (j = 0; j < SinogramPtr->N_r; j++)
  for (k = 0; k < SinogramPtr->N_t; k++)
  {
    ErrorSino[i][j][k] = SinogramPtr->Projection[i][j][k] - ErrorSino[i][j][k] - SinogramPtr->ProjOffset[j][k];
    if (fabs(ErrorSino[i][j][k]*sqrt(TomoInputsPtr->Weight[i][j][k])) < TomoInputsPtr->ErrorSinoThresh)
	SinogramPtr->ProjSelect[i][j][k] = true;
    else
    	SinogramPtr->ProjSelect[i][j][k] = false;
    /*	if (ErrorSino[i][j][k]*sqrt(TomoInputsPtr->Weight[i][j][k]) < -30)
    TomoInputsPtr->Weight[i][j][k] = 0;*/
    avg+=ErrorSino[i][j][k];
  }
  avg = avg/(SinogramPtr->N_r*SinogramPtr->N_t*SinogramPtr->N_p);
  check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Average of error sinogram in node %d is %f\n", TomoInputsPtr->node_rank, avg);
  
  dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
  if (TomoInputsPtr->Write2Tiff == 1)
  	flag = WriteMultiDimArray2Tiff (error_file, dimTiff, 0, 3, 1, 2, &(ErrorSino[0][0][0]), 0, TomoInputsPtr->debug_file_ptr);
  for (i = 0; i < ScannedObjectPtr->N_time; i++)
  {
    free(AMatrixPtr[i].values);
    free(AMatrixPtr[i].index);
  }
  free (AMatrixPtr);
  multifree(SinogramPtr->Projection,3);
  return (flag);
}


/*Updates the variance parameter \sigma*/
void update_variance_parameter (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino)
{
  int32_t k, i, j;
  Real_t temp_acc = 0, temp = 0;
  
  #pragma omp parallel for private(i, j, temp) reduction(+:temp_acc)
  for (k = 0; k < SinogramPtr->N_p; k++)
  for (i = 0; i < SinogramPtr->N_r; i++)
  for (j = 0; j < SinogramPtr->N_t; j++)
  {
    TomoInputsPtr->Weight[k][i][j] = TomoInputsPtr->Weight[k][i][j]*TomoInputsPtr->var_est;
    if (SinogramPtr->ProjSelect[k][i][j] == true)
    temp = ErrorSino[k][i][j]*ErrorSino[k][i][j]*TomoInputsPtr->Weight[k][i][j];
    else
    temp = fabs(ErrorSino[k][i][j])*TomoInputsPtr->ErrorSinoDelta*TomoInputsPtr->ErrorSinoThresh*sqrt(TomoInputsPtr->Weight[k][i][j]*TomoInputsPtr->var_est);
    temp_acc += temp;
  }
  
  MPI_Allreduce(&temp_acc, &temp, 1, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
  TomoInputsPtr->var_est = temp/((Real_t)TomoInputsPtr->node_num*(Real_t)SinogramPtr->N_p*(Real_t)SinogramPtr->N_r*(Real_t)SinogramPtr->N_t);
  #pragma omp parallel for private(i, j)
  for (k = 0; k < SinogramPtr->N_p; k++)
  for (i = 0; i < SinogramPtr->N_r; i++)
  for (j = 0; j < SinogramPtr->N_t; j++)
  {
    TomoInputsPtr->Weight[k][i][j] /= TomoInputsPtr->var_est;
    if (fabs(ErrorSino[k][i][j]*sqrt(TomoInputsPtr->Weight[k][i][j])) < TomoInputsPtr->ErrorSinoThresh)
    SinogramPtr->ProjSelect[k][i][j] = true;
    else
    SinogramPtr->ProjSelect[k][i][j] = false;
  }
}

void update_d_offset_rect_patch_constraint (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino)
{
    Real_t sign, **b, **Lambda, temp;
    Real_arr_t **x;
    int32_t i, j, k;
   
    b = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    Lambda = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    x = (Real_arr_t**)multialloc(sizeof(Real_arr_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    memset(&(b[0][0]), 0, SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_t));
    memset(&(Lambda[0][0]), 0, SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_t));
    memset(&(x[0][0]), 0, SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_arr_t));

    #pragma omp parallel for collapse(2) private(i, j, k, temp, sign)
    for (i = 0; i < SinogramPtr->N_r; i++)
    {
    	for (j = 0; j < SinogramPtr->N_t; j++)
    	{
     		b[i][j] = 0;
      		Lambda[i][j] = 0;
      		for (k = 0; k < SinogramPtr->N_p; k++)
      		{
        		temp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[k][i][j]);
        		if (SinogramPtr->ProjSelect[k][i][j] == true)
        		{
          			Lambda[i][j] += TomoInputsPtr->Weight[k][i][j];
          			b[i][j] += (ErrorSino[k][i][j] + SinogramPtr->ProjOffset[i][j])*TomoInputsPtr->Weight[k][i][j];
        		}
        		else
        		{
	  			sign = (ErrorSino[k][i][j] > 0) - (ErrorSino[k][i][j] < 0);
          			Lambda[i][j] += temp/fabs(ErrorSino[k][i][j]);
          			b[i][j] += (ErrorSino[k][i][j] + SinogramPtr->ProjOffset[i][j])*temp/fabs(ErrorSino[k][i][j]);
        		}
      		}
    	}
    }

    constrained_quad_opt (Lambda, b, SinogramPtr->off_constraint, x, SinogramPtr->N_r, SinogramPtr->N_t, SinogramPtr->off_constraint_num, TomoInputsPtr);
      
    #pragma omp parallel for collapse(3) private(i, j, k)
    for (k = 0; k < SinogramPtr->N_p; k++)
    {
    	for (i = 0; i < SinogramPtr->N_r; i++)
    	{	
    		for (j = 0; j < SinogramPtr->N_t; j++)
    		{
			ErrorSino[k][i][j] += SinogramPtr->ProjOffset[i][j] - x[i][j];
        		if (fabs(ErrorSino[k][i][j]*sqrt(TomoInputsPtr->Weight[k][i][j])) < TomoInputsPtr->ErrorSinoThresh)
        			SinogramPtr->ProjSelect[k][i][j] = true;
        		else
        			SinogramPtr->ProjSelect[k][i][j] = false;
      		}
    	} 
    }

    memcpy(&(SinogramPtr->ProjOffset[0][0]),&(x[0][0]),SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_arr_t));
    
    multifree(b,2);
    multifree(Lambda,2);
    multifree(x,2);  
}
  
/*Updates the projection offset error parameter d_i*/
void update_d_offset_zero_mean_constraint (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino)
{
    Real_t sign, **numerator, num_sum = 0, temp, **denominator, den_sum = 0, gamma = 0;
    int32_t i, j, k;
   
    numerator = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    denominator = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    #pragma omp parallel for private(j, k, temp, sign) reduction(+:num_sum, den_sum)
    for (i = 0; i < SinogramPtr->N_r; i++)
    for (j = 0; j < SinogramPtr->N_t; j++)
    {
      numerator[i][j] = 0;
      denominator[i][j] = 0;
      for (k = 0; k < SinogramPtr->N_p; k++)
      {
        temp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[k][i][j]);
        if (SinogramPtr->ProjSelect[k][i][j] == true)
        {
          numerator[i][j] += ErrorSino[k][i][j]*TomoInputsPtr->Weight[k][i][j];
          denominator[i][j] += TomoInputsPtr->Weight[k][i][j];
        }
        else
        {
	  sign = (ErrorSino[k][i][j] > 0) - (ErrorSino[k][i][j] < 0);
          numerator[i][j] += temp*sign;
          denominator[i][j] += temp/fabs(ErrorSino[k][i][j]);
        }
      }
      num_sum += SinogramPtr->ProjOffset[i][j] + (numerator[i][j]/denominator[i][j]);
      den_sum += 1.0/denominator[i][j];
    }
    gamma = num_sum/den_sum;

    #pragma omp parallel for private(j, k)
    for (i = 0; i < SinogramPtr->N_r; i++)
    for (j = 0; j < SinogramPtr->N_t; j++)
    {
      SinogramPtr->ProjOffset[i][j] = SinogramPtr->ProjOffset[i][j] + (numerator[i][j]-gamma)/denominator[i][j];
      for (k = 0; k < SinogramPtr->N_p; k++)
      {
	ErrorSino[k][i][j] -= (numerator[i][j]-gamma)/denominator[i][j];
        if (fabs(ErrorSino[k][i][j]*sqrt(TomoInputsPtr->Weight[k][i][j])) < TomoInputsPtr->ErrorSinoThresh)
        SinogramPtr->ProjSelect[k][i][j] = true;
        else
        SinogramPtr->ProjSelect[k][i][j] = false;
      }
    }
    
   multifree(numerator,2);
   multifree(denominator,2);
  }

/*Updates the projection offset error parameter d_i*/
void update_d_offset_unconstrained (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino)
{
    Real_t sign, **numerator, temp, **denominator;
    int32_t i, j, k;
   
    numerator = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    denominator = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    #pragma omp parallel for private(j, k, temp, sign)
    for (i = 0; i < SinogramPtr->N_r; i++)
    for (j = 0; j < SinogramPtr->N_t; j++)
    {
      numerator[i][j] = 0;
      denominator[i][j] = 0;
      for (k = 0; k < SinogramPtr->N_p; k++)
      {
        temp = TomoInputsPtr->ErrorSinoThresh*TomoInputsPtr->ErrorSinoDelta*sqrt(TomoInputsPtr->Weight[k][i][j]);
        if (SinogramPtr->ProjSelect[k][i][j] == true)
        {
          numerator[i][j] += ErrorSino[k][i][j]*TomoInputsPtr->Weight[k][i][j];
          denominator[i][j] += TomoInputsPtr->Weight[k][i][j];
        }
        else
        {
	  sign = (ErrorSino[k][i][j] > 0) - (ErrorSino[k][i][j] < 0);
          numerator[i][j] += temp*sign;
          denominator[i][j] += temp/fabs(ErrorSino[k][i][j]);
        }
      }
    }

    #pragma omp parallel for private(j, k)
    for (i = 0; i < SinogramPtr->N_r; i++)
    for (j = 0; j < SinogramPtr->N_t; j++)
    {
      SinogramPtr->ProjOffset[i][j] = SinogramPtr->ProjOffset[i][j] + (numerator[i][j])/denominator[i][j];
      for (k = 0; k < SinogramPtr->N_p; k++)
      {
	ErrorSino[k][i][j] -= (numerator[i][j])/denominator[i][j];
        if (fabs(ErrorSino[k][i][j]*sqrt(TomoInputsPtr->Weight[k][i][j])) < TomoInputsPtr->ErrorSinoThresh)
        SinogramPtr->ProjSelect[k][i][j] = true;
        else
        SinogramPtr->ProjSelect[k][i][j] = false;
      }
    }
    
   multifree(numerator,2);
   multifree(denominator,2);
  }

void update_Sinogram_Offset (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr, Real_arr_t*** ErrorSino)
{
	if (TomoInputsPtr->OffsetConstraintType == 1)	
		update_d_offset_unconstrained (SinogramPtr, TomoInputsPtr, ErrorSino);
	else if (TomoInputsPtr->OffsetConstraintType == 2)	
		update_d_offset_zero_mean_constraint (SinogramPtr, TomoInputsPtr, ErrorSino);
	else if (TomoInputsPtr->OffsetConstraintType == 3)
		update_d_offset_rect_patch_constraint (SinogramPtr, TomoInputsPtr, ErrorSino);
}

  /*Implements mutithreaded shared memory parallelization using OpenMP and splits work among
  threads. Each thread gets a certain time slice and z block to update.
  Multithreading is done within the z-blocks assigned to each node.
  ErrorSino - Error sinogram
  Iter - Present iteration number
  MagUpdateMap - Magnitude update map containing the magnitude of update of each voxel
  Mask - If a certain element is true then the corresponding voxel is updated*/
  int updateVoxelsTimeSlices(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr, Real_arr_t** DetectorResponse, /*AMatrixCol* VoxelLineResponse,*/ Real_arr_t*** ErrorSino, int32_t Iter, Real_arr_t**** MagUpdateMap, uint8_t*** Mask)
  {
    Real_t AverageUpdate = 0, tempUpdate, avg_update_percentage, total_vox_mag = 0.0, vox_mag = 0.0;
    int32_t xy_start, xy_end, i, j, K, block, idx, **z_start, **z_stop;
    Real_t tempTotPix = 0, total_pix = 0;
    long int **zero_count, total_zero_count = 0;
    int32_t** thread_num = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
    MPI_Request *send_reqs, *recv_reqs;
    send_reqs = (MPI_Request*)get_spc(ScannedObjectPtr->N_time, sizeof(MPI_Request));
    recv_reqs = (MPI_Request*)get_spc(ScannedObjectPtr->N_time, sizeof(MPI_Request));
    z_start = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
    z_stop = (int32_t**)multialloc(sizeof(int32_t), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
    
    randomly_select_x_y (ScannedObjectPtr, TomoInputsPtr, Mask);
    
    zero_count = (long int**)multialloc(sizeof(long int), 2, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks);
    /*	offset_numerator = (Real_t**)multialloc(sizeof(Real_t), 2, SinogramPtr->N_r, SinogramPtr->N_t);
    memset(&(offset_denominator[0][0]), 0, SinogramPtr->N_r*SinogramPtr->N_t*sizeof(Real_t));
    
    for (k = 0; k < SinogramPtr->N_p; k++)
    for (i = 0; i < SinogramPtr->N_r; i++)
    for (j = 0; j < SinogramPtr->N_t; j++)
    offset_denominator[i][j] += TomoInputsPtr->Weight[k][i][j]; */
    
    memset(&(zero_count[0][0]), 0, ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*sizeof(long int));
    /*	K = ScannedObjectPtr->N_time*ScannedObjectPtr->N_z*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
    K = (K - total_zero_count)/(ScannedObjectPtr->gamma*K);*/
    K = ScannedObjectPtr->NHICD_Iterations;
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Number of NHICD iterations is %d.\n", K);
    for (j = 0; j < K; j++)
    {
      total_vox_mag = 0.0;
      #pragma omp parallel for collapse(2) private(i, block, idx, xy_start, xy_end) reduction(+:total_vox_mag)
      for (i = 0; i < ScannedObjectPtr->N_time; i++)
      for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
      {
        idx = (i % 2 == 0) ? block: block + 1;
        z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);
        z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;
        z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx];
        xy_start = j*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K);
        xy_end = (j + 1)*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K) - 1;
        xy_end = (j == K - 1) ? TomoInputsPtr->UpdateSelectNum[i][idx] - 1: xy_end;
        /*	printf ("Loop 1 Start - j = %d, i = %d, idx = %d, z_start = %d, z_stop = %d, xy_start = %d, xy_end = %d\n", j, i, idx, z_start[i][idx], z_stop[i][idx], xy_start, xy_end);*/
        total_vox_mag += updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], xy_start, xy_end, TomoInputsPtr->x_rand_select[i][idx], TomoInputsPtr->y_rand_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, /*VoxelLineResponse,*/ Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
        thread_num[i][idx] = omp_get_thread_num();
      }
      
      /*check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Send MPI info\n");*/
      MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
      /*	check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "update_Sinogram_Offset: Will compute projection offset error\n");*/
      if (TomoInputsPtr->updateProjOffset > 1)
      update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
      /*	check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "update_Sinogram_Offset: Done computing projection offset error\n");*/
      MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
      #pragma omp parallel for collapse(2) private(i, block, idx, xy_start, xy_end) reduction(+:total_vox_mag)
      for (i = 0; i < ScannedObjectPtr->N_time; i++)
      for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
      {
        idx = (i % 2 == 0) ? block + 1: block;
        z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);
        z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;
        z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx];
        xy_start = j*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K);
        xy_end = (j + 1)*floor(TomoInputsPtr->UpdateSelectNum[i][idx]/K) - 1;
        xy_end = (j == K - 1) ? TomoInputsPtr->UpdateSelectNum[i][idx] - 1: xy_end;
        total_vox_mag += updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], xy_start, xy_end, TomoInputsPtr->x_rand_select[i][idx], TomoInputsPtr->y_rand_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, /*VoxelLineResponse,*/ Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
        thread_num[i][idx] = omp_get_thread_num();
        /*	printf ("Loop 2 - i = %d, idx = %d, z_start = %d, z_stop = %d, xy_start = %d, xy_end = %d\n", i, idx, z_start[i][idx], z_stop[i][idx], xy_start, xy_end);*/
      }
      
      MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
      if (TomoInputsPtr->updateProjOffset > 1)
      update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
      MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
      VSC_based_Voxel_Line_Select(ScannedObjectPtr, TomoInputsPtr, MagUpdateMap);
      /*	check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Number of NHICD voxel lines to be updated in iteration %d is %d\n", j, num_voxel_lines);*/
      if (Iter > 1 && TomoInputsPtr->no_NHICD == 0)
      {
        #pragma omp parallel for collapse(2) private(i, block, idx)
        for (i = 0; i < ScannedObjectPtr->N_time; i++)
        for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
        {
          idx = (i % 2 == 0) ? block: block + 1;
          z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);
          z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;
          z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx];
          updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], 0, TomoInputsPtr->NHICDSelectNum[i][idx]-1, TomoInputsPtr->x_NHICD_select[i][idx], TomoInputsPtr->y_NHICD_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, /*VoxelLineResponse,*/ Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
          thread_num[i][idx] = omp_get_thread_num();
          /*	printf ("Loop 1 NHICD - i = %d, idx = %d, z_start = %d, z_stop = %d\n", i, idx, z_start[i][idx], z_stop[i][idx]);*/
        }
        
        MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
        if (TomoInputsPtr->updateProjOffset > 1)
        update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
        MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 0);
        
        #pragma omp parallel for collapse(2) private(i, block, idx)
        for (i = 0; i < ScannedObjectPtr->N_time; i++)
        for (block = 0; block < TomoInputsPtr->num_z_blocks; block = block + 2)
        {
          idx = (i % 2 == 0) ? block + 1: block;
          z_start[i][idx] = idx*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);
          z_stop[i][idx] = (idx + 1)*floor(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks) - 1;
          z_stop[i][idx] = (idx >= TomoInputsPtr->num_z_blocks - 1) ? ScannedObjectPtr->N_z - 1: z_stop[i][idx];
          updateVoxels (i, i, z_start[i][idx], z_stop[i][idx], 0, TomoInputsPtr->NHICDSelectNum[i][idx]-1, TomoInputsPtr->x_NHICD_select[i][idx], TomoInputsPtr->y_NHICD_select[i][idx], SinogramPtr, ScannedObjectPtr, TomoInputsPtr, ErrorSino, DetectorResponse, /*VoxelLineResponse,*/ Iter, &(zero_count[i][idx]), MagUpdateMap[i][idx], Mask[i]);
          thread_num[i][idx] = omp_get_thread_num();
          /*	printf ("Loop 2 NHICD - i = %d, idx = %d, z_start = %d, z_stop = %d\n", i, idx, z_start[i][idx], z_stop[i][idx]);*/
        }
        
        MPI_Send_Recv_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
        if (TomoInputsPtr->updateProjOffset > 1)
        update_Sinogram_Offset (SinogramPtr, TomoInputsPtr, ErrorSino);
        MPI_Wait_Z_Slices (ScannedObjectPtr, TomoInputsPtr, send_reqs, recv_reqs, 1);
      }
    }
    
    if (TomoInputsPtr->updateVar == 1)
    update_variance_parameter (SinogramPtr, TomoInputsPtr, ErrorSino);
    
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Time Slice, Z Start, Z End - Thread : ");
    total_pix = 0;
    for (i=0; i<ScannedObjectPtr->N_time; i++){
      for (block=0; block<TomoInputsPtr->num_z_blocks; block++){
        total_pix += TomoInputsPtr->UpdateSelectNum[i][block]*(ScannedObjectPtr->N_z/TomoInputsPtr->num_z_blocks);
        for (j=0; j<TomoInputsPtr->UpdateSelectNum[i][block]; j++){
          AverageUpdate += MagUpdateMap[i][block][TomoInputsPtr->y_rand_select[i][block][j]][TomoInputsPtr->x_rand_select[i][block][j]];
        }
        total_zero_count += zero_count[i][block];
        check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "%d,%d,%d-%d; ", i, z_start[i][block], z_stop[i][block], thread_num[i][block]);
      }
    }
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "\n");
    
    MPI_Allreduce(&AverageUpdate, &tempUpdate, 1, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_pix, &tempTotPix, 1, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_vox_mag, &vox_mag, 1, MPI_REAL_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
    AverageUpdate = tempUpdate/(tempTotPix);
    AverageUpdate = convert2Hounsfield(AverageUpdate);
    check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Average voxel update over all voxels is %f, total voxels is %f.\n", AverageUpdate, tempTotPix);
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Zero count is %ld.\n", total_zero_count);
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Variance parameter divisor is %f.\n", (Real_t)TomoInputsPtr->node_num*(Real_t)SinogramPtr->N_p*(Real_t)SinogramPtr->N_r*(Real_t)SinogramPtr->N_t);
    
    multifree(zero_count,2);
    multifree(thread_num,2);
    multifree(z_start,2);
    multifree(z_stop,2);
    free(send_reqs);
    free(recv_reqs);
    /*	multifree(offset_numerator,2);
    multifree(offset_denominator,2);*/
    avg_update_percentage = 100*tempUpdate/vox_mag;
    check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Percentage average magnitude of voxel updates is %f.\n", avg_update_percentage);
    
    if (avg_update_percentage < TomoInputsPtr->StopThreshold)
    {
      check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Percentage average magnitude of voxel updates is less than convergence threshold.\n");
      return (1);
    }
    return(0);
  }

  /*ICD_BackProject calls the ICD optimization function repeatedly till the stopping criteria is met.*/
  int ICD_BackProject(Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr)
  {
    #ifndef NO_COST_CALCULATE
    Real_t cost, cost_0_iter, cost_last_iter, percentage_change_in_cost = 0;
    char costfile[100]=COST_FILENAME;
    #endif
    Real_arr_t ***ErrorSino, **H_r, *H_t;
    Real_t x, y;
    int32_t j, flag = 0, Iter, i, k;
    int dimTiff[4];
    char VarEstFile[100] = VAR_PARAM_FILENAME;
    char scaled_error_file[100] = SCALED_ERROR_SINO_FILENAME;
    time_t start;
    char detect_file[100] = DETECTOR_RESPONSE_FILENAME;
    char projselect_file[100] = PROJ_SELECT_FILENAME;
    char MagUpdateMapFile[100] = MAG_UPDATE_FILENAME;
    uint8_t ***Mask;
    /*AMatrixCol *VoxelLineResponse;*/
    #ifdef POSITIVITY_CONSTRAINT
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Enforcing positivity constraint\n");
    #endif
    
    Real_arr_t**** MagUpdateMap = (Real_arr_t****)multialloc(sizeof(Real_arr_t), 4, ScannedObjectPtr->N_time, TomoInputsPtr->num_z_blocks, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);
    H_r = (Real_arr_t **)multialloc(sizeof(Real_arr_t), 2, SinogramPtr->N_p, DETECTOR_RESPONSE_BINS + 1);
    H_t = (Real_arr_t *)get_spc(DETECTOR_RESPONSE_BINS + 1, sizeof(Real_arr_t));
    ErrorSino = (Real_arr_t***)multialloc(sizeof(Real_arr_t), 3, SinogramPtr->N_p, SinogramPtr->N_r, SinogramPtr->N_t);
    Mask = (uint8_t***)multialloc(sizeof(uint8_t), 3, ScannedObjectPtr->N_time, ScannedObjectPtr->N_y, ScannedObjectPtr->N_x);
    
    memset(&(MagUpdateMap[0][0][0][0]), 0, ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x*sizeof(Real_arr_t));
/*    omp_set_num_threads(TomoInputsPtr->num_threads);*/
    check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Number of CPU cores is %d\n", (int)omp_get_num_procs());
    /*	check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "ICD_BackProject: Number of threads is %d\n", TomoInputsPtr->num_threads) ;*/
    for (i = 0; i < ScannedObjectPtr->N_time; i++)
    for (j = 0; j < ScannedObjectPtr->N_y; j++)
    for (k = 0; k < ScannedObjectPtr->N_x; k++){
      x = ScannedObjectPtr->x0 + ((Real_t)k + 0.5)*ScannedObjectPtr->delta_xy;
      y = ScannedObjectPtr->y0 + ((Real_t)j + 0.5)*ScannedObjectPtr->delta_xy;
      if (x*x + y*y < TomoInputsPtr->radius_obj*TomoInputsPtr->radius_obj)
        Mask[i][j][k] = 1;
      else
        Mask[i][j][k] = 0;
    }
    
    DetectorResponseProfile (H_r, H_t, SinogramPtr, ScannedObjectPtr, TomoInputsPtr);
    dimTiff[0] = 1; dimTiff[1] = 1; dimTiff[2] = SinogramPtr->N_p; dimTiff[3] = DETECTOR_RESPONSE_BINS+1;
    sprintf(detect_file, "%s_n%d", detect_file, TomoInputsPtr->node_rank);
    if (TomoInputsPtr->Write2Tiff == 1)
    	if (WriteMultiDimArray2Tiff (detect_file, dimTiff, 0, 1, 2, 3, &(H_r[0][0]), 0, TomoInputsPtr->debug_file_ptr)) goto error;
    start = time(NULL);
    if (initObject(SinogramPtr, ScannedObjectPtr, TomoInputsPtr, MagUpdateMap)) goto error;
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Time taken to read object = %fmins\n", difftime(time(NULL),start)/60.0);
    if (initErrorSinogam(SinogramPtr, ScannedObjectPtr, TomoInputsPtr, H_r, ErrorSino/*, VoxelLineResponse*/)) goto error;
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Time taken to initialize object and compute error sinogram = %fmins\n", difftime(time(NULL),start)/60.0);
    #ifndef NO_COST_CALCULATE
    cost = computeCost(SinogramPtr,ScannedObjectPtr,TomoInputsPtr,ErrorSino);
    cost_0_iter = cost;
    cost_last_iter = cost;
    check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "------------- Iteration 0, Cost = %f------------\n",cost);
    if (TomoInputsPtr->node_rank == 0)
    	Write2Bin (costfile, 1, 1, 1, 1, sizeof(Real_t), &cost, TomoInputsPtr->debug_file_ptr);
    #endif /*Cost calculation endif*/
   
    start=time(NULL);
    for (Iter = 1; Iter <= TomoInputsPtr->NumIter; Iter++)
    {
      flag = updateVoxelsTimeSlices (SinogramPtr, ScannedObjectPtr, TomoInputsPtr, H_r, /*VoxelLineResponse,*/ ErrorSino, Iter, MagUpdateMap, Mask);
      if (TomoInputsPtr->WritePerIter == 1)
      	if (write_ObjectProjOff2TiffBinPerIter (SinogramPtr, ScannedObjectPtr, TomoInputsPtr)) goto error;
      #ifndef NO_COST_CALCULATE
      cost = computeCost(SinogramPtr,ScannedObjectPtr,TomoInputsPtr,ErrorSino);
      percentage_change_in_cost = ((cost - cost_last_iter)/(cost - cost_0_iter))*100.0;
      check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Percentage change in cost is %f.\n", percentage_change_in_cost);
      check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Variance parameter estimate = %f.\n", TomoInputsPtr->var_est);
      check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "------------- Iteration = %d, Cost = %f, Time since start of ICD = %fmins ------------\n",Iter,cost,difftime(time(NULL),start)/60.0);
      if (TomoInputsPtr->node_rank == 0)
	Append2Bin (costfile, 1, 1, 1, 1, sizeof(Real_t), &cost, TomoInputsPtr->debug_file_ptr);
      check_error(cost > cost_last_iter, TomoInputsPtr->node_rank == 0, TomoInputsPtr->debug_file_ptr, "Cost value increased.\n");
      cost_last_iter = cost;
      /*if (percentage_change_in_cost < TomoInputsPtr->cost_thresh && flag != 0 && Iter > 1){*/
      if (flag != 0 && Iter > 1){
        check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Convergence criteria is met.\n");
        break;
      }
      #else
      check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Variance parameter estimate = %f\n",TomoInputsPtr->var_est);
      check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "-------------ICD_BackProject: ICD Iter = %d, time since start of ICD = %fmins------------.\n",Iter,difftime(time(NULL),start)/60.0);
      if (flag != 0 && Iter > 1){
        check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Convergence criteria is met.\n");
        break;
      }
      #endif
      flag = fflush(TomoInputsPtr->debug_file_ptr);
      if (flag != 0)
      	check_warn(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Cannot flush buffer.\n");
    }
   
    for (i = 0; i < SinogramPtr->N_p; i++)
    for (j = 0; j < SinogramPtr->N_r; j++)
    for (k = 0; k < SinogramPtr->N_t; k++)
    	ErrorSino[i][j][k] *= sqrt(TomoInputsPtr->Weight[i][j][k]);
    
    if (TomoInputsPtr->node_rank == 0)
    	Write2Bin (VarEstFile, 1, 1, 1, 1, sizeof(Real_t), &(TomoInputsPtr->var_est), TomoInputsPtr->debug_file_ptr);
    int32_t size = ScannedObjectPtr->N_time*TomoInputsPtr->num_z_blocks*ScannedObjectPtr->N_y*ScannedObjectPtr->N_x;
    if (write_SharedBinFile_At (MagUpdateMapFile, &(MagUpdateMap[0][0][0][0]), TomoInputsPtr->node_rank*size, size, TomoInputsPtr->debug_file_ptr)) goto error;
    sprintf(scaled_error_file, "%s_n%d", scaled_error_file, TomoInputsPtr->node_rank);
    sprintf(projselect_file, "%s_n%d", projselect_file, TomoInputsPtr->node_rank);
    dimTiff[0] = 1; dimTiff[1] = SinogramPtr->N_p; dimTiff[2] = SinogramPtr->N_r; dimTiff[3] = SinogramPtr->N_t;
    if (TomoInputsPtr->Write2Tiff == 1)
    {
      if (WriteMultiDimArray2Tiff (scaled_error_file, dimTiff, 0, 3, 1, 2, &(ErrorSino[0][0][0]), 0, TomoInputsPtr->debug_file_ptr)) goto error;
      if (WriteBoolArray2Tiff (projselect_file, dimTiff, 0, 3, 1, 2, &(SinogramPtr->ProjSelect[0][0][0]), 0, TomoInputsPtr->debug_file_ptr)) goto error;
    }
    
    multifree(ErrorSino,3);
    multifree(H_r,2);
    free(H_t);
    multifree(Mask,3);
    multifree(MagUpdateMap, 4);
    
    check_debug(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Finished running ICD_BackProject.\n");
    flag = fflush(TomoInputsPtr->debug_file_ptr);
    if (flag != 0 )
       check_warn(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "Cannot flush buffer.\n");
    
    check_info(TomoInputsPtr->node_rank==0, TomoInputsPtr->debug_file_ptr, "The estimated value of variance parameter is %f.\n", TomoInputsPtr->var_est);
    return(0);

error:
    multifree(ErrorSino,3);
    multifree(H_r,2);
    free(H_t);
    multifree(Mask,3);
    multifree(MagUpdateMap, 4);
    return(-1);
  }
