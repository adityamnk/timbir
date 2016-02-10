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


#ifndef XT_IOMISC_H
#define XT_IOMISC_H

#include "XT_Structures.h"
#include <stdlib.h>
Real_arr_t convert_HU2um(Real_arr_t val);
Real_arr_t convert_um2HU (Real_arr_t val);
int32_t Write2Bin (char *filename, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t datatype_size, void* img, FILE *debug_file_ptr);
int32_t Append2Bin (char *filename, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t datatype_size, void* img, FILE *debug_file_ptr);
int32_t Read4mBin (char *filename, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t datatype_size, void* img, FILE *debug_file_ptr);
int32_t WriteMultiDimArray2Tiff (char *filename, int dim[4], int dim2loop_1, int dim2loop_2, int dim2write_1, int dim2write_2, Real_arr_t* img, int hounsfield_flag, FILE* debug_file_ptr);
int32_t write_ObjectProjOff2TiffBinPerIter (Sinogram* SinogramPtr, ScannedObject* ScannedObjectPtr, TomoInputs* TomoInputsPtr);
int32_t WriteBoolArray2Tiff (char *filename, int dim[4], int dim2loop_1, int dim2loop_2, int dim2write_1, int dim2write_2, bool* imgin, int hounsfield_flag, FILE* debug_file_ptr);
int32_t WriteInt32Tiff(char* filename, int height, int width, int32_t** imgin, int hounsfield_flag, FILE *debug_file_ptr);
int32_t Write2Tiff(char* filename, int height, int width, Real_arr_t** img, int hounsfield_flag, FILE *debug_file_ptr);
void Arr1DToArr2D (Real_arr_t* arr1d, Real_arr_t** arr2d, int32_t N_y, int32_t N_x);
void Arr2DToArr1D (Real_arr_t* arr1d, Real_arr_t** arr2d);
void Arr1DToArr3D (Real_arr_t* arr1d, Real_arr_t*** arr3d, int32_t N_z, int32_t N_y, int32_t N_x);
void Arr3DToArr1D (Real_arr_t* arr1d, Real_arr_t*** arr3d);
Real_arr_t**** Arr1DToArr4D (Real_arr_t* arr1d, int32_t N_t, int32_t N_z, int32_t N_y, int32_t N_x);
Real_arr_t* Arr4DToArr1D (Real_arr_t**** arr4d);

#endif
