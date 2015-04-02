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



#ifndef XT_CONSTANTS_H
#define XT_CONSTANTS_H

/*#define NO_COST_CALCULATE*/
/*#define EXTRA_DEBUG_MESSAGES*/

typedef double Real_t;
typedef float Real_arr_t; /*Don't change to 'double' without first changing the floats to doubles in XT_Engine.c*/
#define MPI_REAL_DATATYPE MPI_DOUBLE
#define MPI_REAL_ARR_DATATYPE MPI_FLOAT

#define ZERO_SKIPPING
#define PROJECTIONS_FILENAME "projections"
#define WEIGHTS_FILENAME "weights"
#define OBJECT_FILENAME "object"
#define INIT_OBJECT_FILENAME "init_object"
#define PROJ_OFFSET_FILENAME "proj_offset"
#define MAG_UPDATE_FILENAME "mag_update"
#define RUN_STATUS_FILENAME "status"
#define VAR_PARAM_FILENAME "variance_estimate"
#define SCALED_ERROR_SINO_FILENAME "scaled_errorsino"
#define DETECTOR_RESPONSE_FILENAME "detector_response"
#define PROJ_SELECT_FILENAME "proj_select"
#define COST_FILENAME "cost"

#define OBJECT_INIT_VAL 0
#define MRF_P 1.2
#define MRF_Q 2.0

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923132169163975144   /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4         0.785398163397448309615660845819875721  /* pi/4 */
#endif

#define PROFILE_RESOLUTION 1536
#define BEAM_RESOLUTION 512
#define DETECTOR_RESPONSE_BINS 64

#define NHOOD_Z_MAXDIM 3
#define NHOOD_Y_MAXDIM 3
#define NHOOD_X_MAXDIM 3
#define NHOOD_TIME_MAXDIM 3

#define MAX_NUM_ITERATIONS 1000
#define OVER_RELAXATION_FACTOR 1.5
#define ENABLE_TIFF_WRITES 1 /*To disable generating tiff images use '0'*/
#define COST_CONVG_THRESHOLD 0.1
#define PROJ_OFFSET_INIT 0
#define NO_NHICD 0
#define WRITE_EVERY_ITER 0
#define ZINGER_ENABLE_PARAM_T 4.0
#define ZINGER_ENABLE_PARAM_DELTA 0.2
#define VAR_PARAM_INIT 1
#define COMPUTE_RMSE_CONVG 0
#define ZINGER_DISABLE_PARAM_T 100000
#define ZINGER_DISABLE_PARAM_DELTA 1
#define MIN_XY_RECON_RES 64
#define MAX_MULTRES_NUM 8
#define MIN_ROWS_PER_NODE 2
#define MIN_PROJECTION_ROWS 4

#define HFIELD_UNIT_CONV_CONST 0.0001
#define AIR_MASS_ATT_COEFF 0.496372335005353 /*in cm^2/g. computed using cubic interpolation*/
#define WATER_MASS_ATT_COEFF 0.521225397034623

#define WATER_DENSITY 1.0 /*in g/cm^3*/ 
#define AIR_DENSITY 0.001205
#define HOUNSFIELD_WATER_MAP 1000
#define HOUNSFIELD_AIR_MAP 0
#ifndef HOUNSFIELD_MAX
	#define HOUNSFIELD_MAX 60000
#endif
#ifndef HOUNSFIELD_MIN
	#define HOUNSFIELD_MIN 10000
#endif

#define PHANTOM_XY_SIZE 1024
#define PHANTOM_Z_SIZE 4

#endif /*#ifndef XT_CONSTANTS_H*/
