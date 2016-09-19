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
#include <stdint.h>
#include <mpi.h>
#include <getopt.h>
#include <stdlib.h>
#include "mbir4d.h"
#include <math.h>
#include "XT_Views.h"
#include "XT_HDFIO.h"
#include "XT_Constants.h"
#include <string.h>
/*Function prototype definitions which will be defined later in the file.*/
void read_command_line_args (int32_t argc, char **argv, char path2data[], char path2whites[], char path2darks[], int32_t *proj_rows, int32_t *datafile_row0, int32_t *proj_cols, int32_t *proj_start, int32_t *proj_num, int32_t *K, int32_t *N_theta, int32_t *r, float *min_acquire_time, float *rotation_speed, float *vox_wid, float *rot_center, float *sig_s, float *sig_t, float *c_s, float *c_t, float *convg_thresh, int32_t *remove_rings, int32_t *quad_convex, float *huber_delta, float *huber_T, uint8_t *restart, FILE* debug_msg_ptr);

/*The main function which reads the command line arguments, reads the data,
  and does the reconstruction.*/
int main(int argc, char **argv)
{
	uint8_t restart;
	char path2data[10000], path2whites[10000], path2darks[10000];
	int32_t proj_rows, datafile_row0, proj_cols, proj_start, proj_num, K, N_theta, r, remove_rings, quad_convex, nodes_num, nodes_rank, total_projs;
	float *object, *projections, *weights, *proj_angles, *proj_times, *recon_times, min_acq_time, rot_speed, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, huber_delta, huber_T;
	FILE *debug_msg_ptr;

	/*initialize MPI process.*/	
	MPI_Init(&argc, &argv);
	/*Find the total number of nodes.*/
	MPI_Comm_size(MPI_COMM_WORLD, &nodes_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &nodes_rank);
	
	/*All messages to help debug any potential mistakes or bugs are written to debug.log*/
/*	debug_msg_ptr = fopen("debug.log", "w");*/
	debug_msg_ptr = stdout;
	/*Read the command line arguments to determine the reconstruction parameters*/
	read_command_line_args (argc, argv, path2data, path2whites, path2darks, &proj_rows, &datafile_row0, &proj_cols, &proj_start, &proj_num, &K, &N_theta, &r, &min_acq_time, &rot_speed, &vox_wid, &rot_center, &sig_s, &sig_t, &c_s, &c_t, &convg_thresh, &remove_rings, &quad_convex, &huber_delta, &huber_T, &restart, debug_msg_ptr);
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Number of nodes is %d and command line input argument values are - path2data = %s, path2whites = %s, path2darks = %s, proj_rows = %d, datafile_row0 = %d, proj_cols = %d, proj_start = %d, proj_num = %d, K = %d, N_theta = %d, r = %d, min_acq_time = %f, rot_speed = %f, vox_wid = %f, rot_center = %f, sig_s = %f, sig_t = %f, c_s = %f, c_t = %f, convg_thresh = %f, remove_rings = %d, quad_convex = %d, huber_delta = %f, huber_T = %f, restart = %d\n", nodes_num, path2data, path2whites, path2darks, proj_rows, datafile_row0, proj_cols, proj_start, proj_num, K, N_theta, r, min_acq_time, rot_speed, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, quad_convex, huber_delta, huber_T, restart);	

	rot_speed = rot_speed*M_PI/180.0;
	/*Allocate memory for data arrays used for reconstruction.*/
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Allocating memory for data ....\n");
	projections = (float*)calloc ((proj_num*proj_rows*proj_cols)/nodes_num, sizeof(float));
	weights = (float*)calloc ((proj_num*proj_rows*proj_cols)/nodes_num, sizeof(float));
	proj_angles = (float*)calloc (proj_num, sizeof(float));
	proj_times = (float*)calloc (proj_num, sizeof(float));

        if (projections == NULL || weights == NULL || proj_angles == NULL || proj_times == NULL)
                fprintf(debug_msg_ptr, "ERROR: In function main(), calloc() returned NULL!\n");

	/*Read data*/
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reading data ....\n");
	if (read_data (path2data, path2whites, path2darks, projections, weights, datafile_row0, proj_rows, proj_cols, proj_start, proj_num, debug_msg_ptr)) 
	{
		fprintf(debug_msg_ptr, "ERROR: Cannot read data files\n");
		return (-1);
	}
	total_projs = gen_interlaced_views_0_to_Inf(proj_angles, proj_times, K, N_theta, proj_start, proj_num, min_acq_time, rot_speed, debug_msg_ptr);
	fprintf(debug_msg_ptr, "The effective number of view angles after deletion is %d\n", total_projs);
	
	float step = N_theta/r;
	int32_t i, recon_num = floor(((float)(proj_times[proj_num-1] - proj_times[0]))/step + 0.5);
	fprintf (debug_msg_ptr, "Number of reconstruction time steps are %d.\n", recon_num);
	recon_times = (float*)calloc (recon_num + 1, sizeof(float));
	recon_times[0] = proj_times[0];
	for (i = 1; i <= recon_num; i++)
		recon_times[i] = recon_times[i-1] + step;
		
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reconstructing the data ....\n");
	/*Run the reconstruction*/
	reconstruct (&object, projections, weights, proj_angles, proj_times, recon_times, proj_rows, proj_cols, proj_num, recon_num, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, quad_convex, huber_delta, huber_T, restart, debug_msg_ptr);
/*	free(object);
	free(projections);
	free(weights);
	free(proj_angles);
	free(proj_times);
	free(recon_times);*/

	fclose (debug_msg_ptr); 
	MPI_Finalize();
	return (0);
}


/*Function which parses the command line input to the C code and initializes several variables.*/
void read_command_line_args (int32_t argc, char **argv, char path2data[], char path2whites[], char path2darks[], int32_t *proj_rows, int32_t *datafile_row0, int32_t *proj_cols, int32_t *proj_start, int32_t *proj_num, int32_t *K, int32_t *N_theta, int32_t *r, float *min_acquire_time, float *rotation_speed, float *vox_wid, float *rot_center, float *sig_s, float *sig_t, float *c_s, float *c_t, float *convg_thresh, int32_t *remove_rings, int32_t *quad_convex, float *huber_delta, float *huber_T, uint8_t *restart, FILE* debug_msg_ptr)
{
	int32_t option_index;
	char c;
	static struct option long_options[] =
        {
               {"path2data",  required_argument, 0, 'a'}, /*Path to the input HDF file.*/
               {"path2whites",  required_argument, 0, 'b'}, /*Path to the input HDF file.*/
               {"path2darks",  required_argument, 0, 'c'}, /*Path to the input HDF file.*/
               {"proj_rows",  required_argument, 0, 'd'}, 
               {"datafile_row0",  required_argument, 0, 'e'}, 
/*Number of rows (or slices) in the projection image. Typically, it is the number of detector bins in the axial direction.*/
               {"proj_cols",  required_argument, 0, 'f'}, 
/*Number of columns in the projection image. Typically, it is the number of detector bins in the cross-axial direction.*/
               {"proj_start",  required_argument, 0, 'g'},
/*Starting projection used for reconstruction*/ 
               {"proj_num",  required_argument, 0, 'h'}, 
/*Total number of 2D projections used for reconstruction.*/
               {"K",  required_argument, 0, 'i'},
/*Number of interlaced sub-frames in one frame.*/ 
               {"N_theta",  required_argument, 0, 'j'}, 
/*Number of projections in one frame.*/
               {"r",  required_argument, 0, 'k'}, 
/*Number of reconstruction time steps in one frame.*/
               {"min_acquire_time",  required_argument, 0, 'l'}, 
/*The minimum time between views. If the time between views is less than 'min_acquire_time', then that the 2nd view is deleted from the list.*/
               {"rotation_speed",  required_argument, 0, 'm'}, 
/*Rotation speed of the object in units of degrees per second.*/
               {"vox_wid",  required_argument, 0, 'n'}, /*Side length of a cubic voxel in inverse units of linear attenuation coefficient of the object. 
		For example, if units of "vox_wid" is mm, then attenuation coefficient will have units of mm^-1, and vice versa.
		Note that attenuation coefficient is what we are trying to reconstruct.*/
               {"rot_center",    required_argument, 0, 'o'}, /*Center of rotation of object, in units of detector pixels. 
		For example, if center of rotation is exactly at the center of the object, then rot_center = proj_num_cols/2.
		If not, then specify as to which detector column does the center of rotation of the object projects to. */
               {"sig_s",  required_argument, 0, 'p'}, /*Spatial regularization parameter of the prior model.*/
               {"sig_t",  required_argument, 0, 'q'}, /*Temporal regularization parameter of the prior model.*/
               {"c_s",  required_argument, 0, 'r'}, 
		/*parameter of the spatial qGGMRF prior model. 
 		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_s.*/ 
               {"c_t",  required_argument, 0, 's'}, 
		/*parameter of the temporal qGGMRF prior model. 
  		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_t.*/ 
               {"convg_thresh",    required_argument, 0, 't'}, /*Used to determine when the algorithm is converged at each stage of multi-resolution.
		If the ratio of the average magnitude of voxel updates to the average voxel value expressed as a percentage is less
		than "convg_thresh" then the algorithm is assumed to have converged and the algorithm stops.*/
               {"remove_rings",    required_argument, 0, 'u'}, /*If specified, it models the detector non-uniformities which result in unknown offset error in the projections. '0' means no ring correction. '1' enables ring correction. '2' uses the improved ring correction but might introduce a mean shift in the reconstruction. 
		The ring artifacts in the reconstruction should reduce.*/
               {"quad_convex",    no_argument, 0, 'v'}, 
		/*If this flag is passed when passing the executable, then the algorithm uses a convex quadratic forward model. This model does not account for the zinger measurements which causes streak artifacts in the reconstruction. If '0', then the algorithm uses a generalized Huber function which models the effect of zingers. This reduces streak artifacts in the reconstruction. Also, using '1' disables estimation of variance parameter 'sigma' and '0' enables it.*/
		{"huber_delta", required_argument, 0, 'x'},
		/*The parameter \delta of the generalized Huber function which models the effect of zingers. Legal values are in the range 0 to 1.*/
		{"huber_T", required_argument, 0, 'y'}, 
		/*The threshold parameter T of the generalized Huber function. All positive values are legal values.*/
               {"restart",    no_argument, 0, 'w'}, /*If the reconstruction gets killed due to any unfortunate reason (like exceeding walltime in a super-computing cluster), use this flag to restart the reconstruction from the beginning of the current multi-resolution stage. Don't use restart if WRITE_EVERY_ITER  is 1.*/
               {0, 0, 0, 0}
         };

	*restart = 0;
	*quad_convex = 0;
	*huber_delta = GEN_HUBER_PARAM_DELTA;
	*huber_T = GEN_HUBER_PARAM_T;
	while(1)
	{		
	   c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:vwx:y:", long_options, &option_index);
           /* Detect the end of the options. */
          if (c == -1) break;
	  switch (c) { 
		case  0 : fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Argument not recognized\n");		break;
		case 'a': strcpy(path2data, optarg);				break;
		case 'b': strcpy(path2whites, optarg);				break;
		case 'c': strcpy(path2darks, optarg);				break;
		case 'd': *proj_rows = (int32_t)atoi(optarg);			break;
		case 'e': *datafile_row0 = (int32_t)atoi(optarg);			break;
		case 'f': *proj_cols = (int32_t)atoi(optarg);			break;
		case 'g': *proj_start = (int32_t)atoi(optarg);			break;
		case 'h': *proj_num = (int32_t)atoi(optarg);			break;
		case 'i': *K = (int32_t)atoi(optarg);			break;
		case 'j': *N_theta = (int32_t)atoi(optarg);			break;
		case 'k': *r = (int32_t)atoi(optarg);			break;
		case 'l': *min_acquire_time = (int32_t)atof(optarg);			break;
		case 'm': *rotation_speed = (int32_t)atof(optarg);			break;
		case 'n': *vox_wid = (float)atof(optarg);			break;
		case 'o': *rot_center = (float)atof(optarg);			break;
		case 'p': *sig_s = (float)atof(optarg);			break;
		case 'q': *sig_t = (float)atof(optarg);			break;
		case 'r': *c_s = (float)atof(optarg);				break;
		case 's': *c_t = (float)atof(optarg);				break;
		case 't': *convg_thresh = (float)atof(optarg);			break;
		case 'u': *remove_rings = (int32_t)atoi(optarg);		break;
		case 'v': *quad_convex = 1;		break;
		case 'x': *huber_delta = (float)atof(optarg);		break;
		case 'y': *huber_T = (float)atof(optarg);		break;
		case 'w': *restart = 1;		break;
		case '?': fprintf(debug_msg_ptr, "ERROR: Cannot recognize argument %s\n",optarg); break;
		}
	}

	if(argc-optind > 0)
		fprintf(debug_msg_ptr, "ERROR: Argument list has an error\n");
}


