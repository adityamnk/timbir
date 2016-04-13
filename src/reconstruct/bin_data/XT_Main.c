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
#include <stdlib.h>
#include <getopt.h>
#include "XT_Main.h"
#include "mbir4d.h"

/*Function prototype definitions which will be defined later in the file.*/
void read_data (float *projections, float *weights, float *proj_angles, float *proj_times, float *recon_times, int32_t proj_rows, int32_t proj_cols, int32_t proj_num, int32_t recon_num, FILE* debug_file_ptr);
void read_command_line_args (int32_t argc, char **argv, int32_t *proj_rows, int32_t *proj_cols, int32_t *proj_num, int32_t *recon_num, float *vox_wid, float *rot_center, float *sig_s, float *sig_t, float *c_s, float *c_t, float *convg_thresh, int32_t *remove_rings, int32_t *quad_convex, float *huber_delta, float *huber_T, uint8_t *restart, FILE* debug_msg_ptr);

/*The main function which reads the command line arguments, reads the data,
  and does the reconstruction.*/
int main(int argc, char **argv)
{
	uint8_t restart;
	int32_t proj_rows, proj_cols, proj_num, recon_num, remove_rings, quad_convex, nodes_num, nodes_rank;
	float *object, *projections, *weights, *proj_angles, *proj_times, *recon_times, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, huber_delta, huber_T;
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
	read_command_line_args (argc, argv, &proj_rows, &proj_cols, &proj_num, &recon_num, &vox_wid, &rot_center, &sig_s, &sig_t, &c_s, &c_t, &convg_thresh, &remove_rings, &quad_convex, &huber_delta, &huber_T, &restart, debug_msg_ptr);
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Number of nodes is %d and command line input argument values are proj_rows = %d, proj_cols = %d, proj_num = %d, recon_num = %d, vox_wid = %f, rot_center = %f, sig_s = %f, sig_t = %f, c_s = %f, c_t = %f, convg_thresh = %f, remove_rings = %d, quad_convex = %d, huber_delta = %f, huber_T = %f, restart = %d\n", nodes_num, proj_rows, proj_cols, proj_num, recon_num, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, quad_convex, huber_delta, huber_T, restart);	
	
	/*Allocate memory for data arrays used for reconstruction.*/
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Allocating memory for data ....\n");
	projections = (float*)calloc ((proj_num*proj_rows*proj_cols)/nodes_num, sizeof(float));
	weights = (float*)calloc ((proj_num*proj_rows*proj_cols)/nodes_num, sizeof(float));
	proj_angles = (float*)calloc (proj_num, sizeof(float));
	proj_times = (float*)calloc (proj_num, sizeof(float));
	recon_times = (float*)calloc (recon_num + 1, sizeof(float));

	/*Read data*/
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reading data ....\n");
	read_data (projections, weights, proj_angles, proj_times, recon_times, proj_rows, proj_cols, proj_num, recon_num, debug_msg_ptr);
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reconstructing the data ....\n");
	/*Run the reconstruction*/
	reconstruct (&object, projections, weights, proj_angles, proj_times, recon_times, proj_rows, proj_cols, proj_num, recon_num, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, quad_convex, huber_delta, huber_T, restart, debug_msg_ptr);
	/*free(projections);
	free(weights);
	free(proj_angles);
	free(proj_times);
	free(recon_times);
	free(object);	*/

	fclose (debug_msg_ptr); 
	MPI_Finalize();
	return (0);
}

void read_BinFile (char filename[100], float* data, int32_t offset, int32_t size, FILE* debug_file_ptr)
{
	MPI_File fh;
	MPI_Status status;
	char BinFilename[100];
	int32_t len;

    	sprintf(BinFilename, "%s.bin", filename);
	MPI_File_open(MPI_COMM_WORLD, BinFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at(fh, offset*sizeof(float), data, size, MPI_FLOAT, &status);
	MPI_Get_count(&status, MPI_FLOAT, &len);
    	if(len == MPI_UNDEFINED || len != size)
	{
		fprintf (debug_file_ptr, "ERROR: read_BinFile: Read %d number of elements from the file %s.bin at an offset of %d bytes.\n. However, required number of elements is %d.", len, filename, offset, size);
		exit(1);
	}
	MPI_File_close(&fh);
}

void read_data (float *projections, float *weights, float *proj_angles, float *proj_times, float *recon_times, int32_t proj_rows, int32_t proj_cols, int32_t proj_num, int32_t recon_num, FILE* debug_file_ptr)
{
	char projections_filename[] = PROJECTIONS_FILENAME;
	char weights_filename[] = WEIGHTS_FILENAME;
	char proj_angles_filename[] = PROJ_ANGLES_FILENAME;
	char proj_times_filename[] = PROJ_TIMES_FILENAME;
	char recon_times_filename[] = RECON_TIMES_FILENAME;
	int32_t i, offset, size, rank, num_nodes;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
	for (i = 0; i < proj_num; i++)
	{
		size = proj_rows*proj_cols/num_nodes;
		offset = i*proj_rows*proj_cols + rank*size;
		read_BinFile (projections_filename, projections, offset, size, debug_file_ptr);
		read_BinFile (weights_filename, weights, offset, size, debug_file_ptr);
		projections = projections + size;
		weights = weights + size;
	}
	read_BinFile (proj_angles_filename, proj_angles, 0, proj_num, debug_file_ptr);
	read_BinFile (proj_times_filename, proj_times, 0, proj_num, debug_file_ptr);
	read_BinFile (recon_times_filename, recon_times, 0, recon_num + 1, debug_file_ptr);
}

/*Function which parses the command line input to the C code and initializes several variables.*/
void read_command_line_args (int32_t argc, char **argv, int32_t *proj_rows, int32_t *proj_cols, int32_t *proj_num, int32_t *recon_num, float *vox_wid, float *rot_center, float *sig_s, float *sig_t, float *c_s, float *c_t, float *convg_thresh, int32_t *remove_rings, int32_t *quad_convex, float *huber_delta, float *huber_T, uint8_t *restart, FILE* debug_msg_ptr)
{
	int32_t option_index;
	char c;
	static struct option long_options[] =
        {
               {"proj_rows",  required_argument, 0, 'a'}, /*Number of rows in the projection image. Typically, it is the number of detector bins in the axial direction.*/
               {"proj_cols",  required_argument, 0, 'b'}, /*Number of columns in the projection image. Typically, it is the number of detector bins in the cross-axial direction.*/
               {"proj_num",  required_argument, 0, 'c'}, /*Total number of 2D projections used for reconstruction.*/
               {"recon_num",  required_argument, 0, 'd'}, /*Number of 3D time samples in the 4D reconstruction. For 3D reconstructions, this value should be set to 1.*/
               {"vox_wid",  required_argument, 0, 'e'}, /*Side length of a cubic voxel in inverse units of linear attenuation coefficient of the object. 
		For example, if units of "vox_wid" is mm, then attenuation coefficient will have units of mm^-1, and vice versa.
		Note that attenuation coefficient is what we are trying to reconstruct.*/
               {"rot_center",    required_argument, 0, 'f'}, /*Center of rotation of object, in units of detector pixels. 
		For example, if center of rotation is exactly at the center of the object, then rot_center = proj_cols/2.
		If not, then specify as to which detector column does the center of rotation of the object projects to. */
               {"sig_s",  required_argument, 0, 'g'}, /*Spatial regularization parameter of the prior model.*/
               {"sig_t",  required_argument, 0, 'h'}, /*Temporal regularization parameter of the prior model.*/
               {"c_s",  required_argument, 0, 'i'}, 
		/*parameter of the spatial qGGMRF prior model. 
 		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_s.*/ 
               {"c_t",  required_argument, 0, 'j'}, 
		/*parameter of the temporal qGGMRF prior model. 
  		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_t.*/ 
               {"convg_thresh",    required_argument, 0, 'k'}, /*Used to determine when the algorithm is converged at each stage of multi-resolution.
		If the ratio of the average magnitude of voxel updates to the average voxel value expressed as a percentage is less
		than "convg_thresh" then the algorithm is assumed to have converged and the algorithm stops.*/
               {"remove_rings",    required_argument, 0, 'l'}, /*If specified, it models the detector non-uniformities which result in unknown offset error in the projections. '0' means no ring correction. '1' enables ring correction. '2' uses the improved ring correction but might introduce a mean shift in the reconstruction. 
		The ring artifacts in the reconstruction should reduce.*/
               {"quad_convex",    no_argument, 0, 'm'}, 
		/*Legal values are '0' and '1'. If '1', then the algorithm uses a convex quadratic forward model. This model does not account for the zinger measurements which causes streak artifacts in the reconstruction. If '0', then the algorithm uses a generalized Huber function which models the effect of zingers. This reduces streak artifacts in the reconstruction. Also, using '1' disables estimation of variance parameter 'sigma' and '0' enables it.*/
               {"huber_delta", optional_argument, 0, 'o'},
		/*The parameter \delta of the generalized Huber function which models the effect of zingers. Legal values are in the range 0 to 1.*/
		{"huber_T", optional_argument, 0, 'p'},
		/*The threshold parameter T of the generalized Huber function. All positive values are legal values.*/
		{"restart",    no_argument, 0, 'n'}, /*If the reconstruction gets killed due to any unfortunate reason (like exceeding walltime in a super-computing cluster), use this flag to restart the reconstruction from the beginning of the current multi-resolution stage. Don't use restart if WRITE_EVERY_ITER  is 1.*/
               {0, 0, 0, 0}
         };

	*restart = 0;
	*quad_convex = 0;
	*huber_delta = GEN_HUBER_PARAM_DELTA;
	*huber_T = GEN_HUBER_PARAM_T;
	while(1)
	{		
	   c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:mno::p::", long_options, &option_index);
           /* Detect the end of the options. */
          if (c == -1) break;
	  switch (c) { 
		case  0 : fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Argument not recognized\n");		break;
		case 'a': *proj_rows = (int32_t)atoi(optarg);			break;
		case 'b': *proj_cols = (int32_t)atoi(optarg);			break;
		case 'c': *proj_num = (int32_t)atoi(optarg);			break;
		case 'd': *recon_num = (int32_t)atoi(optarg);			break;
		case 'e': *vox_wid = (float)atof(optarg);			break;
		case 'f': *rot_center = (float)atof(optarg);			break;
		case 'g': *sig_s = (float)atof(optarg);			break;
		case 'h': *sig_t = (float)atof(optarg);			break;
		case 'i': *c_s = (float)atof(optarg);				break;
		case 'j': *c_t = (float)atof(optarg);				break;
		case 'k': *convg_thresh = (float)atof(optarg);			break;
		case 'l': *remove_rings = (int32_t)atoi(optarg);		break;
		case 'm': *quad_convex = 1;		break;
		case 'o': if(optarg) *huber_delta = (float)atof(optarg);		break;
		case 'p': if(optarg) *huber_T = (float)atof(optarg);		break;
		case 'n': *restart = 1;		break;
		case '?': fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Cannot recognize argument %s\n",optarg); break;
		}
	}

	if(argc-optind > 0)
		fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Argument list has an error\n");
}


