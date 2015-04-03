#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <stdlib.h>
#include <getopt.h>
#include "XT_Main.h"
#include "mbir4d.h"
#include "string.h"
#include "XT_HDFIO.h"

/*Function prototype definitions which will be defined later in the file.*/
int32_t read_data (char data_filename[], float *projections, float *weights, float *proj_angles, float *proj_times, int32_t datafile_row0, int32_t proj_rows, int32_t proj_cols, int32_t proj_start, int32_t proj_num, FILE* debug_file_ptr);
void read_command_line_args (int32_t argc, char **argv, char path2data[], int32_t *datafile_row0, int32_t *proj_rows, int32_t *proj_cols, int32_t *proj_start, int32_t *proj_num, int32_t *recon_num, float *vox_wid, float *rot_center, float *sig_s, float *sig_t, float *c_s, float *c_t, float *convg_thresh, int32_t *remove_rings, int32_t *remove_streaks, uint8_t *restart, FILE* debug_msg_ptr);

/*The main function which reads the command line arguments, reads the data,
  and does the reconstruction.*/
int main(int argc, char **argv)
{
	uint8_t restart;
	int32_t proj_rows, proj_cols, proj_num, recon_num, remove_rings, remove_streaks, nodes_num, nodes_rank, datafile_row0, proj_start, i;
	float *object, *projections, *weights, *proj_angles, *proj_times, *recon_times, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, delta;
	FILE *debug_msg_ptr;
	char path2data[10000];

	/*initialize MPI process.*/	
	MPI_Init(&argc, &argv);
	/*Find the total number of nodes.*/
	MPI_Comm_size(MPI_COMM_WORLD, &nodes_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &nodes_rank);
	
	/*All messages to help debug any potential mistakes or bugs are written to debug.log*/
	debug_msg_ptr = fopen("debug.log", "w");
	debug_msg_ptr = stdout;
	
	/*Read the command line arguments to determine the reconstruction parameters*/
	read_command_line_args (argc, argv, path2data, &datafile_row0, &proj_rows, &proj_cols, &proj_start, &proj_num, &recon_num, &vox_wid, &rot_center, &sig_s, &sig_t, &c_s, &c_t, &convg_thresh, &remove_rings, &remove_streaks, &restart, debug_msg_ptr);
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Number of nodes is %d and command line input argument values are proj_rows = %d, datafile_row0 = %d, proj_cols = %d, proj_start = %d, proj_num = %d, recon_num = %d, vox_wid = %f, rot_center = %f, sig_s = %f, sig_t = %f, c_s = %f, c_t = %f, convg_thresh = %f, remove_rings = %d, remove_streaks = %d, restart = %d\n", nodes_num, proj_rows, datafile_row0, proj_cols, proj_start, proj_num, recon_num, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, remove_streaks, restart);	
	
	/*Allocate memory for data arrays used for reconstruction.*/
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Allocating memory for data ....\n");
	projections = (float*)calloc ((proj_num*proj_rows*proj_cols)/nodes_num, sizeof(float));
	weights = (float*)calloc ((proj_num*proj_rows*proj_cols)/nodes_num, sizeof(float));
	proj_angles = (float*)calloc (proj_num, sizeof(float));
	proj_times = (float*)calloc (proj_num, sizeof(float));
	recon_times = (float*)calloc (recon_num + 1, sizeof(float));

	/*Read data*/
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reading data ....\n");
	if (read_data (path2data, projections, weights, proj_angles, proj_times, datafile_row0, proj_rows, proj_cols, proj_start, proj_num, debug_msg_ptr)) {goto error;}

	recon_times[0] = proj_times[0];
	recon_times[recon_num] = proj_times[proj_num-1];
	delta = (proj_times[proj_num-1] - proj_times[0])/recon_num;
	for (i = 1; i < recon_num; i++)
		recon_times[i] = recon_times[i-1] + delta;

	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reconstructing the data ....\n");
	/*Run the reconstruction*/
	reconstruct (&object, projections, weights, proj_angles, proj_times, recon_times, proj_rows, proj_cols, proj_num, recon_num, vox_wid, rot_center, sig_s, sig_t, c_s, c_t, convg_thresh, remove_rings, remove_streaks, restart, debug_msg_ptr);
	if (nodes_rank == 0) fprintf(debug_msg_ptr, "main: Reconstruction complete!\n");
	if (object)
		free(object);
	else
		fprintf(debug_msg_ptr, "ERROR: main: Reconstruction values not returned!\n");
		
	free(projections);
	free(weights);
	free(proj_angles);
	free(proj_times);
	free(recon_times);
	fclose (debug_msg_ptr); 
	MPI_Finalize();
	return (0);

error:
	free(projections);
	free(weights);
	free(proj_angles);
	free(proj_times);
	free(recon_times);
	fclose (debug_msg_ptr); 
	MPI_Finalize();
	return (-1);
}


int32_t read_data (char data_filename[], float *projections, float *weights, float *proj_angles, float *proj_times, int32_t datafile_row0, int32_t proj_rows, int32_t proj_cols, int32_t proj_start, int32_t proj_num, FILE* debug_file_ptr)
{
	if (read_ProjWeightData (data_filename, projections, weights, datafile_row0, proj_rows, proj_cols, proj_start, proj_num, debug_file_ptr)) return (-1);
	if (read_AngleTimeReconList (data_filename, proj_angles, proj_times, proj_start, proj_num, debug_file_ptr)) return (-1);
	return (0);
}

/*Function which parses the command line input to the C code and initializes several variables.*/
void read_command_line_args (int32_t argc, char **argv, char path2data[], int32_t *datafile_row0, int32_t *proj_rows, int32_t *proj_cols, int32_t *proj_start, int32_t *proj_num, int32_t *recon_num, float *vox_wid, float *rot_center, float *sig_s, float *sig_t, float *c_s, float *c_t, float *convg_thresh, int32_t *remove_rings, int32_t *remove_streaks, uint8_t *restart, FILE* debug_msg_ptr)
{
	int32_t option_index;
	char c;
	static struct option long_options[] =
        {
               {"datafile_row0",  required_argument, 0, 'a'}, /*Starting row in projection data at which to do the reconstruction.*/
               {"proj_rows",  required_argument, 0, 'b'}, /*Number of columns in the projection image. Typically, it is the number of detector bins in the cross-axial direction.*/
               {"proj_cols",  required_argument, 0, 'c'}, /*Number of rows (or slices) in the projection image. Typically, it is the number of detector bins in the axial direction.*/
               {"proj_start",  required_argument, 0, 'd'}, /*Projection index at which to start the reconstruction.*/
               {"proj_num",  required_argument, 0, 'e'}, /*Total number of 2D projections used for reconstruction.*/
               {"recon_num",  required_argument, 0, 'f'}, /*Number of 3D time samples in the 4D reconstruction. For 3D reconstructions, this value should be set to 1.*/
               {"vox_wid",  required_argument, 0, 'g'}, /*Side length of a cubic voxel in inverse units of linear attenuation coefficient of the object. 
		For example, if units of "vox_wid" is mm, then attenuation coefficient will have units of mm^-1, and vice versa.
		Note that attenuation coefficient is what we are trying to reconstruct.*/
               {"rot_center",    required_argument, 0, 'h'}, /*Center of rotation of object, in units of detector pixels. 
		For example, if center of rotation is exactly at the center of the object, then rot_center = proj_num_cols/2.
		If not, then specify as to which detector column does the center of rotation of the object projects to. */
               {"sig_s",  required_argument, 0, 'i'}, /*Spatial regularization parameter of the prior model.*/
               {"sig_t",  required_argument, 0, 'j'}, /*Temporal regularization parameter of the prior model.*/
               {"c_s",  required_argument, 0, 'k'}, 
		/*parameter of the spatial qGGMRF prior model. 
 		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_s.*/ 
               {"c_t",  required_argument, 0, 'l'}, 
		/*parameter of the temporal qGGMRF prior model. 
  		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_t.*/ 
               {"convg_thresh",    required_argument, 0, 'm'}, /*Used to determine when the algorithm is converged at each stage of multi-resolution.
		If the ratio of the average magnitude of voxel updates to the average voxel value expressed as a percentage is less
		than "convg_thresh" then the algorithm is assumed to have converged and the algorithm stops.*/
               {"remove_rings",    required_argument, 0, 'n'}, /*If specified, it models the detector non-uniformities which corrects the ring artifacts. '0' means no ring correction. '1' enables ring correction. '2' uses the improved ring correction but might introduce a mean shift in the reconstruction. 
		The ring artifacts in the reconstruction should reduce.*/
		       {"remove_streaks",    required_argument, 0, 'o'}, /*If specified, it models the effect of anamalous measurements (also called zingers). The streak artifacts in the reconstruction should reduce.*/
		       {"restart",    no_argument, 0, 'p'}, /*If the reconstruction gets killed due to any unfortunate reason (like exceeding walltime in a super-computing cluster), use this flag to restart the reconstruction from the beginning of the current multi-resolution stage. Don't use restart if WRITE_EVERY_ITER  is 1.*/
		       {"path2data",  required_argument, 0, 'q'}, /*Path to the HDF dataset containing the projection and weight data.*/
		       {0, 0, 0, 0}
		 };

		while(1)
		{		
		   c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:", long_options, &option_index);
		   /* Detect the end of the options. */
		  if (c == -1) break;
		  switch (c) { 
			case  0 : fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Argument not recognized\n");		break;
			case 'a': *datafile_row0 = (int32_t)atoi(optarg);			break;
			case 'b': *proj_rows = (int32_t)atoi(optarg);			break;
			case 'c': *proj_cols = (int32_t)atoi(optarg);			break;
			case 'd': *proj_start = (int32_t)atoi(optarg);			break;
			case 'e': *proj_num = (int32_t)atoi(optarg);			break;
			case 'f': *recon_num = (int32_t)atoi(optarg);			break;
			case 'g': *vox_wid = (float)atof(optarg);			break;
			case 'h': *rot_center = (float)atof(optarg);			break;
			case 'i': *sig_s = (float)atof(optarg);			break;
			case 'j': *sig_t = (float)atof(optarg);			break;
			case 'k': *c_s = (float)atof(optarg);				break;
			case 'l': *c_t = (float)atof(optarg);				break;
			case 'm': *convg_thresh = (float)atof(optarg);			break;
			case 'n': *remove_rings = (int32_t)atoi(optarg);		break;
			case 'o': *remove_streaks = (int32_t)atoi(optarg);		break;
			case 'p': *restart = (uint8_t)atoi(optarg);		break;
			case 'q': strcpy(path2data, optarg);			break;
			case '?': fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Cannot recognize argument %s\n",optarg); break;
			}
		}

		if(argc-optind > 0)
			fprintf(debug_msg_ptr, "ERROR: read_command_line_args: Argument list has an error\n");
	}


