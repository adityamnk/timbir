TO DO:
- Make desired changes to the Makefile
	- Change the compiler. Default is 'mpicc'. 
	- Set code optimization levels using -ON flag, where N is the optimization level
	- Make more changes if necessary
- Make sure the MPI compiler and the openmp libraries are available
- Run 'make' in the command line
- Make generates the executable named XT_Main

	
INPUT ARGUMENTS : 
	datafile_row0 - Starting row in the projection data at which the reconstruction is done.
	proj_rows - Number of rows (or slices) of the projection used for reconstruction. It is the number of detector bins in the axial direction (i.e., axis of rotation).
	proj_cols - Total number of columns in the projection image. Typically, it is the number of detector bins in the cross-axial direction (i.e., perpendicular to axis of rotation).
	proj_start - Projection index at which to start the reconstruction.
	proj_num - Total number of 2D projections used for reconstruction.
	K - Number of interlaced sub-frames in the data.
	N_theta - Total number of projections in a frame. 
	r - Number of reconstruction time samples in one frame.
	min_acquire_time - The minimum time between views. If the time between views is less than 'min_acquire_time', then that the 2nd view is deleted from the list.
	rotation_speed - Rotation speed of the object in degrees per second.
	vox_wid - Side length of a cubic voxel in inverse units of linear attenuation coefficient of the object. 
		For example, if units of "vox_wid" is mm, then attenuation coefficient will have units of mm^-1, and vice versa.
		Note that attenuation coefficient is what we are trying to reconstruct.
	rot_center - Center of rotation of object, in units of detector pixels. 
		For example, if center of rotation is exactly at the center of the object, then rot_center = proj_cols/2.
		If not, then specify as to which detector column does the center of rotation of the object projects to.
	sig_s - Spatial regularization parameter of the qGGMRF prior model. 'sig_s' has to be varied to achieve the optimum reconstruction quality. Reducing 'sig_s' will make the reconstruction smoother and increasing it will make it sharper but also noisier.
	sig_t - Temporal regularization parameter of the qGGMRF prior model. 'sig_t' has to be varied to achieve best quality. Reducing it will increase temporal smoothness which can improve quality. However, excessive smoothing along time might introduce artifacts.
	c_s - A parameter of the spatial qGGMRF prior model. It should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sig_s. For instance, choose c_s < 0.01*D/sig_s where 'D' is a rough estimate for the maximum change in value of the reconstruction along an edge in space.
	c_t - A parameter of the temporal qGGMRF prior model. It should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_t. For instance, choose c_t < 0.01*D/sig_t where 'D' is a rough estimate for the maximum change in value of the reconstruction along a temporal edge.
	convg_thresh - Used to determine when the algorithm has converged at each stage of multi-resolution. It is expressed as a percentage (chosen in the range of 0 to 100). If the ratio of the average magnitude of voxel updates to the average voxel value expressed as a percentage is less than "convg_thresh" then the algorithm is assumed to have converged and the algorithm stops.
	remove_rings - If specified, it models the detector non-uniformities and reduces the ring artifacts in the reconstruction. '0' means no ring correction. '1' enables ring correction. '2' uses the improved ring correction by enforcing a zero mean constraint on the offset errors. '3' uses a more advanced ring correction by enforcing a zero constraint on the weighted mean of offset errors over overlapping rectangular patches. If used, the ring artifacts in the reconstruction should reduce.
	restart - If the reconstruction gets killed due to any unfortunate reason (like exceeding walltime in a super-computing cluster), use this flag to restart the reconstruction from the beginning of the last run multi-resolution stage. Don't use restart if WRITE_EVERY_ITER is 1.
	path2data - Path to the dataset in HDF format containing the measured data (photon count data). The data must be stored at /exchange/data in uint16 format.  
	path2whites - Path to the dataset in HDF format containing the measured data in the absence of sample (bright field photon count data). The data must be stored at /exchange/data_white in uint16 format. 
	path2darks - Path to the dataset in HDF format containing the measured data in the absence of the X-ray beam (dark field photon count data). The data must be stored at /exchange/data_dark in uint16 format. 
	quad_convex - If used, then the algorithm uses a convex quadratic forward model. This model does not account for the zinger measurements which causes streak artifacts in the reconstruction. If not specified, then the algorithm uses a generalized Huber function which models the effect of zingers. This reduces streak artifacts in the reconstruction. Also, using this argument disables estimation of variance parameter 'sigma'.
	huber_delta - The parameter \delta of the generalized Huber function which models the effect of zingers. Legal values are in the range 0 to 1.
	huber_T - The threshold parameter T of the generalized Huber function. All positive values are legal values.


OUTPUT FILES :
	object_time_N.bin - Binary file containing the N th time sample of the reconstruction in floating point format. The number of elements in each file is recon_num x proj_rows x proj_cols x proj_cols. 
	proj_offset.bin - Binary file containing the estimated values of the projection offset (used to correct for ring artifacts). The number of elements is proj_rows x proj_cols.
	variance_estimate.bin - Binary file with one element in floating point format. It contains the estimated value of the variance parameter. 
