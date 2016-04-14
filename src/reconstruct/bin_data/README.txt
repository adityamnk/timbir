TO DO:
- Install following dependencies
	- MPI compiler (openmpi, intel mpi, etc.)
	- Open MP 
	- Make utility
- Make desired changes to the Makefile
	- Change the compiler. Default is 'mpicc'. 
	- Set code optimization levels using -ON flag, where N is the optimization level
	- Make more changes if necessary
- Run 'make' in the command line
- Make utility generates the executable named XT_Main


INPUT FILES :
	projections.bin - Binary file containing the projection data. The projections must be in the form of a 1D floating point array in raster order (row major order) of size proj_num x proj_cols x proj_rows, where proj_num is the number of projections, proj_rows is the number of projection rows, and proj_cols is the number of projection columns. A projection is typically computed as the logarithm of the ratio of the measurement the absence of the sample of the measurement with the sample. The values are stored in row major order i.e, the fastest changing dimension is along a row, followed by columns and then the projection angles. 
	weights.bin - Binary file containing the weight data. The weights must be in the form of a 1D floating point array in raster order of size proj_num x proj_cols x proj_rows. Every entry of 'weights' is used to appropriately weigh each term of the 'projections' array in the likelihood term of MBIR. The ordering of values is same as that for projections.bin. Typically, weights is set to be equal to the raw detector measurements.
	proj_angles.bin - Binary file containing the list of angles, in radians and floating point format, at which the projections are acquired. The angular range can progressively increase, from 0 to infinity. This occurs if the object is rotated continuously, from 0 to pi radians, then to 2 x pi radians, and so on.
	proj_times.bin - Binary file containing the list of times, in floating point format, at which the projections are acquired.
	recon_times.bin - Binary file containing the list of reconstruction time steps in floating point format. The reconstruction is assumed to be peicewise constant with fixed/varying step-sizes. Thus, the binary file should contain the array with times at which the steps occur. For example, to do a single 3D reconstruction, we use a array of two elements, first element being the time of the first projection and the second element being the time of the last projection. To do a 4D reconstruction with two time samples, the binary file should contain an array with three elements with the first element being the time at which the reconstruction starts, the second element being the time at which the first reconstruction time sample ends, and the last element being the time at which the second time sample ends. 

OUTPUT FILES :
	object_time_N.bin - Binary file containing the N th time sample of the reconstruction in floating point format. The number of elements in each file is recon_num x proj_rows x proj_cols x proj_cols. 
	proj_offset.bin - Binary file containing the estimated values of the projection offset (used to correct for ring artifacts). The number of elements is proj_rows x proj_cols.
	variance_estimate.bin - Binary file with one element in floating point format. It contains the estimated value of the variance parameter. 
	
INPUT ARGUMENTS : 
	proj_rows - Number of rows (or slices) in the projection image. Typically, it is the number of detector bins in the axial direction (i.e., axis of rotation).
	proj_cols - Number of columns in the projection image. Typically, it is the number of detector bins in the cross-axial direction (i.e., perpendicular to axis of rotation).
	proj_num - Total number of 2D projections (view angles) in the binary file 'projections.bin'.
	recon_num - Total number of 3D time samples in the 4D reconstruction. For 3D reconstructions, this value should be set to 1.
	vox_wid - Side length of a cubic voxel in inverse units of linear attenuation coefficient of the object. 
		For example, if units of "vox_wid" is mm, then attenuation coefficient will have units of mm^-1, and vice versa.
		Note that attenuation coefficient is the physical quantity that is being reconstructed.
	rot_center - Center of rotation of object, in units of detector pixels. 
		For example, if center of rotation is exactly at the center of the object, then rot_center = proj_cols/2.
		If not, then specify as to which detector column does the center of rotation of the object lie.
	sig_s - Spatial regularization parameter of the qGGMRF prior model. Spatial regularization parameter to be varied to achieve the optimum reconstruction quality. Reducing 'sig_s' will make the reconstruction smoother and increasing it will make it sharper but also noisier.
	sig_t - Temporal regularization parameter of the qGGMRF prior model. Temporal regularization parameter to be varied to achieve best quality. Reducing it will increase temporal smoothness which can improve quality. However, excessive smoothing along time might introduce artifacts.
	c_s - parameter of the spatial qGGMRF prior model. 
		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sig_s. For instance, choose c_s < 0.01*D/sig_s where 'D' is a rough estimate for the maximum change in value of the reconstruction along an edge in space.
	c_t - parameter of the temporal qGGMRF prior model. 
  		Should be fixed to be much lesser (typically 0.01 times) than the ratio of voxel difference over an edge to sigma_t. For instance, choose c_t < 0.01*D/sig_t where 'D' is a rough estimate for the maximum change in value of the reconstruction along a temporal edge.
	convg_thresh - Used to determine when the algorithm is converged at each stage of multi-resolution. It is expressed as a percentage (chosen in the range of 0 to 100). If the ratio of the average magnitude of voxel updates to the average voxel value expressed as a percentage is less than "convg_thresh" then the algorithm is assumed to have converged and the algorithm stops.
	remove_rings - If specified, it models the detector non-uniformities which result in unknown offset error in the projections. '0' means no ring correction. '1' enables ring correction. '2' uses the improved ring correction by enforcing a zero mean constraint on the offset errors. '3' uses a more improved ring correction by enforcing a zero constraint on the weighted mean of offset errors over overlapping rectangular patches. If used, the ring artifacts in the reconstruction should reduce.
	quad_convex - Legal values are '0' and '1'. If '1', then the algorithm uses a convex quadratic forward model. This model does not account for the zinger measurements which causes streak artifacts in the reconstruction. If '0', then the algorithm uses a generalized Huber function which models the effect of zingers. This reduces streak artifacts in the reconstruction. Also, using '1' disables estimation of variance parameter 'sigma' and '0' enables it.
	huber_delta - The parameter \delta of the generalized Huber function which models the effect of zingers. Legal values are in the range 0 to 1.
	huber_T - The threshold parameter T of the generalized Huber function. All positive values are legal values.
	restart - If the reconstruction gets killed due to any unfortunate reason (like exceeding walltime in a super-computing cluster), use this flag to restart the reconstruction from the beginning of the last run multi-resolution stage. Don't use restart if WRITE_EVERY_ITER is 1.


