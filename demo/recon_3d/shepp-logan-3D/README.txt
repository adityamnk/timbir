
Reconstruction of Shepp-Logan phantom

** SUMMARY **
This is a simple example to demonstrate the reconstruction of a Shepp-Logan phantom from simulated data.

** HOW? **
- Compile code by running 'make' on the command line in /src/MBIR_4D/ and /src/reconstruct/bin_data/ 
- Read /src/MBIR_4D/README and /src/reconstruct/bin_data/README for instructions to run the reconstruction code and to understand the command line options.
- Run the script run_unix.sh if running on a unix/linux/Mac computer. Use run_cluster.sh if running on the Purdue rice super-computing cluster.
	It uses 'mpiexec' as the MPI launcher. Otherwise, use the appropriate MPI launcher in your environment.
	NOTE: THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
- The matlab code which generates the binary data files is gen_data_shepplogan.m	
