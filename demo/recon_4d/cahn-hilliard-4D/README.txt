
Reconstruction of a phantom generated using the Cahn-Hilliard equation

** SUMMARY **
This is a simple example to demonstrate 4D reconstruction of a time varying phantom using simulated data.

** HOW? **
- Compile code by running 'make' on the command line in /src/MBIR_4D/ and /src/reconstruct/bin_data/ 
- Read /src/MBIR_4D/README and /src/reconstruct/bin_data/README for instructions to run the reconstruction code and to understand the command line options.
- Run the script run_unix.sh if running on a unix/linux/Mac computer. Use run_cluster.sh if running on the Purdue rice super-computing cluster.
	It uses 'mpiexec' as the MPI launcher. Otherwise, use the appropriate MPI launcher in your environment.
	NOTE: THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
