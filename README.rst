TIMBIR
######

TIMBIR is a method for 4D time-space reconstruction of data acquired using synchrotron X-ray computed tomography.
This software can be used to reconstruct any tomographic data in 3D or 4D.

-----------------
Documentation:
-----------------
`http://timbir.readthedocs.org <http://timbir.readthedocs.org>`_

------------
Dependencies
------------
- MPI compiler
- Open MP
- Make utility
- HDF library (optional)

-------------
Demo/Examples
-------------
- Unix/Linux/Mac OS
	- 3D Reconstruction: Change directory to demo/recon_3d/shepp-logan-3D/ and run the script run_unix.sh
	- 4D Reconstruction: Change directory to demo/recon_4d/cahn-hilliard-4D/ and run the script run_unix.sh
	
- Supercomputing cluster (Rice cluster at Purdue)
	- 3D Reconstruction: Change directory to demo/recon_3d/shepp-logan-3D/ and run the script run_cluster.sh
	- 4D Reconstruction: Change directory to demo/recon_4d/cahn-hilliard-4D/ and run the script run_cluster.sh

----------------
Compiling TIMBIR
----------------
- To compile the MBIR algorithm code, run the following commands in a terminal:
	- git clone https://github.com/adityamnk/timbir.git timbir
	- cd timbir/src/MBIR_4D
	- make clean
	- make
This generates library files in timbir/src/lib. For more information, read the README in timbir/src/MBIR_4D.

------------------------------------
Running the reconstruction algorithm
------------------------------------
- If the input data format is a standard binary, compile and run the code in timbir/src/reconstruct/bin_data. 
	- cd timbir/src/reconstruct/bin_data
	- make clean
	- make
	- # This generates executables in the same folder.
	- # For more information on binary data format and running the code, read the README in timbir/src/reconstruct/bin_data:

- If the input data is in HDF format used at APS, compile and run the code in timbir/src/reconstruct/aps_data. 
	- cd timbir/src/reconstruct/aps_data
	- make clean
	- make
	- # This generates executables in the same folder.
	- # For more information on APS data format and running the code, read the README in timbir/src/reconstruct/aps_data:


If the input data is in standard HDF format, compile and run the code in timbir/src/reconstruct/std_data. 
	- cd timbir/src/reconstruct/std_data
	- make clean
	- make 
	- # Generates executables in the same folder
	- # For more informantion on standard HDF data format and running the code, read the README in timbir/src/reconstruct/std_data:

---------
Citation
---------
- \K. A. Mohan, S. V. Venkatakrishnan, J. W. Gibbs, E. B. Gulsoy, X. Xiao, M. De Graef, P. W. Voorhees, and C. A. Bouman, "TIMBIR: A method for time-space reconstruction from interlaced views," IEEE Transactions on Computational Imaging, vol. 1, no. 2, pp. 96.111, June 2015. 
