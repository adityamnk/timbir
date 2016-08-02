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

-------------
Demo/Examples
-------------
- Unix/Linux/Mac OS
	- 3D Reconstruction: Change directory to demo/recon_3d/shepp-logan-3D/ and run the script run_unix.sh
	- 4D Reconstruction: Change directory to demo/recon_4d/cahn-hilliard-4D/ and run the script run_unix.sh
	
- Supercomputing cluster (Rice cluster at Purdue)
	- 3D Reconstruction: Change directory to demo/recon_3d/shepp-logan-3D/ and run the script run_cluster.sh
	- 4D Reconstruction: Change directory to demo/recon_4d/cahn-hilliard-4D/ and run the script run_cluster.sh

---------
Citation
---------
- \K. A. Mohan, S. V. Venkatakrishnan, J. W. Gibbs, E. B. Gulsoy, X. Xiao, M. De Graef, P. W. Voorhees, and C. A. Bouman, "TIMBIR: A method for time-space reconstruction from interlaced views," IEEE Transactions on Computational Imaging, vol. 1, no. 2, pp. 96.111, June 2015. 
