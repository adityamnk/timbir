==================
Install directions
==================

This section covers the basics of how to download and install TIMBIR.

.. contents:: Contents:
   :local:

Software Dependencies
=====================

    MPI compiler (Intel or Open MPI) 
    OpenMP
    make utility

Compiling the code
==================

1. Compile the MBIR algorithm code in ./Source_Code/MBIR_4D.
   For more info, read the README in ./Source_Code/MBIR_4D.
	> cd ./Source_Code/MBIR_4D
	> make	#Generates library files in ./Source_Code/lib

2. Run the reconstruction algorithm 
	a. If the input data format is standard binary, run the code in ./Source_Code/reconstruct/basic.
	   For more info on data format and running the code, read the README in ./Source_Code/reconstruct/basic.
		> cd ./Source_Code/reconstruct/basic
		> make #Generates executables in the same folder
	b. If the input data is in HDF format used at APS, run the code in ./Source_Code/reconstruct/read_data/APS_Data_Format
	   For more info on data format and running the code, read the README in ./Source_Code/reconstruct/read_data/APS_Data_Format
		> cd ./Source_Code/reconstruct/real_data/APS_Data_Format
		> make #Generates executables in the same folder
	c. If the input data is in standard HDF format, run the code in ./Source_Code/reconstruct/read_data/Standard_Data_Format
	   For more info on data format and running the code, read the README in ./Source_Code/reconstruct/read_data/Standard_Data_Format
		> cd ./Source_Code/reconstruct/real_data/Standard_Data_Format
		> make #Generates executables in the same folder
