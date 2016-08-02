==================
Install directions
==================

This section covers the basics of how to download and install TIMBIR.

.. contents:: Contents:
   :local:

Dependencies
============

- MPI compiler (Intel or Open MPI) 
- OpenMP
- make utility
- HDF library (optional)

Dependencies Install
=====================
- Install MPI ( Here I use openMPI as a example)::

   $ wget https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.2.tar.gz
   $ tar xvzf openmpi-1.10.2.tar.gz
   $ cd <openmpi path>
   $ ./configure --prefix=<your mpi install path>
   $ make all install
   
- Install HDF5 library. First go to https://www.hdfgroup.org/HDF5/release/obtainsrc.html#conf to find the appropriate hdf5 library for your platform (here we use hdf5-1.8.16.tar as example::

   $ tar xvf hdf5-1.8.16.tar
   $ cd <hdf5 path>
   $ ./configure --prefix=/clhome/KYUE/lib/hdf5 --enable-fortran --enable-cxx
   $ make
   $ make install

- Set your library path with HDF5 library and MPI library (here we use bash as example)::

   $ vi env.sh (create a bash script)
   $ export HDF5_BASE=<hdf5 full path>
   $ export MPI_BASE=<MPI full path>
   $ export PATH = ${MPI_BASE}/bin:${HDF5_BASE}/bin:$PATH
   $ export LD_LIBRARY_PATH= = ${MPI_BASE}/lib:${HDF5_BASE}/lib64:$PATH
   $ source env.sh

Compiling TIMBIR
================

To compile the MBIR algorithm code::

   $ git clone https://github.com/adityamnk/timbir.git timbir
   $ cd timbir/src/MBIR_4D
   $ make clean
   $ make

This generates library files in timbir/src/lib. For more information, read the README in timbir/src/MBIR_4D.

Running the reconstruction algorithm
====================================
 
If the input data format is a standard binary, compile and run the code in timbir/src/reconstruct/bin_data.
For more information on data format and running the code, read the README in timbir/src/reconstruct/bin_data::

   $ cd timbir/src/reconstruct/bin_data
   $ make clean
   $ make

This generates executables in the same folder.

If the input data is in HDF format used at APS, compile and run the code in timbir/src/reconstruct/aps_data.
For more information on data format and running the code, read the README in timbir/src/reconstruct/aps_data::

   $ cd timbir/src/reconstruct/aps_data
   $ make clean
   $ make

This generates executables in the same folder.

If the input data is in standard HDF format, compile and run the code in timbir/src/reconstruct/std_data.
For more informantion on data format and running the code, read the README in timbir/src/reconstruct/std_data::

   $ cd timbir/src/reconstruct/std_data
   $ make clean
   $ make #Generates executables in the same folder
