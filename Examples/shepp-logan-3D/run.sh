#!/bin/tcsh 
# Example file for 3DMBIR code. Reconstruction a few slices from a data set from BL 8.3.2. by sub-setting the views. 

#PBS -q regular
#PBS -l mppwidth=96
#PBS -l walltime=4:00:00
#PBS -N PCT1
#PBS -e PCT1.$PBS_JOBID.err
#PBS -o PCT1.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
setenv OMP_NUM_THREADS 24
setenv CRAY_ROOTFS DSL
module load PrgEnv-intel
module load python/2.7.3
module load h5py
module load pil
module load mpi4py

# THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
aprun -j 2 -n 1 -N 1 -d 24 -S 1 -ss -cc numa_node ../../Source_Code/reconstruct/basic/XT_Main --proj_rows 16 --proj_cols 256 --proj_num 256 --recon_num 1 --vox_wid 1 --rot_center 128 --sig_s 0.01 --sig_t 1 --c_s 0.00001 --c_t 0.00001 --convg_thresh 1 --remove_rings 0 
