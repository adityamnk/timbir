module load devel
export PARALLEL=1
export OMP_NUM_THREADS=32
uniq < $PBS_NODEFILE > nodefile

# THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
mpiexec -n 1 -machinefile nodefile ../../Source_Code/reconstruct/basic/XT_Main --proj_rows 4 --proj_cols 256 --proj_num 1024 --recon_num 32 --vox_wid 1 --rot_center 132 --sig_s 0.02 --sig_t 0.001 --c_s 0.00001 --c_t 0.00001 --convg_thresh 1 --remove_rings 3 --remove_streaks 1
