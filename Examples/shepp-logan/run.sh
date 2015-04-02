module load devel
export PARALLEL=1
export OMP_NUM_THREADS=32
uniq < $PBS_NODEFILE > nodefile

mpiexec -n 1 -machinefile nodefile ../../Source_Code/reconstruct/basic/XT_Main --proj_rows 24 --proj_cols 512 --proj_num 256 --recon_num 1 --vox_wid 1 --rot_center 256 --sig_s 0.01 --sig_t 1 --c_s 0.00001 --c_t 0.00001 --convg_thresh 1 --remove_rings 0 --remove_streaks 0
