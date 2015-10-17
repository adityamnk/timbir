cd $PBS_O_WORKDIR
module load devel
export PARALLEL=1
export OMP_NUM_THREADS=32
uniq < $PBS_NODEFILE > nodefile

mpiexec -n 1 -machinefile nodefile ./XT_Main --path2data /scratch/conte/m/mohank/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_14.hdf --path2whites /scratch/conte/m/mohank/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --path2darks /scratch/conte/m/mohank/Argonne_Datasets/APS_Beamtime_042014/APS14_AlCu_19.hdf --proj_rows 8 --datafile_row0 4 --proj_cols 768 --proj_start 33000 --proj_num 1536 --K 32 --N_theta 1536 --r 32 --min_acquire_time 0.0047 --rotation_speed 720 --vox_wid 1.3 --rot_center 385 --sig_s 0.010413 --sig_t 0.000104 --c_s 0.000001 --c_t 0.000001 --convg_thresh 1 --remove_rings 2 --quad_convex
