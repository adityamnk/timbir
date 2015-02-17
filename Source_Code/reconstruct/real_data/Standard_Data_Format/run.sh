#mpiexec -n 1 ./XT_Main --proj_rows 16 --proj_cols 768 --proj_num 6144 --recon_num 128 --vox_wid 1.3 --rot_center 385 --sig_s 0.010413 --sig_t 0.000104 --c_s 0.000001 --c_t 0.000001 --convg_thresh 1 --remove_rings 1 --remove_streaks 1   
cd $PBS_O_WORKDIR
module load devel
export PARALLEL=1
export OMP_NUM_THREADS=32
uniq < $PBS_NODEFILE > nodefile

mpiexec -n 2 -machinefile nodefile ./XT_Main --path2data /scratch/conte/m/mohank/Argonne_Datasets/APS_Beamtime_042014/AlCu_14_slice_0_31_Std.hdf --proj_rows 8 --datafile_row0 4 --proj_cols 768 --proj_start 0 --proj_num 1536 --recon_num 32 --vox_wid 1.3 --rot_center 385 --sig_s 0.010413 --sig_t 0.000104 --c_s 0.000001 --c_t 0.000001 --convg_thresh 1 --remove_rings 2 --remove_streaks 1
