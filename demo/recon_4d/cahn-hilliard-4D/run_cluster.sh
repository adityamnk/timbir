module load devel
#module load valgrind
cd $PBS_O_WORKDIR
export PARALLEL=1
export OMP_NUM_THREADS=32
uniq < $PBS_NODEFILE > nodefile

cd ../../../src/MBIR_4D/
make

cd ../reconstruct/basic/
make

cd ../../../demo/recon_3d/shepp-logan-3D/

# THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
mpiexec -n 1 -machinefile nodefile ../../../src/reconstruct/basic/XT_Main --proj_rows 4 --proj_cols 256 --proj_num 256 --recon_num 8 --vox_wid 1 --rot_center 132 --sig_s 0.02 --sig_t 0.0005 --c_s 0.00001 --c_t 0.00001 --convg_thresh 1 --remove_rings 3 
