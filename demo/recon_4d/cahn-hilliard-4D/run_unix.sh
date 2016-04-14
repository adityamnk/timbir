
cd ../../../src/MBIR_4D/
make clean
make

cd ../reconstruct/bin_data/
make clean
make

cd ../../../demo/recon_4d/cahn-hilliard-4D/

# THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
mpiexec -n 1 ../../../src/reconstruct/bin_data/XT_Main --proj_rows 4 --proj_cols 256 --proj_num 256 --recon_num 8 --vox_wid 1 --rot_center 132 --sig_s 0.02 --sig_t 0.0005 --c_s 0.00001 --c_t 0.00001 --convg_thresh 1 --remove_rings 3 
