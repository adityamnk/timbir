
cd ../../../src/MBIR_4D/
make

cd ../reconstruct/bin_data/
make

cd ../../../demo/recon_3d/shepp-logan-3D/

# THE REGULARIZATION PARAMETERS USED IN THIS EXAMPLE ARE NOT OPTIMAL. IT IS JUST A WORKING CASE.
mpiexec -n 1 ../../../src/reconstruct/basic/XT_Main --proj_rows 16 --proj_cols 256 --proj_num 256 --recon_num 1 --vox_wid 1 --rot_center 128 --sig_s 0.01 --sig_t 1 --c_s 0.00001 --c_t 0.00001 --convg_thresh 1 --remove_rings 0 --quad_convex 
