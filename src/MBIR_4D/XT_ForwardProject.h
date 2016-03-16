
void forward_project_voxel (Sinogram* SinogramPtr, Real_t voxel_val, Real_arr_t*** ErrorSino, AMatrixCol* AMatrixPtr, /*AMatrixCol* VoxelLineResponse,*/ int32_t sino_idx, int32_t slice);
void forward_project_voxel_AMat1D (Sinogram* SinogramPtr, Real_t voxel_val, Real_arr_t*** ErrorSino, AMatrixCol* AMatrixPtr, /*AMatrixCol* VoxelLineResponse,*/ int32_t sino_idx, int32_t slice);
void forward_project_voxel_AMat2D (Sinogram* SinogramPtr, Real_t voxel_val, Real_arr_t*** ErrorSino, AMatrixCol* AMatrixPtr, /*AMatrixCol* VoxelLineResponse,*/ int32_t sino_idx, int32_t slice);
