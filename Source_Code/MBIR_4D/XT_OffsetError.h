#ifndef XT_OFFSETERROR_H
#define XT_OFFSETERROR_H
void gen_offset_constraint_windows (Sinogram* SinogramPtr, TomoInputs* TomoInputsPtr);
void constrained_quad_opt (Real_t** Lambda, Real_t** b, Real_arr_t*** A, Real_arr_t** x, int32_t Nr, int32_t Nt, int32_t M, TomoInputs* TomoInputsPtr);
void compute_d_constraint (Real_arr_t*** A, Real_arr_t **d, int32_t Nr, int32_t Nt, int32_t M, FILE* debug_file_ptr);
#endif
