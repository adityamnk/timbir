#ifndef XT_HDF_IO_H
#define XT_HDF_IO_H

#include "hdf5.h"
int32_t read_data (char data_filename[], char whites_filename[], char darks_filename[], float *projections, float *weights, int32_t datafile_row0, int32_t proj_rows, int32_t proj_cols, int32_t proj_start, int32_t proj_num, FILE* debug_file_ptr);

#endif
