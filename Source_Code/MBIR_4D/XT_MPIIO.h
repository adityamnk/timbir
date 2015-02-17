#ifndef XT_MPIIO_H

#define XT_MPIIO
int32_t read_SharedBinFile_At (char filename[], Real_arr_t* data, int32_t offset, int32_t size, FILE* debug_file_ptr);
/*void read_SharedBinFile (char filename[], Real_t* data, int32_t size);*/
int32_t write_SharedBinFile_At (char filename[], Real_arr_t* data, int32_t offset, int32_t size, FILE* debug_file_ptr);
/*void write_SharedBinFile (char filename[], Real_t* data, int32_t size);*/
#endif
