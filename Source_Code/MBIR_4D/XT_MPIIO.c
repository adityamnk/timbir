#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "XT_Structures.h"
/*Reads 'size' number of elements from the binary data file with name 'filename' 
 starting at 'offset'. The data read is stored in 'data'. */
int32_t read_SharedBinFile_At (char filename[100], Real_arr_t* data, int32_t offset, int32_t size, FILE* debug_file_ptr)
{
	MPI_File fh;
	MPI_Status status;
	char BinFilename[100];
	int32_t len;

    	sprintf(BinFilename, "%s.bin", filename);
	MPI_File_open(MPI_COMM_WORLD, BinFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at(fh, offset*sizeof(Real_arr_t), data, size, MPI_REAL_ARR_DATATYPE, &status);
	MPI_Get_count(&status, MPI_REAL_ARR_DATATYPE, &len);
	MPI_File_close(&fh);
    	if(len == MPI_UNDEFINED || len != size)
	{
		fprintf (debug_file_ptr, "ERROR: read_SharedBinFile_At: Read %d number of elements from the file %s at an offset of %d bytes.\n. However, required number of elements is %d.", len, filename, offset, size);
		return(-1);
	}
	return(0);
/*	else
		fprintf (debug_file_ptr, "read_SharedBinFile_At: Read %d number of elements from the file %s at an offset of %d bytes.\n", size, filename, offset);*/
}

/*Reads 'size' number of elements from the binary data file with name 'filename' 
 starting at the offset rank*size. The data read is stored in 'data'. */
/*void read_SharedBinFile (char filename[], Real_t* data, int32_t size)
{
	int32_t rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	read_SharedBinFile_At (filename, data, rank*size, size);
}*/

/*Writes 'size' number of elements to the binary data file with name 'filename' 
 starting at 'offset'. The data is written from the array 'data'. */
int32_t write_SharedBinFile_At (char filename[], Real_arr_t* data, int32_t offset, int32_t size, FILE* debug_file_ptr)
{
	MPI_File fhw;
	MPI_Status status;
	char BinFilename[100];
	int32_t len;
    	sprintf(BinFilename, "%s.bin", filename);
	MPI_File_open(MPI_COMM_WORLD, BinFilename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
	MPI_File_write_at(fhw, offset*sizeof(Real_arr_t), data, size, MPI_REAL_ARR_DATATYPE, &status);
	MPI_Get_count(&status, MPI_REAL_ARR_DATATYPE, &len);
	MPI_File_close(&fhw);
    	if(len == MPI_UNDEFINED || len != size)
	{
		fprintf (debug_file_ptr, "ERROR: write_SharedBinFile_At: Wrote %d number of elements to the file %s at an offset of %d bytes.\n. However, actual number of elements to be written is %d.", len, filename, offset, size);
		return(-1);
	}
	return(0);
/*	else
		fprintf (debug_file_ptr, "read_SharedBinFile_At: Read %d number of elements from the file %s at an offset of %d bytes.\n", size, filename, offset);*/
}

