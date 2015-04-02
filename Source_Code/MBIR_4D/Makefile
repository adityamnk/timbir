CC = mpicc #MPI compiler
CFLAGS = -Wall -fopenmp -O3 -ansi -lm #Compiler flags. Can change the code optimization levels.
MACROS = -DNO_COST_CALCULATE -DNO_DEBUG #Macros
OBJS = XT_Engine.o allocate.o randlib.o tiff.o XT_AMatrix.o XT_ForwardProject.o XT_ICD_update.o XT_Init.o XT_IOMisc.o XT_MPI.o XT_MPIIO.o XT_NHICD.o XT_Profile.o XT_VoxUpdate.o invert.o nrutil.o XT_OffsetError.o
SRCS = XT_Engine.c allocate.c randlib.c tiff.c XT_AMatrix.c XT_ForwardProject.c XT_ICD_update.c XT_Init.c XT_IOMisc.c XT_MPI.c XT_MPIIO.c XT_NHICD.c XT_Profile.c XT_VoxUpdate.c invert.c nrutil.c XT_OffsetError.c
HEADERS = allocate.h tiff.h TiffWrite.h XT_AMatrix.h XT_ForwardProject.h XT_Init.h XT_MPI.h XT_NHICD.h XT_Structures.h randlib.h TiffUtilities.h typeutil.h XT_Constants.h XT_ICD_update.h XT_IOMisc.h XT_MPIIO.h XT_Profile.h XT_VoxUpdate.h invert.h nrutil.h XT_OffsetError.h
LIB = ../lib/libmbir4d.a #The generated library which contains the MBIR reconstruction routines 
$(LIB): $(OBJS)
	ar rs $(LIB) $(OBJS)
$(OBJS): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) $(MACROS) -c $(SRCS)
clean:
	rm -f *.o
	rm -f *.gch
	rm -f $(LIB)
