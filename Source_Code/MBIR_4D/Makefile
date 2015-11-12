 # ============================================================================
 # Copyright (c) 2015 K. Aditya Mohan (Purdue University)
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without modification,
 # are permitted provided that the following conditions are met:
 #
 # Redistributions of source code must retain the above copyright notice, this
 # list of conditions and the following disclaimer.
 #
 # Redistributions in binary form must reproduce the above copyright notice, this
 # list of conditions and the following disclaimer in the documentation and/or
 # other materials provided with the distribution.
 #
 # Neither the name of K. Aditya Mohan, Purdue
 # University, nor the names of its contributors may be used
 # to endorse or promote products derived from this software without specific
 # prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 # FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 # DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 # SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 # CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 # USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

CC = mpicc #MPI compiler
CFLAGS = -Wall -O3 -fopenmp -ansi -lm #Compiler flags. Can change the code optimization levels.
#CFLAGS = -g -Wall -fopenmp -ansi -lm #Compiler flags. Can change the code optimization levels.
#MACROS = -DNO_COST_CALCULATE -DNO_DEBUG #Macros
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
