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


CC = mpicc
CFLAGS = -Wall -fopenmp -O3 -ansi
EXEC = XT_Main
OBJS = XT_Main.o XT_Views.o XT_HDFIO.o allocate.o
SRCS = XT_Main.c XT_Views.c XT_HDFIO.c allocate.c
LIBFLAGS = -lmbir4d -lm -lhdf5
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJS) $(LIBFLAGS) -L../../../lib/
$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -c $(SRCS) -I../../../include/
clean:
	rm -f $(EXEC)
	rm -f *.o
	rm -f *.gch
