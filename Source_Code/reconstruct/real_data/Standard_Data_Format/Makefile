CC = mpicc
CFLAGS = -Wall -fopenmp -O3 -ansi
EXEC = XT_Main
OBJS = XT_Main.o XT_HDFIO.o allocate.o
SRCS = XT_Main.c XT_HDFIO.c allocate.c
LIBFLAGS = -lmbir4d -lm -lhdf5
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJS) $(LIBFLAGS) -L../../../lib/
$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -c $(SRCS) -I../../../include/
clean:
	rm -f $(EXEC)
	rm -f *.o
	rm -f *.gch
