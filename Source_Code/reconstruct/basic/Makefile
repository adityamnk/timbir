CC = mpicc
CFLAGS = -Wall -fopenmp -O3 -ansi
EXEC = XT_Main
OBJS = XT_Main.o
SRCS = XT_Main.c
HEADERS = XT_Main.h
LIBFLAGS = -lmbir4d -lm
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJS) $(LIBFLAGS) -L../../lib/
$(OBJS): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRCS) -I../../include/
clean:
	rm -f $(EXEC)
	rm -f *.o
	rm -f *.gch
