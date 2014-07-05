CC=gcc
#Leave out -fopenmp for now until tested
CFLAGS=-Wall -O3
CFLAGSDBG=-Wall -O0 -g
LDFLAGS=-lm

all: trapFind

trapFind: trapFind.o checkPaths.o datastruct.o
	$(CC) $(CFLAGS) trapFind.o checkPaths.o datastruct.o -o trapFind $(LDFLAGS) 

#Why do we have problems in the makefile but not the command line?
#Maybe try a different make system
trapFind.o: trapFind.c trapFind.h checkPaths.h checkPaths.c
	$(CC) $(CFLAGS) trapFind.c trapFind.h checkPaths.h checkPaths.c $(LDFLAGS)

checkPaths.o: checkPaths.c checkPaths.h
	$(CC) $(CFLAGS) checkPaths.c checkPaths.h $(LDFLAGS)

datastruct.o: datastruct.c datastruct.h
	$(CC) $(CFLAGS) datastruct.c datastruct.h $(LDFLAGS)
