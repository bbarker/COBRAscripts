CC=gcc
CFLAGS=-Wall -O3
CFLAGSDBG=-Wall -O0 -g
LDFLAGS=-lm

all: trapFind

trapFind: trapFind.o checkPaths.o datastruct.o
	$(CC) $(CFLAGS) $(LDFLAGS) trapFind.o checkPaths.o datastruct.o -o trapFind

#Why do we have problems in the makefile but not the command line?
#Maybe try a different make system
trapFind.o: trapFind.c trapFind.h checkPaths.h
	$(CC) $(CFLAGS) trapFind.c 

checkPaths.o: checkPaths.c checkPaths.h
	$(CC) $(CFLAGS) checkPaths.c

datastruct.o: datastruct.c datastruct.h
	$(CC) $(CFLAGS) datastruct.c
