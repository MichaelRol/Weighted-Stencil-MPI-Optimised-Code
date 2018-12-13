#
# Makefile to build example MPI programs 
#

EXE=stencil.exe

EXES=$(EXE)

CFLAGS=-Wall -g -DDEBUG -std=c99 -Ofast -march=native

all: $(EXES)

$(EXES): %.exe : %.c
	mpicc $(CFLAGS) -o $@ $^

.PHONY: clean all

clean:
	\rm -f $(EXES) 
	\rm -f *.o
