#
# Makefile to build example MPI programs 
#

EXE=skeleton2-heated-plate.exe

EXES=$(EXE1)

CFLAGS=-Wall -g -DDEBUG

all: $(EXES)

$(EXES): %.exe : %.c
	mpicc $(CFLAGS) -o $@ $^

.PHONY: clean all

clean:
	\rm -f $(EXES) 
	\rm -f *.o
