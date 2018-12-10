#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

#define NROWS 1024
#define NCOLS 1024
#define MASTER 0

int calc_nrows_from_rank(int rank, int size);

int main(int argc, char* argv[]) {
     
    int right;
    int left;
    int size;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_nrows = calc_nrows_from_rank(rank, size);
    int local_ncols = NCOLS;
    double **grid;
    double **grid_next;
    double *sendbuf;
    double *recvbuf;
        

    left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
    right = (rank + 1) % size;
    
     
    printf("RANK: %d\nNROWS: %d\nNCOLS: %d\nLROWS: %d\nLCOLS: %d\n", rank, NROWS, NCOLS, local_nrows, local_ncols);

    MPI_Finalize();

    return EXIT_SUCCESS;
    

}

int calc_nrows_from_rank(int rank, int size) {
    int nrows;
    
    nrows = NROWS/size;
    if ((NROWS % size) != 0) {
        nrows += NROWS % size;
    }

  return nrows;

}
