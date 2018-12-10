#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

#define NROWS 1024
#define NCOLS 1024
#define MASTER 0

int calc_nrows_from_rank(int rank, int size);

int main(int argc, char* argv[]) {
     
    int x, y; //rows and columns 
    int right;
    int left;
    int size;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_nrows = calc_nrows_from_rank(rank, size);
    int local_ncols = NCOLS;
    float *image;
    float *tmp_image;
    float *sendbuf;
    float *recvbuf;
    //int remote_nrows = calc_nrows_from_rank(size - 1, size);
    //float *printbuf;

    left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
    right = (rank + 1) % size;
    
    image = malloc(sizeof(float) * NCOLS * NROWS);
    tmp_image = malloc(sizeof(float) * NCOLS * NROWS);
    
    sendbuf = malloc(sizeof(float) * local_nrows);
    recvbuf = malloc(sizeof(float) * local_nrows);
    
    //printbuf = (double*)malloc(sizeof(double) * (remote_nrows + 2));

    //Init to 0
    for (y = 0; y < NROWS; y++) {
        for (x = 0; x < NCOLS; y++) {
            image[y*NROWS+x] = 0.0f;
            tmp_image[y*NROWS+x] = 0.0f;
        }
    }

    //Init to checkboard
    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
           for (int jj = j*NROWS/8; jj < (j+1)*NROWS/8; jj++) {
               for (int ii = i*NCOLS/8; ii < (i+1)*NCOLS/8; ii++) {
                   if((i+j)%2) image[jj+ii*NROWS] = 100.0f;
               }
           }
       }
    }



    printf("RANK of Node: %d\nNROWS: %d\nNCOLS: %d\nLROWS: %d\nLCOLS: %d\n", rank, NROWS, NCOLS, local_nrows, local_ncols);

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
