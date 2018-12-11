#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

#define MASTER 0
#define OUTPUT_FILE "stencil.pgm"

int calc_nrows_from_rank(int rank, int size, int ny);
void output_image(const char * file_name, const int nx, const int ny, float * restrict image);
void init_image(const int nx, const int ny, float * restrict image, float * restrict tmp_image);

int main(int argc, char* argv[]) {

    // Check usage
    if (argc != 4) {
        fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int niters = atoi(argv[3]);
    int x, y; //rows and columns 
    int right;
    int left;
    int size;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_nrows = calc_nrows_from_rank(rank, size, ny);
    int local_ncols = nx;
    float *image;
    float *tmp_image;
    float *sendbuf;
    float *recvbuf;
    int firstcell = rank * local_nrows;
    int lastcell = firstcell + local_nrows - 1;
    
    // if (rank == size - 1) {
    //     lastcell = ny - 1;
    // } else {
    //     lastcell = (rank + 1) * local_nrows - 1;
    // }
    //int remote_nrows = calc_nrows_from_rank(size - 1, size, ny);
    //float *printbuf;

    left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
    right = (rank + 1) % size;
    
    printf("My Rank is: %d. Rows per rank: %d, Start: %d, End: %d\n\n", rank, local_nrows, firstcell, lastcell);
    
    image = malloc(sizeof(float) * nx * ny);
    tmp_image = malloc(sizeof(float) * nx * ny);
    
    sendbuf = malloc(sizeof(float) * local_ncols);
    recvbuf = malloc(sizeof(float) * local_ncols);
    
    //printbuf = (double*)malloc(sizeof(double) * (remote_nrows + 2));

    init_image(nx, ny, image, tmp_image);

    if (rank == MASTER) output_image(OUTPUT_FILE, nx, ny, image);
    MPI_Finalize();
    return EXIT_SUCCESS;
     

}

// Create the input image
void init_image(const int nx, const int ny, float * restrict image, float * restrict tmp_image) {

    //Init to 0
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            image[y*ny+x] = 0.0f;
            tmp_image[y*ny+x] = 0.0f;
        }
    }

    //Init to checkboard
    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
           for (int jj = j*ny/8; jj < (j+1)*ny/8; jj++) {
               for (int ii = i*nx/8; ii < (i+1)*nx/8; ii++) {
                   if((i+j)%2) image[jj+ii*ny] = 100.0f;
              }
           }
       }
    }
}

int calc_nrows_from_rank(int rank, int size, int ny) {
    int nrows;
    
    nrows = ny/size;
    if ((ny % size) != 0) {
        if (rank == size - 1) nrows += ny % size;
    }

    return nrows;

}

void output_image(const char * file_name, const int nx, const int ny, float * restrict image) {

  // Open output file
    FILE *fp = fopen(file_name, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
        exit(EXIT_FAILURE);
    }

    // Ouptut image header
    fprintf(fp, "P5 %d %d 255\n", nx, ny);

    // Calculate maximum value of image
    // This is used to rescale the values
    // to a range of 0-255 for output
    float  maximum = 0.0f;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (image[j+i*ny] > maximum)
              maximum = image[j+i*ny];
        }
    }

    // Output image, converting to numbers 0-255
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fputc((char)(255.0*image[j+i*ny]/maximum), fp);
        }
    }

    // Close the file
    fclose(fp);

}
