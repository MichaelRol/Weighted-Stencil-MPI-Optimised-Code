#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h> 
#include "mpi.h"

#define MASTER 0
#define OUTPUT_FILE "stencil.pgm"

int calc_nrows_from_rank(int rank, int size, int ny);
void output_image(const char * file_name, const int nx, const int ny, float * restrict image);
void init_image(const int nx, const int ny, float * restrict image, float * restrict tmp_image);
void stencil(const int nx, const int ny, float *  image, float * tmp_image, int firstcell, int lastcell, float * sendbuf, float * recvbuf, int above, int below, MPI_Status status);
double wtime(void);

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
    int above;
    int below;
    int size;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    int local_nrows = calc_nrows_from_rank(rank, size, ny);
    int local_ncols = nx;
    float *image;
    float *tmp_image;
    float *sendbuf;
    float *recvbuf;
    int firstcell = rank * local_nrows;
    int lastcell;

    //calculate the starting and ending cells 
    if (rank == size - 1) {
        lastcell = ny - 1;
    } else {
        lastcell = firstcell + local_nrows - 1;
    }
    
    //calculate the rank above and below
    above = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
    below = (rank + 1) % size;
    
    //printf("My Rank is: %d. Rows per rank: %d, Start: %d, End: %d\n\n", rank, local_nrows, firstcell, lastcell);
    
    //allocate memory for images and bufferers
    image = malloc(sizeof(float) * nx * ny);
    tmp_image = malloc(sizeof(float) * nx * ny);
    
    sendbuf = malloc(sizeof(float) * local_ncols);
    recvbuf = malloc(sizeof(float) * local_ncols);

    init_image(nx, ny, image, tmp_image);
    double tic = wtime();
    //for (int t = 0; t < niters; t++) {
        stencil(nx, ny, image, tmp_image, firstcell, lastcell, sendbuf, recvbuf, above, below, status);
      //  stencil(nx, ny, tmp_image, image, firstcell, lastcell, sendbuf, recvbuf, above, below, status);
   // }
    double toc = wtime();

    if (rank == MASTER) output_image(OUTPUT_FILE, nx, ny, image);
    MPI_Finalize();
    return EXIT_SUCCESS;
     

}

void stencil(const int nx, const int ny, float *  image, float * tmp_image, int firstcell, int lastcell, float * sendbuf, float * recvbuf, int above, int below, MPI_Status status) {
    //send top row above
    for (int i = 0; i < nx; i++) {
        sendbuf[i] = image[firstcell * nx + i];
    }
    MPI_Send(sendbuf, nx, MPI_FLOAT, above, 123, MPI_COMM_WORLD);
    MPI_Recv(recvbuf, nx, MPI_FLOAT, below, 123, MPI_COMM_WORLD, &status);

    for (int x = 0; x < nx; x++) {
        printf("Row: %d   Value: %f\n", x, recvbuf[x]);
    }

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

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
