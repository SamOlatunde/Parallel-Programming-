//**************************************************************
// Assignment #2
// Name: Samuel Olatunde, and Michelle Orru
// Parallel Programming Date: 03/29/2024
//***************************************************************
// This is an mpi program that simulates the fft algorithm to 
// compute the real and imaginary values. The routine is run and
// timed 3 times an the average rn time is recorded 
//
// To run: sumbit the script for file using sbatch. The script is 
// in the same folder as the file 
//*****************************************************************
#include<stdio.h>
#define _USE_MATH_DEFINES 
#include<math.h>
#include<mpi.h>

#define N 16384

struct complexNum
{
  double real;
  double imag;
};


void calcCooleyTukey(struct complexNum *FFT, struct complexNum *timeDomain, int Size, int offset);

int main(void) {
    int comm_sz; // number of processes
    int my_rank; // a process's rank 
    double start = 0.0;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Datatype MPI_COMPLEXNUM;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_COMPLEXNUM);
    MPI_Type_commit(&MPI_COMPLEXNUM);

    int local_N = N / comm_sz;
    struct complexNum local_FFT[local_N];
    struct complexNum FFT[N];
    struct complexNum timeDomain[N];

   //Hard Coding first 8 entries of time domain 
    timeDomain[0].real = 3.6;
    timeDomain[0].imag = 2.6;
    timeDomain[1].real = 2.9;
    timeDomain[1].imag = 6.3;
    timeDomain[2].real = 5.6;
    timeDomain[2].imag = 4;
    timeDomain[3].real = 4.8;
    timeDomain[3].imag = 9.1;
    timeDomain[4].real = 3.3;
    timeDomain[4].imag = 0.4;
    timeDomain[5].real = 5.9;
    timeDomain[5].imag = 4.8;
    timeDomain[6].real = 5;
    timeDomain[6].imag = 2.6;
    timeDomain[7].real = 4.3;
    timeDomain[7].imag = 4.1;

    // Initializing the other entries of time domain to 0    
    for (int i = 8; i < N; i++) {
        timeDomain[i].real = 0.0;
        timeDomain[i].imag = 0.0;
    }

    if (my_rank == 0) {
        // Start timer before computation
        start = MPI_Wtime();
    }

    calcCooleyTukey(local_FFT, timeDomain, local_N, my_rank * local_N);

    MPI_Gather(local_FFT, local_N, MPI_COMPLEXNUM, FFT, local_N, MPI_COMPLEXNUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        // Timing variables 
        double finish = 0.0; 

        // Take time after we perform last computation 
        finish = MPI_Wtime();

        // Print time elapsed 
        printf("Time Elapsed = %f microseconds\n\n\n", (finish-start)*1000000);

        printf("TOTAL PROCESSED SAMPLES: %d\n", N);
        printf("===========================================\n");
        printf("XR[0]: %f\t", FFT[0].real);
        printf("XI[0]: %f\n", FFT[0].imag);
        printf("===========================================\n");
        for(int i = 1; i <= 10; i++) {
            printf("XR[%d]: %f\tXI[%d]: %f\n", i, FFT[i].real, i, FFT[i].imag);
        }
        printf("===========================================\n\n");
    }

    MPI_Type_free(&MPI_COMPLEXNUM);
    MPI_Finalize();
    
    return 0; 
}




//*******************************************************************
// Name::calcCooleyTukey()
// Parameters: 2 complexNum Pointers
// This function computes the real and imaginary components of the 
// of the fft cooeficients 
//********************************************************************
void calcCooleyTukey(struct complexNum *FFT, struct complexNum *timeDomain, int Size, int offset) {
    double real, imag, theta;

    for (int k = 0; k < Size; k++) {
        real = imag = 0.0;
       

        for(int n = 0; n < N; n ++) {
            theta = ((2 * M_PI) * n * (k + offset)) / N;
            real += (timeDomain[n].real * cos(theta)) + (timeDomain[n].imag * sin(theta));
            imag += (-(timeDomain[n].real) * sin(theta)) + (timeDomain[n].imag * cos(theta));


        }

        FFT[k].real = real;
        FFT[k].imag = imag;

    }
}