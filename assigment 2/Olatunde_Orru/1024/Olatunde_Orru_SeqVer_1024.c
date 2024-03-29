//**************************************************************
// Assignment #2
// Name: Samuel Olatunde, and Michelle Orru
// Parallel Programming Date: 03/29/2024)
//***************************************************************
// This is a serial program that simulates the fft algorithm to 
// compute the real and imagin
//
// To run: sumbit the script for file using sbatch. The script is 
// in the same folder as the file 
//*****************************************************************
#include<stdio.h>
#define _USE_MATH_DEFINES 
#include<math.h>
#include<mpi.h>

#define N 1024

struct complexNum
{
  double real;
  double imag;
};


void calcCooleyTukey(struct complexNum * FFT, struct complexNum * timeDomain);


int main()
{
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


    // initializing the other entries of time domain to 0    
    for (int i = 8; i < N; i++)
    {
        timeDomain[i].real = 0;
        timeDomain[i].imag = 0;
    }
 
    
    return 0; 
}


void calcCooleyTukey(struct complexNum * FFT, struct complexNum * timeDomain)
{

 

}