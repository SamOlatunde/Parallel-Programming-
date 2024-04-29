//*****************************************************************
// End of Semester Project
// Name: Samuel Olatunde , and Sunil Rasaily
// GPU Programming Date: Date of Submission (11/28/2022)
//******************************************************************
// This solution uses global synchronization. The down side is that
// we are limited in the number execution resources available
// because we need all blocks/threads to be in the gpu at the time 
// of sync so as to avoid deadlock 
//******************************************************************
#include<iostream>
#include<cuda.h>
#include<cstdlib>
#include<cuda_runtime.h>
#include"timer.h"
#include <cooperative_groups.h>

namespace cg = cooperative_groups;
//error tolerance
const float eT  = 0.00001;

// data size
#define N  32

// limit for the max number of iterations
#define limit 100


#define blocks 1
#define threads 32

dim3 gridSize(blocks, blocks, 1);
dim3 blockSize(threads, threads, 1);


void print(float * a)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f\t", a[i* N + j]);
        }

        printf("\n\n");
    }
}

double checkSum(float * h)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum += h[i*N+j];
        }
    }

    return sum;
}

__device__ int Continue = 0;

// Function Prototypes
__global__ void calcIntTempDistribution(float * h,float *g);
__device__ int converged (float newValue, float oldValue );
void initMetalPlate(float *h, float * g,float edgeTemp);

int main()
{
   // variable declarations
   float * h, *g ;

   //Allocate dynamic memory in host
   h = (float *) malloc ((N*N) * sizeof(float));
   g = (float *) malloc((N*N) * sizeof(float));

   float edgeTemp =70.5;//300;
   double tStart = 0.0, tStop = 0.0, tElapsed = 0.0;

   //initialize matrix
   initMetalPlate(h,g, edgeTemp);

   //device variables
   float *  hd,  * gd;
   long long int size = (N*N) * sizeof(float);

   // allocate space, and copy data
   cudaMalloc((void**) & hd, size);
   cudaMemcpy(hd, h, size, cudaMemcpyHostToDevice);
   cudaMalloc((void**) & gd, size);
   cudaMemcpy(gd, g, size, cudaMemcpyHostToDevice);

   //set grid dimensions
   /*dim3 dimGrid(8,8,1);
   dim3 dimBlock(32,32,1);*/

   //int THREADS;
   //int BLOCKS; 
   
   //cudaOccupancyMaxPotentialBlockSize(&BLOCKS,&THREADS, calcIntTempDistribution, 0,0);
   //printf("Number of Blocks %d\n", BLOCKS);
   //printf("Number of Threads %d\n", THREADS);
  
   // kernel launch and timing
   GET_TIME(tStart);
   calcIntTempDistribution<<<gridSize, blockSize >>>(hd,gd);
   cudaDeviceSynchronize();
   GET_TIME(tStop);

   cudaMemcpy(h, hd, size, cudaMemcpyDeviceToHost);

   // Compute how long it took
   tElapsed = tStop - tStart;

   printf("The code to be timed took %e seconds\n", tElapsed);
   printf("checkSum: %f\n", checkSum(h));
   //print(h);

   //free global memory
   cudaFree(hd);
   cudaFree(gd);

   // Dellocate dynamic memory
   free(h);
   free(g);

    return 0;
}



//*******************************************************************
// Name::calcIntTempDistribution()
// Parameters: 2 float pointers
//
//********************************************************************
__global__ void calcIntTempDistribution(float * h,float *g)
{
   int iteration = 0;
   int row = blockDim.y * blockIdx.y + threadIdx.y;
   int col = blockDim.x * blockIdx.x + threadIdx.x;

   //number of threads in grid
   int totalThreadsX = blockDim.x * gridDim.x;
   int totalThreadsY = blockDim.y * gridDim.y;
   
   //auto grid = cg::this_grid();

   do
   {
        for (int i = row; i < N; i += totalThreadsY)
        {
            for (int j = col; j < N; j += totalThreadsX)
            {
                //Takes care of boundary points
                if (i != 0 && i != (N - 1) && j != 0 && j != (N - 1))
                {
                    g[i * N + j] = 0.25 * (h[(i - 1) * N + j] + h[(i + 1) * N + j] +
                                            h[i * N + j - 1] + h[i * N + j + 1]);
                }

            }
        }
        
        //synchronize the whole grid 
       // grid.sync();

        Continue = 0;

        for (int i = row; i < N; i += totalThreadsY)
        {
            for (int j = col; j < N; j += totalThreadsX)
            {
                if (converged(g[i * N + j], h[i * N + j]) == 0)
                {
                    Continue = 1;
                }

                h[i * N + j] = g[i * N + j];
            }
        }
        
        //sync grid again to avoid race condition
        // grid.sync();

     iteration++;
   }while(Continue == 1 && iteration < limit);
   //printf("Blah");
   //printf("%d\n", iteration);
}

//*******************************************************************
// Name::converged()
// Parameters: 2 floats
// Tests for convergence of two points. Returns true if the error is
// within error tolerance; false otherwise
//********************************************************************
// bool converged (float newValue, float oldValue )
 __device__ int converged (float newValue, float oldValue )
{
    float er = (newValue-oldValue)/newValue;
    //printf("er %f\n", er);
    if (er < 0) er = -er;

    if (er <= eT)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


//*******************************************************************
// Name::initMetalPlate()
// Parameters: 1 2d float array, 1 float
// Initializes the metal sheet with the intial values of the edges
// and guess values for interior points
//********************************************************************
void initMetalPlate(float *h, float * g,float edgeTemp)
{
   //we reduce the temparture by this value with every
   // outer loop iteration
   float reduceFactor = edgeTemp/N;

   int row = 0;
   int col;

   for( int i = 0; i < (N/2); i++)
   {

        row = i;
        for (col = i; col < N-i; col++ )
        {
            h[row * N + col] = edgeTemp;
            g[row * N + col] = edgeTemp;
        }

        col--;

        for (row = row+1; row < N-i; row++)
        {
            h[row * N + col] = edgeTemp;
            g[row * N + col] = edgeTemp;
        }

        row--;
        col--;

        for(col = col; col >=i; col--)
        {
           h[row * N + col] = edgeTemp;
           g[row * N + col] = edgeTemp;
        }

        row--;
        col++;

        for(row = row; row >i; row--)
        {
            h[row * N + col] = edgeTemp;
            g[row * N + col] = edgeTemp;
        }

      edgeTemp = edgeTemp - reduceFactor;
    }

    // print(g);
    // printf("\n\n");

}