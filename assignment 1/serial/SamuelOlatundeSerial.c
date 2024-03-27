// Samuel Olatunde 
// Parallel Programming - Spring 24, Dr. Colmenares 
// This is a serial programming that computes the sum of all the
// elements of an array 
#include <stdio.h>
#include"timer.h"
#include<stdlib.h>
#define N  640000
#define noProcesses 1

// function prototypes 
void Sum(int * A, long long int * sum);
void init_routine(int * array);

int main(void) 
{
  // allocate memory
  int * A = (int *)malloc( N * sizeof(int));
   
  // initialize allocated memory
  init_routine(A);

  // timing variables 
  double tstart = 0.0, tend= 0.0, telapsed = 0.0; 
  long long int sum = 0; // accumulator 
 
  GET_TIME(tstart); // get time right before entering routine 
  Sum(A, &sum);
  GET_TIME(tend); // time right after routine 

 // compute time elapsed
 telapsed = tend-tstart; 

 printf("Sum = %lld \n" , sum);
 printf("The routine took %f microseconds", telapsed);

  free(A);

  return 0;
}

// Name: Sum 
// Purpose: Add up the contents of array and store in the
//           memory addrees that "sum" points to
void Sum(int * A, long long int * sum)
{
  for (int i = 0; i< N;  i++)
    {
      *sum += A[i];
    }
}

// Name: init_routine
// Purpose: inititalize array with the first N consecutive numbers 
void init_routine(int * array)
{
 
   for(int i= 0; i < N; i++ )
     {
       array[i] = i+1;
     }
}