// Samuel Olatunde 
// how i complied my code: mpicc SamuelOlatundeParallelV1.c  -o myexe 
// Parallel Programming - Spring 24, Dr. Colmenares 
// This is a parallel mpi programming that computes the sume of 
// integers in an array. The main array is hosted in process 0 and 
// and process 0 distributes approriate chunks of data to other 
// processes to compute the partial sum and the processes send 
// the partial sum back to process zero with stores each of these
// partial sums in an array and them sums up the array at the end. 
#include <stdio.h>
#include<mpi.h>
#include<stdlib.h>
#define N  640000


// function prototypes 
void Sum(int * A, int n, long long int * sum);
void SumLLD(long long int * A, int n, long long int * sum);
void init_routine(int * array);

int main(void) 
{
  int comm_sz; // number of processes
  int my_rank; // a processes rank 

  MPI_Init (NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 // calculate number of elements for each process
 int partition_sz = N / comm_sz;

  if(my_rank == 0)
  {
   // timing variables 
   double start = 0.0, finish = 0.0; 

   // hold the sums sent from individual processes
   long long int localSum = 0;

   // allocate memory
   int * A = (int *)malloc( N * sizeof(int));
    
   // initialize allocated memory
   init_routine(A);

   // allocate memory and initialize to zero
   long long int * Summation  = (long long int *)calloc
                                (comm_sz, sizeof(long long int));
   
   // take time once we start computation
   start = MPI_Wtime();
   // compute the some of process 0 parition and save in Summation[0]
   Sum(A, partition_sz, Summation);

   // send paritions of large arrays to respective processes  
   for (int i = 1; i < comm_sz; i++)
   {
      MPI_Send( A + (i * partition_sz), partition_sz, MPI_INT,
               i, 0, MPI_COMM_WORLD);
   }


   // receive summation of corresponing section for each process
   //  and place it in approrpiate position in summation array 
   for (int i= 1; i <comm_sz; i++)
   {
    MPI_Recv(&localSum, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    Summation[i] = localSum;
   }
   
   long long int tot_sum = 0; 
   
   // add up the partial sums 
   SumLLD(Summation, comm_sz, &tot_sum);
   
   // take time after we perform last computation 
   finish = MPI_Wtime();

   // print sum 
   printf("Sum = %lld \n", tot_sum);
   
   // print time elapsed 
   printf("Time Elapsed = %f microseconds", (finish-start)*1000000);

   // dealllocate memory 
   free(A);
   free(Summation);

  }
  else 
  {
    // allocate memory for local array  and initialize to 0
    int * local_A  = (int *)calloc(partition_sz, sizeof(int));

    // receive from process 0 appropriate partition to compute  
    MPI_Recv(local_A, partition_sz, MPI_INT, 0,0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    long long localSum = 0; 
    

    // perform local/partial sum
    Sum(local_A, partition_sz, &localSum);

   // send partial sum back to zero 
    MPI_Send(&localSum, 1, MPI_LONG_LONG, 0, 0,MPI_COMM_WORLD);
    
    // deallocate used memory 
    free(local_A);
  }
  
   
  MPI_Finalize();

  return 0;
}

// Name: Sum 
// Purpose: Add up the contents of  integer array and store in the
//           memory address that "sum" points to
void Sum(int * A, int n, long long int * sum)
{
  for (int i = 0; i< n;  i++)
    {
      *sum += A[i];
    }
}

// Name: SumLLD
// Purpose: Accumulates the contents of long long int array 
//  and store in the memory addrees that "sum" points to
void SumLLD(long long int * A, int n, long long int * sum)
{
  for (int i = 0; i< n;  i++)
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