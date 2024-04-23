// Samuel Olatunde
// Parallel Programming - Spring 24, Dr. Colmenares
// This program utilizes pthreads for parallel computation on arrays.
// It initializes arrays (a, b, d, e) and creates threads for 
// subtraction and division operations. Each task is split into 
//manageable chunks, with four threads dedicated to each. 
// Subtraction computes the difference between corresponding elements
// of b and a, storing the result in c, while division computes the
// quotient of e and d, storing the result in f.
// Additionally, the program calculates
// the dot product of c and f using multiple threads, updating a
// shared global variable (dotproduct) with the computed values,
// protected by a mutex. Upon completion, it prints the summation
// of elements in c and f, along with the dot product, showcasing
// the efficiency of parallel processing in arrays.
// compile: gcc OlatundeA4Phreads.c -o para.out -lpthread
// execute: ./para.out

#include<stdio.h>
#include<pthread.h>

pthread_mutex_t mutex;

int a[30000];
int b[30000];
int d[30000];
int e[30000];


int c[30000];
int f[30000];                               

// Function prototypes
void * Sub_SectionOfArray(void * tid);
void * Div_SectionOfArray(void * tid);
long int sumArray( int * A, int size);
void * Dot_SectionOfArray(void * tid);


long int dotproduct = 0;

int main(int argc, char * argv[])
{
    // Initializing arrays a, b, d, e with predefined values
    for (int i = 0; i < 30000; i++)
    {
       a[i] = 1;
       b[i] = 2;
       d[i] = 5;     
       e[i] = 10;
    }

    pthread_t subThreads[4];
    pthread_t divThreads[4];
    pthread_t dotThreads[6];
    
    pthread_attr_t attr; 
    
    // Initialize mutex
    pthread_mutex_init(&mutex, NULL);
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

     for (int i =0; i < 4; i ++)
     {
        // Create threads for subtraction and division operations
        pthread_create(&subThreads[i], &attr, Sub_SectionOfArray, (void *)i);
        pthread_create(&divThreads[i], &attr, Div_SectionOfArray, (void *)i);

     }


     for (int i =0; i < 6; i++)
    {
        // Create threads for dot product calculation
       pthread_create(&dotThreads[i], &attr, Dot_SectionOfArray, (void *)i);
    }

    // Destroy thread attribute object
    pthread_attr_destroy(&attr);


    for (int i =0; i < 4; i ++)
     {
       // Join threads for subtraction and division operations
       pthread_join(subThreads[i], NULL);
       pthread_join(divThreads[i], NULL);
     }


     for (int i =0; i < 6; i++)
    {// Join threads for dot product calculation
      pthread_join(dotThreads[i], NULL);
    }
    
    // Print results
    printf("Summation of Ellements in array c is equal to: %ld \n", sumArray(c, 30000));
    printf("Summation of Ellements in array f is equal to: %ld \n", sumArray(f, 30000)); 

    printf("The dot procut between c and f is: %ld \n", dotproduct);

    // Cleanup
    pthread_mutex_destroy(&mutex);
    pthread_exit(NULL);
    return 0;
}

// Thread function to perform subtraction operation
void * Sub_SectionOfArray(void * tid)
{
   long thread_id = (long)tid;

   for (int i = thread_id *7500; i < (thread_id+1) * 7500 ; i++) 
   {
     c[i] = b[i] - a[i];
   }
   
   pthread_exit(NULL);
}

// Thread function to perform division operation
void * Div_SectionOfArray(void * tid)
{
  long thread_id = (long)tid;

  for (int i = thread_id * 7500; i < (thread_id + 1) * 7500; i++)
  {
    f[i] = e[i] /d [i];
  }
  pthread_exit(NULL);
}


// Function to calculate the summation of elements in an array
long int sumArray( int * A, int size)
{
   long int sum =0;

  for ( int i = 0; i < size;i++)
  {
    sum += A[i];
  }

  return sum;
}


// Thread function to calculate the dot product of arrays c and f
void * Dot_SectionOfArray(void * tid)
{
    long thread_id = (long)tid;
    long int loc_sum =0;
    for (int i = thread_id * 5000; i < (thread_id + 1) *5000; i++)
    {
        loc_sum +=c[i] * f[i];
        pthread_mutex_lock(&mutex);
        dotproduct += c[i] * f[i];
        pthread_mutex_unlock(&mutex);
    }
    
    printf("I am thread %d and my local sum is: %d \n", thread_id, loc_sum);
    pthread_exit(NULL);
}
