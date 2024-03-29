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
        timeDomain[i].real = 0.0;
        timeDomain[i].imag = 0.0;
    }
    
   calcCooleyTukey(FFT,timeDomain);

   printf("TOTAL PROCESSED SAMPLES: %d\n", N);
    printf("===========================================\n");
    printf("XR[0]: %f\t", FFT[0].real);
    printf("XI[0]: %f\n", FFT[0].imag);
    printf("===========================================\n");
    for(int i =1; i <= 10; i++)
    {
       printf("XR[%d]: %f\tXI[%d]: %f\n", i, FFT[i].real, i, FFT[i].imag);// %f\t\t", FFT[0].real);
    }
    printf("===========================================\n\n");
    
    return 0; 
}



//*******************************************************************
// Name::calcCooleyTukey()
// Parameters: 2 complexNum Pointers, and 1 int 
// A discussion of what the method/function does and required
// parameters as well as return value.
//CODE MUST HAVE COMMENTS, GOOD AND INFORMATIVE COMMENTS [-10 IF IGNORED]
//********************************************************************
void calcCooleyTukey(struct complexNum * FFT, struct complexNum * timeDomain)
{
  /*
   double even_real, even_imag;
   double odd_real, odd_imag;
   double theta;
   double twiddle_real, twiddle_imag;
   double sine, cosine;

    for (int k = 0; k < N/2; k++)
    {
        // initiallizing componenets
        even_real = even_imag = 0.0;
        odd_real = odd_imag = 0.0;
        theta = 0.0;
    
       
        twiddle_real = (cos ((2 * M_PI * k))/N);
        twiddle_imag = -(sin ((2 * M_PI * k)/N));

        for( int n = 0; n < N/2; n++)
        {
            theta = ((2 * M_PI * n * k)/ (N/2));
            
            cosine = cos(theta);
            sine =  sin (theta);
    
            even_real +=(timeDomain[(2 *n)].real * cosine ) + (timeDomain[(2* n)].imag * sine);
            even_imag += ((-(timeDomain[(2 * n) ].real)) * sine) + (timeDomain[2 * n].imag * cosine);

            odd_real +=(timeDomain[(2 *n + 1)].real * cosine ) + (timeDomain[(2* n+1)].imag * sine);
            odd_imag += ((-(timeDomain[(2 * n+1) ].real)) * sine) + (timeDomain[2 * n +1].imag * cosine);

           
        }   
        
        // compute coefficients for 1st half 
        FFT[k].real = even_real + (( twiddle_real * odd_real) - (twiddle_imag * odd_imag));
        FFT[k].imag = even_imag + (( twiddle_real * odd_imag) + (twiddle_imag * odd_real));


        // compute coefficients for 2st half 
        FFT[k + (N/2)].real = even_real - (( twiddle_real * odd_real) - (twiddle_imag * odd_imag));
        FFT[k + (N/2)].imag = even_imag - (( twiddle_real * odd_imag) + (twiddle_imag * odd_real));
        
    }
    */

    double real;
    double imag;
    double theta;
    float real2ndHalf = 0.0, imag2ndHalf = 0.0, theta2ndHalf = 0.0;

    for (int k = 0; k < (N/2); k++)
    {
        real = 0;
        imag = 0;
        theta = 0;
        real2ndHalf = 0.0, imag2ndHalf = 0.0, theta2ndHalf = 0.0;
        for(int n = 0; n < (N); n ++)
        {
            //theta = ((2 * M_PI) * 2 * n * k) / N;
            theta = ((2 * M_PI) * n * k) / N;
            
            real += (timeDomain[(n)].real * cos(theta) ) + (timeDomain[(n)].imag * sin(theta));
            imag += (-(timeDomain[(n) ].real) * sin(theta)) + (timeDomain[n ].imag * cos(theta));
            // dft of even indexed time domain
            /*real += (timeDomain[(2*n)].real * cos(theta) ) + (timeDomain[(2*n)].imag * sin(theta));
            imag += (-(timeDomain[(2*n)].real) * sin(theta)) + (timeDomain[(2*n)].imag * cos(theta));
            
            theta = ((2 * M_PI) * ((2 * n)+ 1 ) * k) / N;

            // dft of even indexed time domain
            real += (timeDomain[(2*n) + 1].real * cos(theta) ) + (timeDomain[(2*n) + 1].imag * sin(theta));
            imag += (-(timeDomain[(2*n) + 1].real) * sin(theta)) + (timeDomain[(2*n) + 1].imag * cos(theta));*/

            theta2ndHalf = (2 * M_PI * n * (k+(N/2)))/ N;

           real2ndHalf += (timeDomain[(n)].real * cos(theta2ndHalf)) +
                                (timeDomain[(n)].imag * sin(theta2ndHalf));

            imag2ndHalf +=  (-(timeDomain[(n)].real)* sin(theta2ndHalf)) + 
                                (timeDomain[(n)].imag * cos(theta2ndHalf)); 

        }


      FFT[k].real = real;
      FFT[k].imag = imag;

    
     FFT[k + (N/2)].real = real2ndHalf;
    FFT[k + (N/2)].imag = imag2ndHalf;
    }

}