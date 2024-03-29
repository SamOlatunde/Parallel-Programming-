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

#define N 16384