/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define N 1024
#define M 1.0
#define W 1.0
#define DELTA 0.5
#define a 0.0625
#define Nconf 50000000
#define Bin 5000
#define Nacor 10


#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

EXTERN double xx[N];

#undef EXTERN

#endif
