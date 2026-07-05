
/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define N 128
#define M 1.0
#define W 1.0
#define DELTA 1.0
#define a 0.5
#define Nconf 10000000
#define Bin 100
#define Nacor 100


#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

EXTERN double xx[N];

#undef EXTERN

#endif

