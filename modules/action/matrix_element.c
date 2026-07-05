/*******************************************************************************
*
* File matrix_element.c
* 
* 
*
*******************************************************************************/
#define MATRIX_ELEMENT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"


double *matrix_element(double *corr_mean, double *dE){

    int t;
    double *Me = malloc(N * sizeof(double));

    for (t = 0; t < N; t++)
    {
        if (dE[t] > 0.0)
        {
            double n = corr_mean[t]*exp(dE[t]*N*a/2.0);
            double d = 2*cosh( dE[t]*a* (N/2 - t) );

            Me[t] = n/d;
        }
        else       
            Me[t]=0;
        
    }
    
    return Me;
    
}