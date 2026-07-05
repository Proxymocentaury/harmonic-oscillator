/*******************************************************************************
*
* File autocorrelatore.c
* 
* Calcolo dell'autocorrelatore
*
*
* 
*******************************************************************************/

#define AUTOCORRELATORE_C

#include <stdlib.h>
#include <stdio.h>
#include "global.h"


void autocorrelatore(double corr[], int n, char *filename){

    double mean = 0;
    double sigma = 0;
    double sum;
    int i;
    int t;

    double *gamma;
    gamma = malloc( Nacor * sizeof(double));

    for (i = 0; i < n; i++) { mean += corr[i]; }     
    mean /= n;

    for (i = 0; i < n; i++) { sigma += corr[i]*corr[i]; }
    sigma = sigma/n - mean*mean;



    for (t = 0; t < Nacor; t++)
    {
        sum=0.0;
        for (i = 0; i < n-t; i++)
        {
            sum += corr[i]*corr[i+t] - mean*mean;
        }
        sum /= (n-t);
        gamma[t] = sum/sigma;
        
    }
    
    FILE *f;
    f = fopen(filename, "w");

    for (t = 0; t < Nacor; t++)
    {
        fprintf(f, "%d %.10f\n", t, gamma[t]);
    }

    fclose(f);
    free(gamma);        
}