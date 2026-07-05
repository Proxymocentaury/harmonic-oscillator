/*******************************************************************************
*
* File correlatore.c
* 
* Calcolo del correlatore
*
*
* 
*******************************************************************************/

#define CORRELATORE_C

#include <stdlib.h>
#include <stdio.h>
#include "global.h"


double *correlatore(void){
    
    double *corr;
    int t;
    int k;
    corr = malloc(N * sizeof(double));

    for ( t = 0; t < N; t++)
    {
        corr[t]=0;
    }
    

    for ( t = 0; t < N; t++)
    {
        for (k = 0; k < N; k++)
        {
            corr[t] += xx[k]*xx[(k+t)%N];
        }

        corr[t]= corr[t]/N;
        
    }

    return corr;

}
