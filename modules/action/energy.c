/*******************************************************************************
*
* File energy.c
* 
* 
*
*******************************************************************************/
#define ENERGY_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"


double *energy(double medie_G[]){

    int t;
    double *dE = malloc(N * sizeof(double));
    
  
       
    for (t = 1; t < N-1; t++)
        {
        double arg = (medie_G[t+1] + medie_G[t-1]) / (2.0 * medie_G[t]);
        
        if (arg >= 1.0)
            dE[t] = acosh(arg) / a;
        else
            dE[t] = 0.0;
        
        }
                
                    
    return dE;
}