/*******************************************************************************
*
* File metropolis.c
*
* calcola nuovi cammini di feynmann tramite l'algoritmo di metropolis
*
*******************************************************************************/

#define METROPOLIS_C

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "action.h"
#include "random.h"

void metropolis(void) {

    int i;
    double delta_S;
    double dxx;

    double random_v[2*N];

    ranlxd(random_v, 2*N); /* genera numeri casuali in [0,1) */    

    for ( i = 0; i < N; i++)
    {
        dxx= 2.0 * DELTA * (random_v[i] - 0.5);

        delta_S = ds(i, dxx);        

        if (delta_S < 0.0 || exp(-delta_S) > random_v[i + N])
        {
            xx[i] += dxx;   
        }
                
    }
    
}

