/*******************************************************************************
*
* File azione_oa.c
* 
* Calcolo dell'azione di un oscillatore armonico
*
* The externally accessible functions are
*
*   double azione_oa(void)
*     calcola l'azione per un oscillatore armonico
* 
* 
*  il mod  (%) serve per i = N-1, i+1 = N ma N mod(N) = 0 
*******************************************************************************/

#define SOA_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

double soa(void) {
    double sum = 0.0;
    double dx;
    int i = 0; 

    for (i = 0; i < N; i++)
    {
        dx = xx[(i+1)%N] - xx[i];
        sum+= 0.5 * M / a * dx * dx + 0.5 * a * M * W * W * xx[i] * xx[i];           
    }
    
    return sum;
}
