/*******************************************************************************
*
* File azione_oa.c
* 
* Calcolo della differenza dell'azione di un oscillatore armonico
*
* The externally accessible functions are
*
* 
*******************************************************************************/

#define DS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"


double ds(int i, double dxx)
{
    int im = (i - 1 + N) % N;  
    int ip = (i + 1) % N;      

   return M/a * dxx * (2*xx[i] - xx[im] - xx[ip] + dxx) 
     + 0.5 * a * M * W * W * dxx * (2*xx[i] + dxx);
}

