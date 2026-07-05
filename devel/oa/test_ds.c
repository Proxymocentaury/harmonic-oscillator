/*******************************************************************************
*
* File test_action.c
*
* controlla la differenza dell'azione
*
*******************************************************************************/

#define MAIN_PROGRAM
#define TEST_DELTA_ACTION_C

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "action.h"
#include "random.h"

int main() {

    int i;
    int k=0;
    double r[1];
    double delta=3.0;
    double dxx;
    double S_old, S_new;
    double deltas_manual, deltas;


    rlxd_init(1, 12345);  /* inizializza generatore random in memoria */
    ranlxd(r, 1); /* genera numero casuale in [0,1) */
    
    dxx   = 2.0 * delta * (r[0] - 0.5);

    
    for (i = 0; i < N; i = i + 1) {
        xx[i] = i;
    }


    while (k<N)
    {
        S_old = soa();


        /* perturbazione e calcolo dell'azione nuova */
        xx[k] = xx[k] + dxx;
        
        S_new = soa();
        deltas_manual = S_new - S_old;

        /* ripristino della configurazione originale */
        xx[k] = xx[k] - dxx;

        deltas = ds(k, dxx);

        /* stampa dei risultati */
        /*printf("Delta brutale = %.6f\n", deltas_manual);
        printf("Delta fast    = %.6f\n", deltas);*/

        printf("Differenza    = %.6e\n", deltas_manual - deltas);

        k +=1;
    }
    
    
    return 0;
}
