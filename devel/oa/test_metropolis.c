/*******************************************************************************
*
* File test_acceptance.c
*
* Calcola il tasso di accettazione del Metropolis
* per verificare che DELTA sia ottimale (~50%)
*
*******************************************************************************/
#define MAIN_PROGRAM
#define TEST_ACCEPTANCE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "action.h"
#include "random.h"

int main(void)
{
    int i, j;
    int Nsweep = 1000;
    int accettate = 0;
    int totale = 0;
    double random_v[2*N];
    double dxx, delta_S;

    rlxd_init(1, 123);

    for (i = 0; i < N; i++) xx[i] = 0.0;

    /* termalizzazione */
    for (i = 0; i < 500; i++)
        metropolis();

    /* misura del tasso di accettazione */
    for (i = 0; i < Nsweep; i++)
    {
        ranlxd(random_v, 2*N);

        for (j = 0; j < N; j++)
        {
            dxx     = 2.0 * DELTA * (random_v[j] - 0.5);
            delta_S = ds(j, dxx);

            if (delta_S < 0.0 || exp(-delta_S) > random_v[j+N])
            {
                xx[j] += dxx;
                accettate++;
            }
            totale++;
        }
    }

    printf("N     = %d\n",   N);
    printf("a     = %.6f\n", (double)a);
    printf("DELTA = %.6f\n", (double)DELTA);
    printf("Tasso di accettazione = %.2f%%\n",
           100.0 * accettate / totale);
    printf("Ottimale: 50%%\n");

    return 0;
}