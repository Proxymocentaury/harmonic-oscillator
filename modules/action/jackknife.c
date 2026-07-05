/*******************************************************************************
*
* File energy_jackknife.c
* 
* Calcolo dell'errore sul Delta E tramite il metodo Jackknife
*
* The externally accessible functions are
*
*   double *energy_jackknife(double *corr_mean, double *dE)
*     Calcola l'errore jackknife su Delta E leggendo la matrice
*     binnata da correlatorebinn.bin
*     Restituisce il vettore sigma_dE[N] — il chiamante fa free()
*
*******************************************************************************/
#define ENERGY_JACKkNIFE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

double *jackknife(double *corr_mean, double *dE)
{
    int k, t;
    int n = Nconf / Bin;

    double *buffer  = malloc(N * sizeof(double));
    double *G_jack  = malloc(N * sizeof(double));
    double *dE_jack = malloc(N * sizeof(double));
    double *acc     = malloc(N * sizeof(double));
    double *sigma   = malloc(N * sizeof(double));  /* restituito al chiamante */

    if (!buffer || !G_jack || !dE_jack || !acc || !sigma)
    {
        perror("malloc");
        return NULL;
    }

    /* inizializzazione accumulatore */
    for (t = 0; t < N; t++) acc[t] = 0.0;

    FILE *f = fopen("correlatorebinn.bin", "rb");
    if (!f) { perror("correlatorebinn.bin"); return NULL; }

    for (k = 0; k < n; k++)
    {
        /* leggi riga k */
        if (fread(buffer, sizeof(double), N, f) != N)
        {
            fprintf(stderr, "Errore lettura bin k=%d\n", k);
            break;
        }

        /* media jackknife G_bar^(k)(t) */
        for (t = 0; t < N; t++)
            G_jack[t] = (n * corr_mean[t] - buffer[t]) / (double)(n - 1);

        /* energia jackknife DeltaE^(k)(t) */
        dE_jack[0]   = 0.0;
        dE_jack[N-1] = 0.0;

        for (t = 1; t < N-1; t++)
        {
            double arg = (G_jack[t+1] + G_jack[t-1]) / (2.0 * G_jack[t]);
            dE_jack[t] = (arg >= 1.0) ? acosh(arg) : 0.0;
        }

        /* accumula (dE_jack[t] - dE[t])^2 */
        for (t = 1; t < N-1; t++)
            acc[t] += (dE_jack[t] - dE[t]) * (dE_jack[t] - dE[t]);
    }

    fclose(f);

    /* calcolo sigma jackknife */
    sigma[0]   = 0.0;
    sigma[N-1] = 0.0;

    for (t = 1; t < N-1; t++)
        sigma[t] = sqrt((double)(n-1) / (double)n * acc[t]);

    free(buffer);
    free(G_jack);
    free(dE_jack);
    free(acc);

    return sigma;
}