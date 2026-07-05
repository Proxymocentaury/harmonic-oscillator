/*******************************************************************************
 *
 * File effective_mass_jackknife.c
 *
 * Calcolo della massa efficace e del relativo errore jackknife.
 *
 * Input:
 *      corr_mean[t]   correlatore medio (dalla media globale)
 *
 * Output:
 *      dE[t]          massa efficace calcolata dalla media globale
 *      sigma_dE[t]    errore jackknife
 *
 ******************************************************************************/
#define EFFECTIVE_MASS_JACKKNIFE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

void effective_mass_jackknife(double *corr_mean,
                               double *dE,
                               double *sigma_dE)
{
    int i, t;
    int nbin = Nconf / Bin;

    double *buffer  = malloc(N * sizeof(double));
    double *G_jack  = malloc(N * sizeof(double));
    double *acc     = malloc(N * sizeof(double));  /* accumulatore */

    if (!buffer || !G_jack || !acc)
    {
        perror("malloc");
        free(buffer); free(G_jack); free(acc);
        return;
    }

    /*
     * Step 1 — dE[t] dalla media globale corr_mean
     */
    dE[0]   = NAN;
    dE[N-1] = NAN;

    for (t = 1; t < N-1; t++)
    {
        double arg = (corr_mean[t+1] + corr_mean[t-1])
                   / (2.0 * corr_mean[t]);

        dE[t] = (arg >= 1.0 && isfinite(arg)) ? acosh(arg) / a : NAN;
    }

    /*
     * Step 2 — inizializza accumulatore
     */
    for (t = 0; t < N; t++) acc[t] = 0.0;

    /*
     * Step 3 — repliche jackknife
     */
    FILE *f = fopen("correlatorebinn.bin", "rb");
    if (!f)
    {
        perror("correlatorebinn.bin");
        free(buffer); free(G_jack); free(acc);
        return;
    }

    for (i = 0; i < nbin; i++)
    {
        if (fread(buffer, sizeof(double), N, f) != N)
        {
            fprintf(stderr, "Errore lettura bin %d\n", i);
            break;
        }

        /* correlatore jackknife G_bar^(k)(t) */
        for (t = 0; t < N; t++)
            G_jack[t] = (nbin * corr_mean[t] - buffer[t])
                      / (double)(nbin - 1);

        /* energia jackknife E^(k)(t) e accumulo */
        for (t = 1; t < N-1; t++)
        {
            if (!isfinite(dE[t])) continue;

            double arg = (G_jack[t+1] + G_jack[t-1])
                       / (2.0 * G_jack[t]);

            if (arg < 1.0 || !isfinite(arg)) continue;

            double E_jack = acosh(arg) / a;
            double diff   = E_jack - dE[t];
            acc[t]       += diff * diff;
        }
    }

    fclose(f);

    /*
     * Step 4 — errore jackknife
     */
    sigma_dE[0]   = NAN;
    sigma_dE[N-1] = NAN;

    for (t = 1; t < N-1; t++)
    {
        if (!isfinite(dE[t]))
        {
            sigma_dE[t] = NAN;
            continue;
        }
        sigma_dE[t] = sqrt((double)(nbin - 1) / (double)nbin * acc[t]);
    }

    free(buffer);
    free(G_jack);
    free(acc);
}