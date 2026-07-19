/*******************************************************************************
 *
 * File matrix_element_jackknife.c
 *
 * Calcolo dell'elemento di matrice e del relativo errore jackknife.
 *
 * Output:
 *      mat_el[t]        valore medio
 *      sigma_mat_el[t]  errore jackknife
 *
 ******************************************************************************/

#define MATRIX_ELEMENT_JACKKNIFE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.h"

void matrix_element_jackknife(double *corr_mean,
                              double *dE,
                              double *mat_el,
                              double *sigma_mat_el)
{
    int i, t;

    int nbin = Nconf / Bin;

    double *buffer   = malloc(N * sizeof(double));
    double *G_jack   = malloc(N * sizeof(double));
    double *M_jack   = malloc(nbin * N * sizeof(double));

    if (!buffer || !G_jack || !M_jack)
    {
        perror("malloc");

        free(buffer);
        free(G_jack);
        free(M_jack);
        return;
    }

    FILE *f = fopen("correlatorebinn.bin", "rb");

    if (!f)
    {
        perror("correlatorebinn.bin");

        free(buffer);
        free(G_jack);
        free(M_jack);
        return;
    }

/***************************************************************************/

    for (i = 0; i < nbin; i++)
    {
        if (fread(buffer, sizeof(double), N, f) != N)
        {
            fprintf(stderr, "Errore lettura bin %d\n", i);
            break;
        }

        /*
         * correlatore jackknife
         */
        for (t = 0; t < N; t++)
        {
            G_jack[t] =
                (nbin * corr_mean[t] - buffer[t])
                / (double)(nbin - 1);
        }


        /**************** elemento di matrice jackknife *******************/
         
         
        for (t = 0; t < N; t++)
        {
            double E = dE[t];

            if (!isfinite(E) || E <= 0.0)
            {
                M_jack[i*N + t] = NAN;
                continue;
            }

            double arg_exp = 0.5 * E * N * a;

            if (arg_exp > 700.0)
            {
                M_jack[i*N + t] = NAN;
                continue;
            }

            double num = G_jack[t] * exp(arg_exp);

            double arg_cosh = E * a * (0.5 * N - t);

            if (fabs(arg_cosh) > 700.0)
            {
                M_jack[i*N + t] = NAN;
                continue;
            }

            double denom = 2.0 * cosh(arg_cosh);

            if (!isfinite(denom) || denom == 0.0)
            {
                M_jack[i*N + t] = NAN;
                continue;
            }

            M_jack[i*N + t] = num / denom;
        }
    }

    fclose(f);

    /*
     * MEDIA JACKKNIFE + ERRORE
     */

    for (t = 0; t < N; t++)
    {
        double sum = 0.0;
        int nvalid = 0;

        for (i = 0; i < nbin; i++)
        {
            double x = M_jack[i*N + t];

            if (isfinite(x))
            {
                sum += x;
                nvalid++;
            }
        }

        if (nvalid == 0)
        {
            mat_el[t] = NAN;
            sigma_mat_el[t] = NAN;
            continue;
        }

        mat_el[t] = sum / (double)nvalid;

        double var = 0.0;

        for (i = 0; i < nbin; i++)
        {
            double x = M_jack[i*N + t];

            if (!isfinite(x))
                continue;

            double diff = x - mat_el[t];
            var += diff * diff;
        }

        sigma_mat_el[t] =
            sqrt((double)(nbin - 1)
                 / (double)nbin
                 * var);
    }

    free(buffer);
    free(G_jack);
    free(M_jack);
}