/*******************************************************************************
*
* File correlatore_sigma.c
* 
* Calcolo deviazione standard del correlatore binnato
*
*******************************************************************************/
#define CORRELATORE_SIGMA_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

double *correlatore_sigma(const char *filename_bin, double *media)
{
    int t, i;
    int nconf = Nconf / Bin;

    FILE *f = fopen(filename_bin, "rb");
    if (!f) { perror("file bin"); return NULL; }

    double *buffer = malloc(N * sizeof(double));
    double *media2 = malloc(N * sizeof(double));
    double *sigma  = malloc(N * sizeof(double));  /* restituito al chiamante */

    if (!buffer || !media2 || !sigma)
    {
        perror("malloc");
        fclose(f);
        return NULL;
    }

    for (t = 0; t < N; t++) media2[t] = 0.0;

    for (i = 0; i < nconf; i++)
    {
        if (fread(buffer, sizeof(double), N, f) != N)
        {
            fprintf(stderr, "Errore lettura file binario\n");
            break;
        }
        for (t = 0; t < N; t++) media2[t] += buffer[t] * buffer[t];
    }

    for (t = 0; t < N; t++)
    {
        double var = media2[t]/nconf - media[t]*media[t];
        if (var < 0.0) var = 0.0;
        sigma[t] = sqrt(var / (nconf - 1));
    }

    free(buffer);
    free(media2);
    fclose(f);

    return sigma;   
}