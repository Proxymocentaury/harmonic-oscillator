/*******************************************************************************
*
* File correlatore_mean.c
* 
* Calcolo media correlatore binnato
*
*******************************************************************************/
#define CORRELATORE_MEAN_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

double *correlatore_mean(const char *filename_bin)
{
    int t, i;
    int nconf = Nconf / Bin;

    FILE *f = fopen(filename_bin, "rb");
    if (!f) { perror("file bin"); return NULL; }

    double *buffer = malloc(N * sizeof(double));
    double *media  = malloc(N * sizeof(double));  
    if (!buffer || !media)
    {
        perror("malloc");
        fclose(f);
        return NULL;
    }

    for (t = 0; t < N; t++) media[t] = 0.0;

    for (i = 0; i < nconf; i++)
    {
        if (fread(buffer, sizeof(double), N, f) != N)
        {
            fprintf(stderr, "Errore lettura file binario\n");
            break;
        }
        for (t = 0; t < N; t++) media[t] += buffer[t];
    }

    for (t = 0; t < N; t++) media[t] /= nconf;

    free(buffer);
    fclose(f);

    return media;   
}