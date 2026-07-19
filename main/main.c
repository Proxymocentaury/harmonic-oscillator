/*******************************************************************************
*
* File main.c
* 
* 
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "action.h"
#include "random.h"

int main() {

    int i,j,k, c=1;
    double action;
    double DE=0.0, sigma_DE=0.0;
    double Mel=0.0, sigma_Mel=0.0;

    double *corr;
    double *corrbin;

    double *coll_corr;
    double *coll_corrbin;

    double *corr_mean;
    double *corr_sigma;

    double *dE;
    double *sigma_dE;
    double *mat_el;
    double *sigma_mat_el;

    
    corrbin = malloc( N * sizeof(double));
    coll_corr = malloc( Nconf * sizeof(double));
    coll_corrbin = malloc( (Nconf/Bin) * sizeof(double));

    rlxd_init(1, (int)time(NULL) );


/*********************** TERMALIZZAZIONE **********************************************/
    FILE *f1; 

    for (i = 0; i < N; i++) xx[i] = 0.0;

    f1 = fopen("termalizzazione.txt", "w");
    
    if (!f1) {
        perror("termalizzazione.txt");
        return 1;
    }

    action = soa();

    fprintf(f1, "%d %.10f\n", 0, action);

    for (i = 1; i < 5000; i++) {
        metropolis();
        action = soa();
        fprintf(f1, "%d %.10f\n", i, action);
    }

    fclose(f1);
    
/*********************** MONTECARLO **********************************************/
            
    i=0;

    FILE *fbin = fopen("correlatorebinn.bin", "wb");

    while (i < Nconf/Bin)
    {

        for ( k = 0; k < N; k++) { corrbin[k]=0.0; }

        for ( j = 0; j < Bin; j++)
        {            
            metropolis();

            corr = correlatore();       
            for ( k = 0; k < N; k++) {corrbin[k] += corr[k]/Bin;} 
            
            coll_corr[i*Bin + j] = corr[c];
            free(corr);
        }

        fwrite(corrbin, sizeof(double), N, fbin);

        coll_corrbin[i] = corrbin[c];
        i++;        

    }/*end while*/
    
    fclose(fbin);

/********************** STATISTIC ***************************************************************/

    autocorrelatore(coll_corr, Nconf, "autocorrelatore.txt");
    autocorrelatore(coll_corrbin, Nconf/Bin, "autocorrelatorebin.txt");


    
    corr_mean = correlatore_mean("correlatorebinn.bin");
    corr_sigma = correlatore_sigma("correlatorebinn.bin", corr_mean);

    FILE *fstats = fopen("medie_sigma_corr.txt", "w");

    if (fstats)
    {
        for (i = 0; i < N; i++)
            fprintf(fstats, "%d %.10f %.10f\n", i, corr_mean[i], corr_sigma[i]);
        fclose(fstats);
    }
    

/*********************** ESTRAZIONE DELLE OSSERVABILI **********************************************/


    dE = malloc(N * sizeof(double));
    sigma_dE = malloc(N * sizeof(double));

    effective_mass_jackknife(corr_mean, dE, sigma_dE);

    mat_el = malloc(N * sizeof(double));
    sigma_mat_el = malloc(N * sizeof(double));
    
    matrix_element_jackknife(corr_mean, dE, mat_el, sigma_mat_el);
                                           
    
    FILE *fdE = fopen("effective_mass.txt", "w");

    if (!fdE) {
    perror("effective_mass.txt");
    return 1;
    }

    for (i = 1; i < N-1; i++) {
        fprintf(fdE, "%d %.10f %.10f\n",
            i, dE[i], sigma_dE[i]);
    }
    fclose(fdE);


    FILE *fM = fopen("matrix_element.txt", "w");

    if (!fM) {
    perror("matrix_element.txt");
    return 1;
    }

    for (i = 1; i < N-1; i++) {
        fprintf(fM, "%d %.10f %.10f\n", i, mat_el[i], sigma_mat_el[i]);         
    }
    fclose(fM);


    for (i = 1; i < N/2; i++)
    {
        if (dE[i] > 0.0 && sigma_dE[i] > 0.0)
        {
            double err_rel = sigma_dE[i] / dE[i];
            if (err_rel <= 0.1)
            {
                DE       += dE[i] / (sigma_dE[i] * sigma_dE[i]);
                sigma_DE += 1.0   / (sigma_dE[i] * sigma_dE[i]);
            }
        }
    }

    DE/= sigma_DE;
    sigma_DE  = 1.0 / sqrt(sigma_DE);

    printf("Delta E = %.6f +- %.6f\n", DE, sigma_DE);


      for (i = 1; i < N/2; i++)
    {
        if (mat_el[i] > 0.0 && sigma_mat_el[i] > 0.0)
        {
            double err_rel = sigma_mat_el[i] / mat_el[i];
            if (err_rel <= 0.1)
            {
                Mel       += mat_el[i] / (sigma_mat_el[i] * sigma_mat_el[i]);
                sigma_Mel += 1.0   / (sigma_mat_el[i] * sigma_mat_el[i]);
            }
        }
    }

    Mel/= sigma_Mel;
    sigma_Mel  = 1.0 / sqrt(sigma_Mel);

    printf("Matrix element = %.6f +- %.6f\n", Mel, sigma_Mel);



    free(corrbin);
    free(coll_corr);
    free(coll_corrbin);
    free(corr_mean);
    free(corr_sigma);
    free(dE);
    free(sigma_dE);
    free(mat_el);
    free(sigma_mat_el);
           
    return 0;
}
