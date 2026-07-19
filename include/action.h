
#ifndef ACTION_H
#define ACTION_H

#ifndef SOA_C
extern double soa(void);
#endif

#ifndef DS_C
extern double ds(int i, double dxx);
#endif

#ifndef METROPOLIS_C
extern void metropolis(void);
#endif

#ifndef CORRELATORE_C
extern double *correlatore(void);
#endif

#ifndef AUTOCORRELATORE_C
extern void autocorrelatore(double corr[], int n, char *filename);
#endif

#ifndef CORRELATORE_MEAN_C
extern double *correlatore_mean(const char *filename_bin);
#endif

#ifndef CORRELATORE_SIGMA_C
extern double *correlatore_sigma(const char *filename_bin, double *media);
#endif


#ifndef ENERGY_JACKKNIFE_C
extern void effective_mass_jackknife(double *corr_mean,double *dE,double *sigma_dE);
#endif


#ifndef MATRIX_ELEMENT_C
extern void matrix_element_jackknife(double *corr_mean, double *dE, double *mat_el, double *sigma_mat_el);                           
#endif

#endif
