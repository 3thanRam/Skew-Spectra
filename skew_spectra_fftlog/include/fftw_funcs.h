#ifndef FFTW_FUNCS
#define FFTW_FUNCS

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include "../include/class_funcs.h"


void fill_Pb(fftw_complex *P_biased,double *krange,struct pars_struct *pars);
void get_fouriercoefs(double *krange,struct pars_struct *pars);
void get_eta(struct pars_struct *pars);
#endif