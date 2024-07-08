#ifndef PLOTTING
#define PLOTTING

#include <stdio.h>
#include <math.h>
#include "../include/class_funcs.h"


void sendPdata(FILE *gnuplotPipe,double *kdata,double *P,int N);
void plot_data(struct pars_struct *pars,int N,double *karray, double **PSn);

#endif