#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#include <stdlib.h> // for atof and atoi
#include <string.h> // for strcmp
#include "../include/ini.h"

#ifndef PARSTRUCT
#define PARSTRUCT
struct pars_struct{
    int N;
    double** ranges;
};
#endif


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifndef RANDOMRANGE
#define RANDOMRANGE

void shuffle(double **array, int n);
int Numbdivs(double **ranges,int N,double *nulist,double dnu);
void reduceincrements(double **arr, int size);
int numbnonzero(double *arr, int size);
void appenddouble(double **array,int *size,double val);
double minincrement(double **ranges,int N);

void genconfig(double **config,int N,int Nc,double min,double max,double dnu);
void findbest(double **ranges,int N,double **configbest);
void printing(double **ranges,int N,double *nuconfig,double dnu);
void mathematicaOutput(double **ranges,int N,double *nuconfig,double dnu);
void findnuvals();//double **ranges,int N

#endif