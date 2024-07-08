#include <math.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <complex.h>

#include <stdlib.h> // for atof and atoi
#include <string.h> // for strcmp
#include "../include/ini.h"
#ifndef PARSTRUCT
#define PARSTRUCT
struct pars_struct{
    int N;
    int Ns;
    int indexn;
    double k;
    double kmax,kmin;
    double qmax,qmin;
    double logqmax,logqmin;
    double kf;
    double nu;
    double mom_units;

    double h;
    double z;
    const char *filePath;
    
    double b1,b2,bg2,f;
    double R;
    
    double complex *biased_eta;
    double complex *fcoefs;

    gsl_vector_complex *v1;
    gsl_vector_complex *v2;
    gsl_vector_complex *vkmu;

    double (*Opn)(int i,struct pars_struct *s);
    double (*W)(struct pars_struct* s,double k);

    gsl_vector_complex **M_vv ;
    gsl_vector_complex **M_vvunsym ;

};
#endif
 
#ifndef CLASS_FUNCS
#define CLASS_FUNCS
void init_cosmoparams(struct pars_struct *pars);
void free_vector_complex_arrays(gsl_vector_complex **vectors, size_t num_vectors);
double Wgauss(struct pars_struct* s,double k);

void init_class(struct pars_struct *pars);
void free_vars(struct pars_struct *pars);

void getkmax(double kmin,double *kmax,double nu);
void read_classdata(double *class_kdata,double *class_pkdata,int size,const char *filePath);
void init_PL(double *class_kdata,double *class_pkdata,int size);
double PL(double k);
void free_PL();

int numb_lines(const char *filePath);
#endif