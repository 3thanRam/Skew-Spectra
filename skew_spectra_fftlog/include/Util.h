#ifndef UTIL
#define UTIL


#include <complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_integration.h>
#include <time.h>


#include "../include/class_funcs.h"


int convind3to1(int n,int n1,int n2,struct pars_struct *pars);
int convind2to1(int n,int n1,struct pars_struct *pars);

void get_krange(double *krange,struct pars_struct *pars);

void fill_vector(gsl_vector_complex *vect,double k,double complex *fcoefs,double complex*biased_eta,int N);

void OperatorSn(int n,double complex nu1,double complex nu2,gsl_vector_complex *vect,struct pars_struct *pars);
void OperatorSn_unsym(int n,double complex nu1,gsl_vector_complex *vect,struct pars_struct *pars);

void get_matrices(struct pars_struct *pars);

double prodv1Mv2(gsl_vector_complex *v1,gsl_matrix_complex *M,gsl_vector_complex *v2,int N);
double sumv1M(gsl_vector_complex *v1,gsl_vector_complex *Munsym,int N);

void kmuvec(double k,double mu,gsl_vector_complex *vkmu);

void get_Mn(gsl_matrix_complex *Mn,gsl_vector_complex *vkmu,struct pars_struct *pars);
void get_UnsymMn(gsl_vector_complex *Mn_unsym,gsl_vector_complex *vkmu,struct pars_struct *pars);

double symintegrand(double x, void * params);
double unsymintegrand(double x, void * params);

void fill_Sn(double *krange,int N,struct pars_struct *pars,double **PSn);
void muk_integrate(double *krange,int N,double **PSn,struct pars_struct *pars);

#endif