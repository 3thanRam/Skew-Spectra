#include "../include/fftw_funcs.h"


void get_eta(struct pars_struct *pars){
    double kmin=pars->kmin;
    double kmax=pars->kmax;
    double nu=pars->nu;
    int N=pars->N;
    double complex biased_eta[N];

    double fact=2. * M_PI *((double)N - 1.)/ ((double)N * log(kmax/kmin));
    for (int m=0; m<N; m++) {
        biased_eta[m] =nu+_Complex_I*  ((double)(m) - (double)N/2.)* fact;
    }

    pars->biased_eta=biased_eta;
}

void fill_Pb(fftw_complex *P_biased,double *krange,struct pars_struct *pars){
    double nu=pars->nu;
    int N=pars->N;

    double delta=log(krange[N-1]/krange[0])/N;
    for (int i = 0; i < N; i++){
        double k=krange[i];
        P_biased[i]=PL(k)*exp(-nu*i*delta);
    }
}

void get_fouriercoefs(double *krange,struct pars_struct *pars){
    
    double kmin=pars->kmin;
    double complex *biased_eta=pars->biased_eta;
    int N=pars->N;


    fftw_complex *P_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fill_Pb(P_b,krange,pars);


    

    double complex fcoefs[N]; 
    fftw_complex *out;
    fftw_plan p;

    
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N, P_b, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); 


    for (int i=0; i < N; i++){
        if (i < N/2){
            fcoefs[i]= cpow(kmin,-biased_eta[i]) *(creal(out[N/2 - i]) - _Complex_I * cimag(out[N/2 - i]))/N;
        }
        else{
            fcoefs[i]= cpow(kmin,-biased_eta[i]) *(creal(out[i - N/2]) + _Complex_I * cimag(out[i - N/2]))/N;
        }
    }
    fcoefs[0]= fcoefs[0]/2.;
    fcoefs[N]= fcoefs[N]/2.;
    fftw_destroy_plan(p);
    fftw_free(P_b); fftw_free(out);
    pars->fcoefs=fcoefs;
}