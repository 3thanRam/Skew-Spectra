#include "../include/Util.h"

int convind3to1(int n,int n1,int n2,struct pars_struct *pars){
    int N=pars->N;
    return N*(n*N+n1)+n2;
}
int convind2to1(int n,int n1,struct pars_struct *pars){
    int N=pars->N;
    return n*N+n1;
}


void get_krange(double *krange,struct pars_struct *pars){
    //double delta=(kmax-kmin)/N;
    double kmax=pars->kmax;
    double kmin=pars->kmin;
    double N=pars->N;

    double delta=log(kmax/kmin)/N;
    for(int i=0;i<N;i++){
    krange[i]=kmin*exp(i*delta);
    //krange[i]=kmin+i*delta;
    }
}




void fill_vector(gsl_vector_complex *vect,double k,double complex *fcoefs,double complex*biased_eta,int N){
    for (int i = 0; i < N; i++){
        double complex z=fcoefs[i]*cpow(k,biased_eta[i]);
        gsl_vector_complex_set(vect, i,gsl_complex_rect(creal(z), cimag(z)) );
    }
}

void OperatorSn(int n,double complex nu1,double complex nu2,gsl_vector_complex *vect,struct pars_struct *pars){
    double complex result;
    result= 1+_Complex_I*1;

    for (int i = 0; i < 8; i++){

    switch (n) {
      case 1:
        result= 1+_Complex_I*1;
        break;
      case 2:
        result= 1+_Complex_I*1;
        break;
      default:
        result= 1+_Complex_I*1;
    }
    gsl_vector_complex_set(vect, i,gsl_complex_rect(creal(result), cimag(result)) );
    }
    
}
void OperatorSn_unsym(int n,double complex nu1,gsl_vector_complex *vect,struct pars_struct *pars){
    double complex result;
    result= 1+_Complex_I*1;

    for (int i = 0; i < 8; i++){

    switch (n) {
      case 1:
        result= 1+_Complex_I*1;
        break;
      case 2:
        result= 1+_Complex_I*1;
        break;
      default:
        result= 1+_Complex_I*1;
    }
    gsl_vector_complex_set(vect, i,gsl_complex_rect(creal(result), cimag(result)) );
    }
    
}

void get_matrices(struct pars_struct *pars){
    int Ns=pars->Ns;
    int N=pars->N;
    double complex *biased_eta=pars->biased_eta;

    for (int n = 0; n < Ns; n++) {
        for (int i = 0; i < N; i++){
            double complex nu1=-0.5*biased_eta[i] ;
            gsl_vector_complex *vect_unsym = gsl_vector_complex_alloc(8); 
            OperatorSn_unsym(n, nu1, vect_unsym, pars);
            pars->M_vv[convind2to1(n, i, pars)] = vect_unsym;

            for (int j = 0; j < N; j++){
                double complex nu2=-0.5*biased_eta[j] ;
                gsl_vector_complex *vect = gsl_vector_complex_alloc(8); 
                OperatorSn(n, nu1, nu2, vect, pars);
                pars->M_vv[convind3to1(n, i, j, pars)] = vect;
            
    }}}
}


double prodv1Mv2(gsl_vector_complex *v1,gsl_matrix_complex *M,gsl_vector_complex *v2,int N){

    gsl_vector_complex * utemp = gsl_vector_complex_alloc (N);// Temp vector to store M_22 * v_m2
    double alpha[2]={1.0,0.0};
    double beta[2]={0.0,0.0};
    double result[2];
    ////Perform the matrix-vector product utemp=M_22 * v_m2 
    cblas_zgemv(CblasRowMajor, CblasNoTrans, N, N, alpha, M->data, N, v2->data, 1, beta, utemp->data, 1);
    // Compute the dot product v_m1^T * utemp
    cblas_zdotu_sub(N, v1->data, 1, utemp->data, 1,result);
    // Free the temporary vector
    gsl_vector_complex_free(utemp);
    return result[0]; //real part
}

double sumv1M(gsl_vector_complex *v1,gsl_vector_complex *Munsym,int N){

    double result[2];
    // Compute the dot product v_m1^T * M13
    cblas_zdotu_sub(N, v1->data, 1, Munsym->data, 1,result);

    return result[0]; //real part sqrt(result[0]*result[0]+result[1]*result[1])
}

void kmuvec(double k,double mu,gsl_vector_complex *vkmu){
    double k1 = k;
    double k2 = k * k;
    double k3 = k2 * k;
    double k4 = k3 * k;
    
    double mu0 = 1.0;
    double mu1 = mu;
    double mu2 = mu * mu;
    double mu3 = mu2 * mu;
    double mu4 = mu3 * mu;

    gsl_vector_complex_set(vkmu, 0, gsl_complex_rect(k1 * mu1, 0.0));
    gsl_vector_complex_set(vkmu, 1, gsl_complex_rect(k2 * mu0, 0.0));
    gsl_vector_complex_set(vkmu, 2, gsl_complex_rect(k2 * mu2, 0.0));
    gsl_vector_complex_set(vkmu, 3, gsl_complex_rect(k3 * mu1, 0.0));
    gsl_vector_complex_set(vkmu, 4, gsl_complex_rect(k3 * mu3, 0.0));
    gsl_vector_complex_set(vkmu, 5, gsl_complex_rect(k4 * mu0, 0.0));
    gsl_vector_complex_set(vkmu, 6, gsl_complex_rect(k4 * mu2, 0.0));
    gsl_vector_complex_set(vkmu, 7, gsl_complex_rect(k4 * mu4, 0.0));

}

void get_Mn(gsl_matrix_complex *Mn,gsl_vector_complex *vkmu,struct pars_struct *pars){
    int N=pars->N;
    int n=pars->indexn;
    gsl_vector_complex **M_vv=pars->M_vv;
    double result[2];

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            gsl_vector_complex *M_v=M_vv[convind3to1(n, i, j,pars)];

            cblas_zdotu_sub(N, M_v->data, 1, vkmu->data, 1,result);
            gsl_matrix_complex_set(Mn, i, j,gsl_complex_rect(result[0], result[1]));
    }}
        
}

void get_UnsymMn(gsl_vector_complex *Mn_unsym,gsl_vector_complex *vkmu,struct pars_struct *pars){
    int N=pars->N;
    int n=pars->indexn;
    gsl_vector_complex **M_vvunsym=pars->M_vvunsym;
    double result[2];

    for (int i = 0; i < N; i++){
        gsl_vector_complex *M_v_unsym=M_vvunsym[convind2to1(n, i,pars)];

        cblas_zdotu_sub(N, M_v_unsym->data, 1, vkmu->data, 1,result);
        gsl_vector_complex_set(Mn_unsym, i, gsl_complex_rect(result[0], result[1]));//gsl_complex_rect(result[0], 0.0));
    }
}

double symintegrand(double x, void * params) {
    
    struct pars_struct *pars = (struct pars_struct *)params; 

    double muk=x*2.-1;
    
    int N=pars->N;
    double k=pars->k;

    gsl_matrix_complex *Mn=gsl_matrix_complex_alloc(N, N);
    gsl_vector_complex *v1=pars->v1;
    gsl_vector_complex *v2=pars->v2;
    gsl_vector_complex *vkmu=pars->vkmu;

    kmuvec(k,muk,vkmu);
    get_Mn(Mn,vkmu,pars);    

    double Sym=2*prodv1Mv2(v1,Mn,v2,N); //*2 from Jacobian
    gsl_matrix_complex_free(Mn);
    return Sym;
}
double unsymintegrand(double x, void * params) {
    
    struct pars_struct *pars = (struct pars_struct *)params; 

    double muk=x*2.-1;
    
    int N=pars->N;
    double k=pars->k;

    gsl_vector_complex *Mnunsym=gsl_vector_complex_alloc(N);
    
    gsl_vector_complex *v1=pars->v1;
    gsl_vector_complex *vkmu=pars->vkmu;

    kmuvec(k,muk,vkmu);
    get_UnsymMn(Mnunsym,vkmu,pars);  

    double UnSym=sumv1M(v1,Mnunsym,N);

    gsl_vector_complex_free(Mnunsym);
    return UnSym;
}


void fill_Sn(double *krange,int N,struct pars_struct *pars,double **PSn){
    double step=2*M_PI/N;

    for(int i=0;i<N;i++){
        PSn[0][i]=cos(i*step);
        PSn[1][i]=sin(i*step);
    }
}


void muk_integrate(double *krange,int N,double **PSn,struct pars_struct *pars){//struct pars_struct *pars

    int Ns=pars->Ns;
    double complex *fcoefs=pars->fcoefs;
    double complex *biased_eta=pars->biased_eta;
    double k;

    double resultsym,resultunsym,error;
    gsl_function F;
    
    F.params = pars;
    gsl_integration_workspace * w= gsl_integration_workspace_alloc(100000);

    for(int i=0;i<N;i++){
        k=krange[i];
        pars->k=k;
        fill_vector(pars->v1,k,fcoefs, biased_eta, N);
        fill_vector(pars->v2,k,fcoefs, biased_eta, N);
        for(int n=0;n<Ns;n++){
        pars->indexn=n;

        F.function = &symintegrand;
        gsl_integration_qags(&F, 0., 1., 0., 1e-7, 100000,w, &resultsym, &error);
        F.function = &unsymintegrand;
        gsl_integration_qags(&F, 0., 1., 0., 1e-7, 100000,w, &resultunsym, &error);
        
        PSn[n][i]=pow(k,3)*(resultsym+PL(k)*resultunsym);
    }
    }
    
    gsl_integration_workspace_free(w);
}


