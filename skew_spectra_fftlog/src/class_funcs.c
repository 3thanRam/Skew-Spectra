#include "../include/class_funcs.h"

static gsl_interp_accel *acc = NULL;
static gsl_spline *spline = NULL;
static int initialized = 0;


static int config_handler(void* user, const char* section, const char* name, const char* value) {
    struct pars_struct* config = (struct pars_struct*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("General", "filePath")) {
        config->filePath = strdup(value); // Remember to free it later
    } else if (MATCH("General", "h")) {
        config->h = atof(value);
    } else if (MATCH("General", "z")) {
        config->z = atof(value);
    } else if (MATCH("General", "nu")) {
        config->nu = atof(value);
    } else if (MATCH("General", "N")) {
        config->N = atoi(value);
    } else if (MATCH("General", "Ns")) {
        config->Ns = atoi(value);
    } else if (MATCH("General", "kmin")) {
        config->kmin = atof(value);
    } else if (MATCH("General", "kmax")) {
        config->kmax = atof(value);
    } else if (MATCH("General", "mom_units")) {
        config->mom_units = atof(value);
    } else if (MATCH("General", "qmin")) {
        config->qmin = config->mom_units*atof(value);
    } else if (MATCH("General", "qmax")) {
        config->qmax = config->mom_units*atof(value);
    } else if (MATCH("General", "kf")) {
        config->kf = config->mom_units*atof(value);
    } else if (MATCH("General", "R")) {
        config->R = atof(value)/config->mom_units;
    } else {
        return 0;  // Unknown section/name, error
    }
    return 1;
}



void free_vector_complex_arrays(gsl_vector_complex **vectors, size_t num_vectors) {
    if (vectors != NULL) {
        for (size_t i = 0; i < num_vectors; i++) {
            if (vectors[i] != NULL) {
                gsl_vector_complex_free(vectors[i]);
            }
        }
        free(vectors);
    }
}

double Wgauss(struct pars_struct* s,double k){
    return exp(-0.5*pow(k*(s->R),2.));
}
void init_cosmoparams(struct pars_struct *pars){
    if (!pars){
        fprintf(stderr, "Main:main:pars->params malloc failure\n");
        exit(EXIT_FAILURE);
    }
    if (ini_parse("cosmoconfig.ini", config_handler, pars) < 0) {
        printf("Can't load 'cosmoconfig.ini'\n");
        exit(EXIT_FAILURE);
    }
    int N=pars->N;
    int Ns=pars->Ns;
    
    pars->logqmax=log(pars->qmax);
    pars->logqmin=log(pars->qmin);    
    pars->W=Wgauss;

    gsl_vector_complex *v1 = gsl_vector_complex_alloc(N);    //N length vector for summing over m1 index
    gsl_vector_complex *v2 = gsl_vector_complex_alloc(N);    //N length vector for summing over m2 index
    gsl_vector_complex *vkmu=gsl_vector_complex_alloc(8);
    pars->v1=v1;
    pars->v2=v2;  
    pars->vkmu=vkmu;  

    int num_vectors=Ns*N*N;
    int num_vectors_unsym=Ns*N;
    gsl_vector_complex **vectors = malloc(num_vectors * sizeof(gsl_vector *));
    gsl_vector_complex **vectors_unsym = malloc(num_vectors_unsym * sizeof(gsl_vector *));
    for (size_t i = 0; i < num_vectors; i++) {
        vectors[i] = gsl_vector_complex_alloc(8);
    }
    for (size_t i = 0; i < num_vectors_unsym; i++) {
        vectors_unsym[i] = gsl_vector_complex_alloc(8);
    }
    pars->M_vv=vectors;
    pars->M_vvunsym=vectors_unsym;
}

void init_class(struct pars_struct *pars){

    double kmin=pars->kmin;
    double kmax=pars->kmax;
    double nu=pars->nu;
    const char *filePath=pars->filePath;

    int class_dsize=numb_lines(filePath)-1;
    double class_kdata[class_dsize], class_pkdata[class_dsize]; //data from Class for extrapolation
    
    read_classdata(class_kdata,class_pkdata,class_dsize,filePath);
    init_PL(class_kdata,class_pkdata,class_dsize);
    getkmax(kmin,&kmax,nu);
    printf("kmax=%lf\n",kmax);
    pars->kmax=kmax;
}

void free_vars(struct pars_struct *pars){
    int Ns=pars->Ns;
    int N=pars->N;

    free((void*)pars->filePath);
    gsl_vector_complex_free(pars->v1);
    gsl_vector_complex_free(pars->v2);
    gsl_vector_complex_free(pars->vkmu);

    free_vector_complex_arrays(pars->M_vv, Ns*N*N);
    free_vector_complex_arrays(pars->M_vvunsym, Ns*N);
    
    free(pars);
    free_PL();
}

void read_classdata(double *class_kdata,double *class_pkdata,int size,const char *filePath){

    FILE *filePointer;
    double value1, value2;
    filePointer = fopen(filePath, "r");
    if (filePointer == NULL) {
        perror("Error opening file");
    }
    int index = 0;
    while (fscanf(filePointer, "%lf %lf", &value1, &value2) != EOF) {
        class_kdata[index] = log(value1);
        class_pkdata[index] = log(value2);
        index++;
    }
    fclose(filePointer);
}

int numb_lines(const char *filePath){
    FILE *filePointer;
    int count = 0;
    char ch;

    filePointer = fopen(filePath, "r");

    if (filePointer == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    while ((ch = fgetc(filePointer)) != EOF) {
        if (ch == '\n') {
            count++;
        }
    }
    // Adding one because the last line may not end with a newline character
    if (ch != '\n' && count != 0) {
        count++;
    }

    fclose(filePointer);
    return count;
}



double linearExtrapolate(double x, double x1, double y1, double x2, double y2) {
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}
void init_PL(double *class_kdata,double *class_pkdata,int size){

    if (!initialized) {
    acc= gsl_interp_accel_alloc ();
    spline= gsl_spline_alloc(gsl_interp_cspline, size); 
    gsl_spline_init(spline,class_kdata,class_pkdata, size);
    initialized = 1;
    }
}
double PL(double k) {
    if (!initialized) {
        fprintf(stderr, "Error: Spline not initialized.\n");
        exit(EXIT_FAILURE);
    }
    return exp(gsl_spline_eval(spline, log(k), acc));
}
void free_PL() {
    if (acc) {
        gsl_interp_accel_free(acc);
        acc = NULL;
    }
    if (spline) {
        gsl_spline_free(spline);
        spline = NULL;
    }
}


void getkmax(double kmin,double *kmax,double nu){
    double kbest,errbest,errtest;
    double k,kbounda,kboundb;

    int N=pow(10,3);
    kbounda=PL(kmin)/pow(kmin,nu);    
    errbest=fabs(PL(*kmax)/pow(*kmax,nu)-kbounda);

    double delta=log(*kmax/kmin)/N; 

    for(int i=N/3; i<N; i++) {
        k=kmin*exp(i*delta);
        kboundb=PL(k)/pow(k,nu);
        errtest=fabs(kboundb-kbounda);
        if(errtest<errbest){
            errbest=errtest;
            kbest=k;
        }
    }
    *kmax=kbest;
}