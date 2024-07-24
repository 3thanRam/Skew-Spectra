#include "../include/randomrange.h"

static int config_handler(void* user, const char* section, const char* name, const char* value) {
    struct pars_struct* config = (struct pars_struct*)user;
    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

    if (MATCH("General", "N")) {
        config->N = atoi(value);
        double** ranges = malloc(config->N * sizeof(double*));
        for (int i = 0; i < config->N; i++) {
            ranges[i] = malloc(4 * sizeof(double));
        }
        config->ranges=ranges;
    }
    if(config->N!=0){
        for (int i = 0; i < config->N; i++) {
            char rangenumb[20];
            snprintf(rangenumb, sizeof(rangenumb), "ranges%dIR", i);
            if (MATCH("General", rangenumb)) {
                char valcopy[40];                
                strcpy(valcopy,value);
                char *ptr = strtok(valcopy, ",");
                config->ranges[i][0] = atof(ptr);
                ptr = strtok(NULL, ",");
                config->ranges[i][1] = atof(ptr);
            }
            snprintf(rangenumb, sizeof(rangenumb), "ranges%dUV", i);
            if (MATCH("General", rangenumb)) {
                char valcopy[40];                
                strcpy(valcopy,value);
                char *ptr = strtok(valcopy, ",");
                config->ranges[i][2] = atof(ptr);
                ptr = strtok(NULL, ",");
                config->ranges[i][3] = atof(ptr);
            }
        }
    }
    return 1;
}


void appenddouble(double **array,int *size,double val){
    if((val!=INFINITY)&&(val!=-INFINITY)){
        double *temp = realloc(*array, (*size + 1) * sizeof(double));
        if (temp == NULL) {
            printf("realloc failed\n");
            return; 
        }
        *array = temp;  
        (*array)[*size] = val; 
        *size = *size + 1; 
    }
}

bool containselem(double *arr, int size, int element) {
    for (int i = 0; i < size; i++) {
        if (arr[i] == element) {
            return true; 
    }}
    return false;
}

void reduceincrements(double **arr, int size){
    double prec=1e-4;
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if((i!=j)&&((fmod((*arr)[i],(*arr)[j])<(*arr)[i]*prec)||fabs(fmod((*arr)[i],(*arr)[j])-(*arr)[j])<(*arr)[i]*prec)){
                (*arr)[i]=0.;
    }}}
}

int numbnonzero(double *arr, int size){
    int zerocount=0;
    for (int i=0; i<size; i++) {
        if(arr[i]!=0.){
            zerocount++;
    }}
    return zerocount;
}

double minincrement(double **ranges,int N){
    double minincr;
    int size = 0;

    double *vals=NULL;
    for (int i=0; i<N; i++) {
        for (int j=0; j<4; j++) {
        appenddouble(&vals,&size,ranges[i][j]);
        //appenddouble(&vals,&size,ranges[i][1]);
    }}

    minincr=INFINITY;
    double d;
    double *increments=NULL;
    int isize=0;
    appenddouble(&increments,&isize,1.);
    
    for (int i=0; i<size; i++){for (int j=0; j<size; j++){
            if(i!=j){
                d=fabs(vals[i]-vals[j]);
                //printf("d:(%d,%d)=%f\n",i,j,d);
                d=d-floor(d);
                if((d>0.)&&(containselem(increments,isize,d)==0)){
                    appenddouble(&increments,&isize,d);
    }}}}
    
    int Cmax=isize;
    int count=0;
    while(numbnonzero(increments, isize)!=1 && count<Cmax){
    reduceincrements(&increments,isize);
    count++;}

    if(numbnonzero(increments, isize)!=1){
        printf("failed to find minimum increment\n");
    }

    for (int i=0; i<isize; i++) {
        //printf("%d:%f \n",i,increments[i]);
        if(increments[i]!=0.){
            minincr=increments[i];
    }}

    free(vals);
    free(increments);
    return minincr;
}

void shuffle(double **array, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1); 
        double temp = (*array)[i];
        (*array)[i] = (*array)[j];
        (*array)[j] = temp;
    }
}


void genconfig(double **config,int N,int Nc,double min,double max,double dnu){
    int Nint=(max-min)/dnu;
    for (int i=0; i<Nc; i++) {
        int r=rand()%Nint;
        (*config)[i]=min+dnu*r;
    }
    for (int i=0; i<N-Nc; i++) {
        (*config)[i+Nc]= (*config)[i%Nc];
    }
    shuffle(config,N);
}



void findbest(double **ranges,int N,double **configbest){
    int Ntests=1e6;
    int Nbest=N;

    double min=INFINITY;
    double max=-INFINITY;
    for (int i=0; i<N; i++){
        for (int j=0; j<4; j++){
            double rangevalij =ranges[i][j];
            if((rangevalij!=-INFINITY)&&(rangevalij<min)){
                min=rangevalij;
            }else if ((rangevalij!=INFINITY)&&(rangevalij>max)){
                max=rangevalij;
            }
    }}
    //printf("min:%f, max:%f\n",min,max);

    double dnu=minincrement(ranges,N)/10;
    min-=1;
    max+=1;
    int Ndivbest=2*N;

    double *testconfig = malloc(N * sizeof(double));
    for (int Nc=1; Nc<N+1; Nc++) {
        
        for(int nt=0;nt<Ntests;nt++){
            genconfig(&testconfig,N,Nc,min,max,dnu);
            int Ndivtest=Numbdivs(ranges, N,testconfig,dnu);
            
            if((Ndivtest<Ndivbest)||((Ndivtest==Ndivbest)&&(Nc<Nbest))){
                Ndivbest=Ndivtest;
                Nbest=Nc;
                free(*configbest);
                *configbest = malloc(N * sizeof(double));
                memcpy(*configbest, testconfig, N * sizeof(double));
    }}}
    if (*configbest != testconfig) {
        free(testconfig);  
    }

    mathematicaOutput(ranges, N,*configbest, dnu);
}


int Numbdivs(double **ranges,int N,double *nulist,double dnu){
    int Ndiv=0;
    for (int i=0; i<N; i++) {
        if((ranges[i][0]-dnu<nulist[i])&&(nulist[i]<ranges[i][1]+dnu)){
            Ndiv++;
        }
        if((ranges[i][2]-dnu<nulist[i])&&(nulist[i]<ranges[i][3]+dnu)){
            Ndiv++;
        }
        
        //if(nulist[i]<ranges[i][0]){
        //    Ndiv++;
        //}
        //if(nulist[i]>ranges[i][1]){
        //    Ndiv++;
        //}
    }
    return Ndiv;
}

void findnuvals(){
    srand (time ( NULL));
    struct pars_struct *pars = malloc(sizeof( *pars));
    pars->N=0;
    if (ini_parse("ranges.ini", config_handler, pars) < 0) {
        printf("Can't load ini file\n");
        exit(EXIT_FAILURE);
    }
    int N=pars->N;
    double **ranges=pars->ranges;
    //for (int i = 0; i < N; i++) {
    //    printf("i:%d ,IR:(%f,%f) UV:(%f,%f)\n",i,ranges[i][0],ranges[i][1],ranges[i][2],ranges[i][3]);
    //}

    double *nuconfig = malloc(N * sizeof(double));
    
    findbest(ranges, N,&nuconfig);
    
    
    free(nuconfig);
    
    for (int i = 0; i < N; i++) {
        free(ranges[i]);
    }
    free(ranges);
}


void printing(double **ranges,int N,double *nuconfig,double dnu){
    int dupl=0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<i; j++) {
            if((i!=j)&&fabs(nuconfig[i]-nuconfig[j])<1e-5){
                dupl++;
            }
        }}
    printf("Number of different values: %d,  Number of divergences: %d\n", N-dupl,Numbdivs(ranges, N,nuconfig, dnu));
    for (int i=0; i<N; i++) {
        printf("nu[%d]= %f  Bound conditions:%d,%d\n", i,nuconfig[i],((ranges[i][0]<nuconfig[i])&&(nuconfig[i]<ranges[i][1])),((ranges[i][2]<nuconfig[i])&&(nuconfig[i]<ranges[i][3])));
    }
}

void mathematicaOutput(double **ranges,int N,double *nuconfig,double dnu){
    int dupl=0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<i; j++) {
            if((i!=j)&&fabs(nuconfig[i]-nuconfig[j])<1e-5){
                dupl++;
            }
        }}
    printf("Nc=%d,  Ndivs=%d\n", N-dupl,Numbdivs(ranges, N,nuconfig,dnu));
    for (int i=0; i<N; i++) {
        printf("nu[[%d]]= %f \n", i,nuconfig[i]);
    }
    printf("insdDiv=[");
    for (int i=0; i<N; i++) {
        if(i!=0){printf(",");}
        printf("[%d,%d]",((ranges[i][0]<nuconfig[i])&&(nuconfig[i]<ranges[i][1])),((ranges[i][2]<nuconfig[i])&&(nuconfig[i]<ranges[i][3])));
    }
    printf("]\n");
}

