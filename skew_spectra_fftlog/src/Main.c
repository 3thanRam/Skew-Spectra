#include "../include/Util.h"
#include "../include/Plotting.h"
#include "../include/fftw_funcs.h"



int main(int argc, char **argv)
{
    
    struct pars_struct *params = malloc(sizeof( *params));
    init_cosmoparams(params);

    int Ns=params->Ns;
    int N=params->N;

    double krange[N]; //Array of N wavenumbers from kmin to kmax with log spacing
    
    double **PSn=malloc(Ns * sizeof(double *));
    for (size_t i = 0; i < Ns; i++) {
        PSn[i] = malloc(N * sizeof(double *));
    }

    struct timespec start, end;
    double elapsed;
    timespec_get(&start, TIME_UTC);
    
    init_class(params);

    get_krange(krange,params);
    get_eta(params);
    
    get_fouriercoefs(krange,params);
    
    get_matrices(params);
    fill_Sn(krange,N,params,PSn);
    //muk_integrate(krange, N,PSn,params);

    timespec_get(&end, TIME_UTC);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Elapsed time: %f seconds\n", elapsed);

    plot_data(params,N,krange,PSn);

    free_vars(params);
    return 0;
}
