#include "../include/Plotting.h"


void sendPdata(FILE *gnuplotPipe, double *kdata, double *P, int N) {
    for (int i = 0; i < N; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", kdata[i], P[i]);
    }
    fprintf(gnuplotPipe, "e\n");
}

void plot_data(struct pars_struct *pars,int N, double *karray, double **PSn) {

    double nu=pars->nu;
    double z=pars->z;
    //int size=pars->size;
    int Ns=pars->Ns;


    const char *colors[] = {"red", "blue", "green", "yellow", "black", "white", "purple"};
    int numColors = sizeof(colors) / sizeof(colors[0]);

    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        fprintf(stderr, "Error opening pipe to Gnuplot.\n");
        return;
    }

    printf("Start plotting\n");

    // Setup PDF output and Multiplot environment
    fprintf(gnuplotPipe, "set terminal pdf enhanced \n");
    fprintf(gnuplotPipe, "set output 'skewspectra_plot.pdf'\n");
    
    fprintf(gnuplotPipe, "set multiplot layout %d,1 title 'Skew-spectra P_{Sn,{/Symbol d}} [(Mpc/h)^3] as a function of k [h/Mpc] for redshift z=%.2f \\& bias {/Symbol n}=%.1f'\n", Ns, z, nu);
    
    for (int i = 0; i < Ns; i++) {
        //fprintf(gnuplotPipe, "set title 'Subplot %d: PS_{%i}'\n", i + 1, i + 1);
        fprintf(gnuplotPipe, "set xlabel 'k'\n");
        fprintf(gnuplotPipe, "set ylabel 'P_{S%i,{/Symbol d}}'\n", i + 1);
        fprintf(gnuplotPipe, "plot '-' using 1:2 title 'PS_{%i}' with lines linecolor '%s'\n", i + 1, colors[i % numColors]);
        sendPdata(gnuplotPipe, karray, PSn[i], N);
    }

    fprintf(gnuplotPipe, "unset multiplot\n");
    pclose(gnuplotPipe);
    printf("Plotting complete.\n");
}