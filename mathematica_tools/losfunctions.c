#include "losfunctions.h"
double complex A1(double complex nu1, double complex nu2)
{
return J(-1 + nu1,nu2)/2. - J(nu1,-1 + nu2)/2. + J(nu1,nu2)/2.;
}

double complex A2(double complex nu1, double complex nu2)
{
return -0.125*J(-2 + nu1,nu2) + J(-1 + nu1,-1 + nu2)/4. + J(-1 + nu1,nu2)/4. - J(nu1,-2 + nu2)/8. + J(nu1,-1 + nu2)/4. - J(nu1,nu2)/8.;
}

double complex B2(double complex nu1, double complex nu2)
{
return (3*J(-2 + nu1,nu2))/8. - (3*J(-1 + nu1,-1 + nu2))/4. + J(-1 + nu1,nu2)/4. + (3*J(nu1,-2 + nu2))/8. - (3*J(nu1,-1 + nu2))/4. + (3*J(nu1,nu2))/8.;
}

double complex A3(double complex nu1, double complex nu2)
{
return (-3*J(-3 + nu1,nu2))/16. + (9*J(-2 + nu1,-1 + nu2))/16. + (3*J(-2 + nu1,nu2))/16. - (9*J(-1 + nu1,-2 + nu2))/16. + (3*J(-1 + nu1,-1 + nu2))/8. + (3*J(-1 + nu1,nu2))/16. + (3*J(nu1,-3 + nu2))/16. - (9*J(nu1,-2 + nu2))/16. + (9*J(nu1,-1 + nu2))/16. - (3*J(nu1,nu2))/16.;
}

double complex B3(double complex nu1, double complex nu2)
{
return (5*J(-3 + nu1,nu2))/16. - (15*J(-2 + nu1,-1 + nu2))/16. + (3*J(-2 + nu1,nu2))/16. + (15*J(-1 + nu1,-2 + nu2))/16. - (9*J(-1 + nu1,-1 + nu2))/8. + (3*J(-1 + nu1,nu2))/16. - (5*J(nu1,-3 + nu2))/16. + (15*J(nu1,-2 + nu2))/16. - (15*J(nu1,-1 + nu2))/16. + (5*J(nu1,nu2))/16.;
}

double complex A4(double complex nu1, double complex nu2)
{
return (3*J(-4 + nu1,nu2))/128. - (3*J(-3 + nu1,-1 + nu2))/32. - (3*J(-3 + nu1,nu2))/32. + (9*J(-2 + nu1,-2 + nu2))/64. + (3*J(-2 + nu1,-1 + nu2))/32. + (9*J(-2 + nu1,nu2))/64. - (3*J(-1 + nu1,-3 + nu2))/32. + (3*J(-1 + nu1,-2 + nu2))/32. + (3*J(-1 + nu1,-1 + nu2))/32. - (3*J(-1 + nu1,nu2))/32. + (3*J(nu1,-4 + nu2))/128. - (3*J(nu1,-3 + nu2))/32. + (9*J(nu1,-2 + nu2))/64. - (3*J(nu1,-1 + nu2))/32. + (3*J(nu1,nu2))/128.;
}

double complex B4(double complex nu1, double complex nu2)
{
return (-15*J(-4 + nu1,nu2))/64. + (15*J(-3 + nu1,-1 + nu2))/16. + (3*J(-3 + nu1,nu2))/16. - (45*J(-2 + nu1,-2 + nu2))/32. + (9*J(-2 + nu1,-1 + nu2))/16. + (3*J(-2 + nu1,nu2))/32. + (15*J(-1 + nu1,-3 + nu2))/16. - (27*J(-1 + nu1,-2 + nu2))/16. + (9*J(-1 + nu1,-1 + nu2))/16. + (3*J(-1 + nu1,nu2))/16. - (15*J(nu1,-4 + nu2))/64. + (15*J(nu1,-3 + nu2))/16. - (45*J(nu1,-2 + nu2))/32. + (15*J(nu1,-1 + nu2))/16. - (15*J(nu1,nu2))/64.;
}

double complex C4(double complex nu1, double complex nu2)
{
return (35*J(-4 + nu1,nu2))/128. - (35*J(-3 + nu1,-1 + nu2))/32. + (5*J(-3 + nu1,nu2))/32. + (105*J(-2 + nu1,-2 + nu2))/64. - (45*J(-2 + nu1,-1 + nu2))/32. + (9*J(-2 + nu1,nu2))/64. - (35*J(-1 + nu1,-3 + nu2))/32. + (75*J(-1 + nu1,-2 + nu2))/32. - (45*J(-1 + nu1,-1 + nu2))/32. + (5*J(-1 + nu1,nu2))/32. + (35*J(nu1,-4 + nu2))/128. - (35*J(nu1,-3 + nu2))/32. + (105*J(nu1,-2 + nu2))/64. - (35*J(nu1,-1 + nu2))/32. + (35*J(nu1,nu2))/128.;
}
