#ifndef UNILAN_H
#define UNILAN_H



/*
    param[0] = n_i_max
    param[1] = b_med
    param[2] = s

*/

double calculate_unilan(double *P, double* param)
{
    double alpha = (1. + param[2]*param[3])/(1.0 + param[3]* *P);
    return  param[0]/2./param[2] * log((1. + param[1] * exp(param[2]) * *P)/ (1. + param[1] * exp(-param[2]) * *P));
}



/*
    param[0] = n_i_max
    param[1] = b_0
    param[2] = E_MED
    param[3] = deltaE

*/
double calculate_unilan_T(double *P, double *T, double* param)
{
    double b_med = unilan_calc_b(param[1], param[2], *T);
    double s = unilan_calc_s(param[3], *T);
    return  param[0]/2./s * log((1. + b_med * exp(s) * *P)/ (1. + b_med * exp(-s) * *P));
}


#endif