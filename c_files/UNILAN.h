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

#endif