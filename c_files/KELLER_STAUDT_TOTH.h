#ifndef KELLER_STAUDT_TOTH_H
#define KELLER_STAUDT_TOTH_H


/*
    DUONG DO BOOK
    param[0] = N_max
    param[1] = b
    param[2] = alpha_m
    param[3] = beta
*/
double calculate_kst(double *P, double* param)
{
    double alpha = (1. + param[2]*param[3])/(1.0 + param[3]* *P);
    return  (param[0] * param[2] * param[1] * *P /pow(1. + pow(param[1] * *P, alpha), 1./alpha));
}

#endif