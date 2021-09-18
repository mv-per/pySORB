#ifndef FREUNDLICH_H
#define FREUNDLICH_H


/*
    param[0] = K
    param[1] = n (heterogeneity)
*/
double calculate_freundlich(double *P, double* param)
{
    return  (param[0] * pow(*P, 1.0/ param[1]));
}


/*
    param[0] = K_0
    param[1] = alpha
    param[2] = A_0
*/
double calculate_freundlich_T(double *P, double *T, double* param)
{
    double inverse_n = R_G * *T / param[2];

    double K = param[0] * exp(- param[1] * inverse_n);

    return K * pow(*P, inverse_n);
}







#endif