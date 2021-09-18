#ifndef TOTH_H
#define TOTH_H


/*

     param[0] = n_max
     param[1] = b
     param[2] = t
*/
double calculate_toth(double *P, double* param)
{
    return  (param[0] * param[1] * *P)/ pow(1. + pow(param[1] * *P, param[2]),1./param[2]);
}


/*
    param = [n_max_0, QSI, b_inf, Q, n_0, alpha]
*/
double calculate_toth_T(double *P, double *T, double* param )
{
    double b = toth_calc_b(param[2], param[3], *T);
    double n_max = toth_n_max_T(param[0], param[1], *T);
    double t = toth_calc_t(param[4], param[5], *T);
    double newparam[] = {n_max, b, t};
    return  calculate_toth(P, newparam);
}

double minimize_toth(double P, double n_exp, double param[]){
    double n_calc = calculate_toth(&P, param);
    // printf(%f, n_calc);
    // return (fabs(n_exp-n_calc)/n_exp)*1000.;
    return fabs(n_exp-n_calc);
}

#endif