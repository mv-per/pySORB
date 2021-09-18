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


double minimize_toth(double P, double n_exp, double param[]){
    double n_calc = calculate_toth(&P, param);
    // printf(%f, n_calc);
    // return (fabs(n_exp-n_calc)/n_exp)*1000.;
    return fabs(n_exp-n_calc);
}

#endif