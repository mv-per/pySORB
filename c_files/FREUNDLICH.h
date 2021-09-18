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


double minimize_freundlich(double P, double n_exp, double param[]){
    double n_calc = calculate_freundlich(&P, param);
    // printf(%f, n_calc);
    // return fabs(n_exp-n_calc);
    return pow(fabs(n_exp-n_calc),2);
}




#endif