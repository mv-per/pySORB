#ifndef SIPS_H
#define SIPS_H

/*
	Returns single n_calc using Sips
    param[0] = n_i_max
    param[1] = b
    param[2] = n
*/
double calculate_sips(double *P, double* param )
{
    
    return  (param[0]*pow(param[1] * *P, 1./param[2]))/(1. + pow(param[1] * *P, 1. / param[2]));
}


double minimize_sips(double P, double n_exp, double* param){
    double n_calc = calculate_sips(&P, param);
    // printf(%f, n_calc);
    // return (fabs(n_exp-n_calc)/n_exp)*1000.;
    return pow(fabs(n_exp-n_calc),2);
    // return (n_exp-n_calc);
}

#endif