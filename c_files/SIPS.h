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


/*
    param = [n_max_0, QSI, b_inf, Q, n_0, alpha]

*/
double calculate_sips_T(double *P, double *T, double* param )
{
    double b = sips_calc_b(param[2], param[3], *T);
    double n_max = sips_n_max_T(param[0], param[1], *T);
    double inverse_n = sips_inverse_n(param[4], param[5], *T);

    // printf("%f \n", (n_max*pow(b * *P, inverse_n))/(1. + pow(b * *P, inverse_n)));
    return  (n_max*pow(b * *P, inverse_n))/(1. + pow(b * *P, inverse_n));
}



double minimize_sips(double P, double n_exp, double* param){
    double n_calc = calculate_sips(&P, param);
    // printf(%f, n_calc);
    // return (fabs(n_exp-n_calc)/n_exp)*1000.;
    return pow(fabs(n_exp-n_calc),2);
    // return (n_exp-n_calc);
}

#endif