#ifndef REDLICH_PETERSON_H
#define REDLICH_PETERSON_H


/*

     param[0] = K_RP
     param[1] = a_RP
     param[2] = g_RP
*/
double calculate_redlichpeterson(double *P, double* param)
{
    return  (param[0] * *P)/ (1. + param[1] *pow(*P, param[2]));
}


double minimize_redlichpeterson(double P, double n_exp, double param[]){
    double n_calc = calculate_redlichpeterson(&P, param);
    // printf(%f, n_calc);
    return (fabs(n_exp-n_calc)/n_exp)*1000.;
}

#endif