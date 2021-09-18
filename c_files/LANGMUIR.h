

// Langmuir


#ifndef LANGMUIR_H
#define LANGMUIR_H

/*
	Returns single n_calc using Langmuir
    param[0] = n_i_max
    param[1] = b
*/
double calculate_langmuir(double *P, double* param)// double &n_i_max, double &b)
{
    return param[0] * param[1] * *P / (1.0 + param[1] * *P);
}

double minimize_langmuir(double P, double n_exp, double* param){
    double n_calc = calculate_langmuir(&P, param);
    // printf(%f, n_calc);
    return (fabs(n_exp-n_calc)/n_exp)*1000.;
}



#endif