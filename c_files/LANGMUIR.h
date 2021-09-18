

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


double calculate_langmuir_T(double *P, double *T, double* param)// double &n_i_max, double &b)
{
    // printf("param = [%f, %f, %f, %f]", param[0], param[1], param[2], param[3]);
    double n_max = langmuir_n_max_T(param[0], param[1], *T);
    double b = langmuir_calc_b(param[2], param[3], *T);
    double np[] = {n_max,b};
    return calculate_langmuir(P, np);
}

double minimize_langmuir(double P, double n_exp, double* param){
    double n_calc = calculate_langmuir(&P, param);
    // printf(%f, n_calc);
    return (fabs(n_exp-n_calc)/n_exp)*1000.;
}



#endif