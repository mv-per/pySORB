#ifndef TEMP_DEPENDENCE_H
#define TEMP_DEPENDENCE_H


/*
    LANGMUIR
*/

// Equation 3.2-18b (DO 1998)
double langmuir_calc_b(double b_inf, double Q, double T){
    return b_inf * exp(Q/R_G/T_0 * (T_0/T - 1.));
}

// Equation 3.2-18d (DO 1998)
double langmuir_n_max_T(double n_max_0, double QSI, double T){
    return n_max_0 * exp(QSI * (1. - T/T_0));
}

/*
    SIPS
*/

// Equation 3.2-18b (DO 1998)
double sips_calc_b(double b_inf, double Q, double T){
    return b_inf * exp(Q/R_G/T_0 * (T_0/T - 1.));
}

// Equation 3.2-18c (DO 1998)
double sips_inverse_n(double n_0, double alpha, double T){
    return 1./n_0 + alpha * (1. - T_0/T);
}

// Equation 3.2-18d (DO 1998)
double sips_n_max_T(double n_max_0, double QSI, double T){
    return n_max_0 * exp(QSI * (1. - T/T_0));
}

/* TOTH
*/

double toth_calc_b(double b_inf, double Q, double T){
    return b_inf * exp(Q/R_G/T_0 * (T_0/T - 1.));
}

double toth_calc_t(double t_0, double alpha, double T){
    return t_0  + alpha * (1. - T_0/T);
}

double toth_n_max_T(double n_max_0, double QSI, double T){
    return n_max_0 * exp(QSI * (1. - T/T_0));
}



/*
    UNILAN
*/

// Equation 3.2-23c (DO 1998)
double unilan_calc_s(double delta_E, double T){
    return delta_E/2./R_G/T;
}

// Equation 3.2-23b (DO 1998)
double unilan_calc_b(double b_inf, double E_MED, double T){
    return b_inf * exp(-E_MED/R_G/T_0 * (1. - T_0/T));
}



#endif