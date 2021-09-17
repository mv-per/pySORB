

// Langmuir


#ifndef LANGMUIR_H
#define LANGMUIR_H

/*
	Returns single n_calc using Langmuir
*/
double CalculateLangmuir(double *P, double &n_i_max, double &b)
{
    return n_i_max * b * *P / (1.0 + b * *P);
}

/*
	Returns partial n_calc and total n_calc using Extended Langmuir
*/
// double * CalculateExtendedLangmuir(vector<vector<double>>& param, vector<double>& y, double *P, int *ncomp)
// {
//     int i;
//     double *res = (double *) malloc( *ncomp+1 * sizeof( double ) );

//     double nT = 0.;
//     double sum1 = 0.;
//     double nmix = 0.;
//     for (i = 0; i < *ncomp; i++){
//         sum1 += param[i][1] * *P * y[i];
//         nmix += param[i][0]*y[i];
//     }

//     for (i = 0; i < *ncomp; i++){
//         res[i] = nmix * param[i][1] * *P * y[i] / (1. + sum1);
//         nT += res[i];
//     }

//     for (i = 0; i < *ncomp; i++){
//         res[i] /= nT;
//     }
    
//     res[*ncomp] = nT;

//     return res;
// }

double MinimizeLANGMUIR(double P, double n_exp, double param[]){
    double n_calc = CalculateLangmuir(&P, param[0],param[1]);
    // printf(%f, n_calc);
    return fabs(n_exp-n_calc);
}



#endif