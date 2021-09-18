#ifndef JENSEN_SEATON_H
#define JENSEN_SEATON_H

/*
   K is the Henry constant
   b is the compressibility of the adsorbed phase and
   c an empirical constant.

*/
double calculate_jensenseaton(double *P, double* param){
    // printf("param = [%f %f %f %f] = %f \n", param[0], param[1], param[2], param[3], param[0] * *P * pow(1. + pow(param[0] * *P / (param[1] * (1. + param[2]* *P)),param[3]),-1./param[3]));
    return param[0] * *P * pow(1. + pow(param[0] * *P / (param[1] * (1. + param[2]* *P)),param[3]),-1./param[3]);
}

#endif


