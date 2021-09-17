#include <malloc.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "nmsimplex.h"


#include "LANGMUIR.h"
#include "GNUPLOT.h"

// #define LEN(arr) 


void optconst(double x[], int n){
    for (size_t i = 0; i < n; i++){
        x[i] = fabs(x[i]);
    }
}




int main()
{   
    double Pfac = 1e6;
    size_t i;
    double ini[] = {1, 1e-5};
    // double param_CO2[] = {7.24921072e+00, 2.49021651e-03, 5.91326737e+00, 1.56288967e+00};
    // double param_N2[] = {4.54307458e+00, 2.70652617e-06, 1.39409355e+00, 9.05687013e-01};
    // double* param[] = {param_CO2, param_N2};  



    vector<double> nexp = {2.11, 2.64, 3.144, 3.746, 4.088, 4.547, 4.883, 5.125, 5.283, 5.445, 5.6, 5.724, 5.862, 5.994};
    vector<double> pexp = {0.003, 0.0058, 0.0105, 0.0218, 0.0332, 0.0618, 0.1021, 0.1518, 0.2002, 0.2823, 0.38, 0.5, 0.7, 1};

    size_t n = pexp.size();

    for( i=0; i<n; i++){
        pexp[i] = pexp[i]*Pfac;
    }
    
    

    static auto OptiminFunc = [&](double x[]) {
        double param[] = {x[0], x[1]};
        double fobj = 0.;
        for (size_t i = 0; i< n; i++){
            fobj += MinimizeLANGMUIR(pexp[i], nexp[i], param);
        }
        // printf("%f\n", fobj);
        return fobj;
    };
    double (*MINFUN)(double *) = [](double x[]) { return OptiminFunc(x); };

    double * new_param = simplex(MINFUN, ini, 4, 1e-6, .5, optconst);

    // double new_param[] = {5.72399534e+00, 1.64655106e-04};

    double * ncalc = (double *) malloc (n * sizeof (double));
    double * pcalc =  (double *) malloc (n * sizeof (double));
    
    #pragma omp parallel for
    printf("P \t nexp \t ncalc\n");
    for( i=0; i<n; i++){
        pcalc[i] = pexp[i];
        ncalc[i] = CalculateLangmuir(&pexp[i],new_param[0],new_param[1]);
        printf("%f \t %f \t %f \n", pexp[i], nexp[i],  ncalc[i]);
    }
    

    gplot plt;
    plt.plot(pcalc,ncalc , &n);
    plt.show();


    return 0;
}
