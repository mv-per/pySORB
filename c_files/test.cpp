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


#include "REDLICH_PETERSON.h"
#include "LANGMUIR.h"
#include "FREUNDLICH.h"
#include "SIPS.h"
#include "TOTH.h"
#include "UNILAN.h"
#include "KELLER_STAUDT_TOTH.h"
#include "ERRORS.h"

// #define LEN(arr) 


void optconst(double x[], int n){
    for (size_t i = 0; i < n; i++){
        x[i] = fabs(x[i]);
    }
}


// TEST ALL

int main()
{   
    double decrease = .95;
    double further_decrease = .75;
    double tesst;

    double Pfac = 1e6;
    size_t i,p,n, RUN;
    double ini_langmuir[2] = {1,1};
    double ini_freundlich[2] = {1,1};
    double ini_sips[3] = {1,1,1};
    double ini_redlichpeterson[3] = {1,1,1};
    double ini_toth[3] = {1,1,1};
    double ini_unilan[3] = {1,1,1};
    double ini_kst[4] = {1,1,1,1};


    RUN = 0;
    size_t NUM_OF_OPTMIZATION = 50;
    double TOL = 1e-12;
    double SCALE = .001;
    string ERR_FUNC = "CHI_2";

    /*
        SSE HYBRID  ARE EABS MPSD SDRE R_S CHI_2

    */
 
    vector<double> nexp_entrada = {2.11, 2.64, 3.144, 3.746, 4.088, 4.547, 4.883, 5.125, 5.283, 5.445, 5.6, 5.724, 5.862, 5.994};
    vector<double> pexp = {0.003, 0.0058, 0.0105, 0.0218, 0.0332, 0.0618, 0.1021, 0.1518, 0.2002, 0.2823, 0.38, 0.5, 0.7, 1};



    n = pexp.size();
    
    double* ncalc_opt = (double *) malloc (n * sizeof (double));
    double* nexp = (double *) malloc (n * sizeof (double));
    // double* pexp = (double *) malloc (n * sizeof (double));

    for( i=0; i<n; i++){
        pexp[i] = pexp[i]*Pfac;
        nexp[i] = nexp_entrada[i];
    }

// static auto error_functions = [&](double* n_exp, double* n_calc, size_t n, size_t p=1,string ERROR_FUNC)
    // double (*ERR_FUNC)(double *) = [](double x[]) { return error_functions(n_exp, n_calc, n, p,ERR_FUNC); };


    start_optimization: 
        static auto OptiminFunc = [&](double x[], string ISOTHERM) {

            if (ISOTHERM == "Freundlich")
            {
                double param[] = {x[0], x[1]};
                p = 2;
                for (i = 0; i< n; i++)
                    {
                    ncalc_opt[i] = calculate_freundlich(&pexp[i], param);
                    }
            }
            else if (ISOTHERM == "Langmuir")
            {
                double param[] = {x[0], x[1]};
                p = 2;
                for (i = 0; i< n; i++){
                    ncalc_opt[i] = calculate_langmuir(&pexp[i], param);
                }
            }
            else if (ISOTHERM == "Sips")
            {
                double param[] = {x[0], x[1], x[2]};
                p = 3;
                for (i = 0; i< n; i++){
                    ncalc_opt[i] = calculate_sips(&pexp[i],  param);
                }
            }

            else if (ISOTHERM == "Unilan")
            {
                double param[] = {x[0], x[1], x[2]};
                p = 3;
                for (i = 0; i< n; i++){
                    ncalc_opt[i] = calculate_unilan(&pexp[i],  param);
                }
            }
            else if (ISOTHERM == "Redlich-Peterson")
            {
                double param[] = {x[0], x[1], x[2]};
                p = 3;
                for (i = 0; i< n; i++)
                {
                    ncalc_opt[i] = calculate_redlichpeterson(&pexp[i], param);
                }
            }

            else if (ISOTHERM == "Toth")
            {
                double param[] = {x[0], x[1], x[2]};
                p = 3;
                for (i = 0; i< n; i++)
                {
                    ncalc_opt[i] = calculate_toth(&pexp[i], param);
                }
            }

            // Keller et al. (BOOK DO)
            else if (ISOTHERM == "kst")
            {
                double param[] = {x[0], x[1], x[2], x[3]};
                p = 3;
                for (i = 0; i< n; i++)
                {
                    ncalc_opt[i] = calculate_kst(&pexp[i], param);
                }
            }

            return get_error(nexp, ncalc_opt, n, p, ERR_FUNC);
        };


        double (*MINFUN)(double *) = [](double x[]) { return OptiminFunc(x, "Freundlich"); };
        double * new_param_freundlich = simplex(MINFUN, ini_freundlich, 2, TOL, SCALE, optconst);

        double (*MINFUN2)(double *) = [](double x[]) { return OptiminFunc(x, "Langmuir"); };
        double * new_param_langmuir = simplex(MINFUN2, ini_langmuir, 2, TOL, SCALE, optconst);

        if (RUN == 0){
            ini_sips[0] = new_param_langmuir[0];
            ini_sips[1] = new_param_langmuir[1];
        }
        double (*MINFUN3)(double *) = [](double x[]) { return OptiminFunc(x, "Sips"); };
        double * new_param_sips = simplex(MINFUN3, ini_sips, 3, TOL, SCALE, optconst);     

        double (*MINFUN4)(double *) = [](double x[]) { return OptiminFunc(x, "Redlich-Peterson"); };
        double * new_param_redlichpeterson = simplex(MINFUN4, ini_redlichpeterson, 3, TOL, SCALE, optconst);
    
        double (*MINFUN5)(double *) = [](double x[]) { return OptiminFunc(x, "Toth"); };
        double * new_param_toth = simplex(MINFUN5, ini_toth, 3, TOL, SCALE, optconst); 

        double (*MINFUN6)(double *) = [](double x[]) { return OptiminFunc(x, "kst"); };
        double * new_param_kst = simplex(MINFUN6, ini_kst, 4, TOL, SCALE, optconst);      

        double (*MINFUN7)(double *) = [](double x[]) { return OptiminFunc(x, "Unilan"); };
        double * new_param_unilan = simplex(MINFUN7, ini_unilan, 3, TOL, SCALE, optconst);
    
        RUN++;

        if (RUN < NUM_OF_OPTMIZATION){
            // if (RUN > 3) {ERR_FUNC = "MPSD";}
            if (RUN % 2){
                tesst = decrease;
            }
            else{ tesst = further_decrease;}
            //2 parameters isotherms
            for (i = 0; i < 2; i++){
                ini_langmuir[i] = new_param_langmuir[i] *tesst;
                ini_freundlich[i] = new_param_freundlich[i]*tesst;
            }
            //3 parameters isotherms
            for (i = 0; i < 3; i++){
                ini_sips[i] = new_param_sips[i]*tesst;
                ini_redlichpeterson[i] = new_param_redlichpeterson[i]*tesst;
                ini_toth[i] = new_param_toth[i]*tesst;
                ini_unilan[i] = new_param_unilan[i]*tesst;
            }

            //4 parameters isotherms
            for (i = 0; i < 4; i++){
                ini_kst[i] = new_param_kst[i]*tesst;
            }

            // printf("%d", RUN);
            goto start_optimization;
        }
        printf("\n");

        // printf("Isotherm model: Freundlich\n K_F = %f\n n = %f \n", new_param_freundlich[0], new_param_freundlich[1]);
        // printf("Isotherm model: Langmuir\n n_max = %f\n b = %e \n", new_param_langmuir[0], new_param_langmuir[1]);  
        // printf("Isotherm model: Sips\n K_S = %f\n a_S = %f\n beta_S = %f \n", new_param_sips[0], new_param_sips[1], new_param_sips[2]);
        // printf("Isotherm model: Redlich-Peterson \n K_R = %f \n a_R = %f \n g_R = %f \n", new_param_redlichpeterson[0], new_param_redlichpeterson[1], new_param_redlichpeterson[2]);
        // printf("Isotherm model: Toth \n K_T = %f \n a_T = %f \n t = %f \n", new_param_toth[0], new_param_toth[1], new_param_toth[2]);

        // double new_param[] = {5.72399534e+00, 1.64655106e-04};

        double * ncalc_langmuir = (double *) malloc (n * sizeof (double));
        double * ncalc_freundlich = (double *) malloc (n * sizeof (double));
        double * ncalc_redlichpeterson = (double *) malloc (n * sizeof (double));
        double * ncalc_sips = (double *) malloc (n * sizeof (double));
        double * ncalc_toth = (double *) malloc (n * sizeof (double));
        double * ncalc_kst = (double *) malloc (n * sizeof (double));
        double * ncalc_unilan = (double *) malloc (n * sizeof (double));

        double * pcalc =  (double *) malloc (n * sizeof (double));
        
        // #pragma omp parallel for
        printf("P \t nexp \t ncalc(freundlich) \t ncalc(Langmuir) \t ncalc(Sips) \t ncalc(Toth) \t ncalc(Unilan) \t ncalc(Redlich-Peterson) \t ncalc(kst)\n");
        for( i=0; i<n; i++){
            pcalc[i] = pexp[i]/Pfac;
            ncalc_langmuir[i] = calculate_langmuir(&pexp[i],new_param_langmuir);
            ncalc_freundlich[i] = calculate_freundlich(&pexp[i],new_param_freundlich);
            ncalc_sips[i] = calculate_sips(&pexp[i],new_param_sips);
            ncalc_toth[i] = calculate_toth(&pexp[i],new_param_toth);
            ncalc_redlichpeterson[i] = calculate_redlichpeterson(&pexp[i],new_param_redlichpeterson);
            ncalc_kst[i] = calculate_kst(&pexp[i],new_param_kst);
            ncalc_unilan[i] = calculate_unilan(&pexp[i],new_param_unilan);

            printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n", pcalc[i], nexp[i],  ncalc_freundlich[i], 
                ncalc_langmuir[i],  ncalc_sips[i], ncalc_toth[i],  ncalc_unilan[i],   ncalc_redlichpeterson[i], ncalc_kst[i]);
        }
    

    // gplot plt;
    // plt.plot(pcalc,ncalc , &n);
    // plt.show();


    return 0;
}







// // TEST FREUNDLICH
// int main()
// {   
//     double Pfac = 1e6;
//     size_t i;
//     double ini[] = {1, 1e-5};
//     // double param_CO2[] = {7.24921072e+00, 2.49021651e-03, 5.91326737e+00, 1.56288967e+00};
//     // double param_N2[] = {4.54307458e+00, 2.70652617e-06, 1.39409355e+00, 9.05687013e-01};
//     // double* param[] = {param_CO2, param_N2};  



//     vector<double> nexp = {2.11, 2.64, 3.144, 3.746, 4.088, 4.547, 4.883, 5.125, 5.283, 5.445, 5.6, 5.724, 5.862, 5.994};
//     vector<double> pexp = {0.003, 0.0058, 0.0105, 0.0218, 0.0332, 0.0618, 0.1021, 0.1518, 0.2002, 0.2823, 0.38, 0.5, 0.7, 1};

//     size_t n = pexp.size();

//     for( i=0; i<n; i++){
//         pexp[i] = pexp[i]*Pfac;
//     }
    
    

//     static auto OptiminFunc = [&](double x[]) {
//         double param[] = {x[0], x[1]};
//         double fobj = 0.;
//         for (size_t i = 0; i< n; i++){
//             fobj += minimize_freundlich(pexp[i], nexp[i], param);
//         }
//         // printf("%f\n", fobj);
//         return fobj;
//     };
//     double (*MINFUN)(double *) = [](double x[]) { return OptiminFunc(x); };

//     double * new_param = simplex(MINFUN, ini, 4, 1e-6, .5, optconst);

//     // double new_param[] = {5.72399534e+00, 1.64655106e-04};

//     double * ncalc = (double *) malloc (n * sizeof (double));
//     double * pcalc =  (double *) malloc (n * sizeof (double));
    
//     #pragma omp parallel for
//     printf("P \t nexp \t ncalc\n");
//     for( i=0; i<n; i++){
//         pcalc[i] = pexp[i];
//         ncalc[i] = calculate_freundlich(&pexp[i],new_param[0],new_param[1]);
//         printf("%f \t %f \t %f \n", pexp[i], nexp[i],  ncalc[i]);
//     }
    

//     gplot plt;
//     plt.plot(pcalc,ncalc , &n);
//     plt.show();


//     return 0;
// }

















// TEST LANGMUIR
// int main()
// {   
//     double Pfac = 1e6;
//     size_t i;
//     double ini[] = {1, 1e-5};
//     // double param_CO2[] = {7.24921072e+00, 2.49021651e-03, 5.91326737e+00, 1.56288967e+00};
//     // double param_N2[] = {4.54307458e+00, 2.70652617e-06, 1.39409355e+00, 9.05687013e-01};
//     // double* param[] = {param_CO2, param_N2};  



//     vector<double> nexp = {2.11, 2.64, 3.144, 3.746, 4.088, 4.547, 4.883, 5.125, 5.283, 5.445, 5.6, 5.724, 5.862, 5.994};
//     vector<double> pexp = {0.003, 0.0058, 0.0105, 0.0218, 0.0332, 0.0618, 0.1021, 0.1518, 0.2002, 0.2823, 0.38, 0.5, 0.7, 1};

//     size_t n = pexp.size();

//     for( i=0; i<n; i++){
//         pexp[i] = pexp[i]*Pfac;
//     }
    
    

//     static auto OptiminFunc = [&](double x[]) {
//         double param[] = {x[0], x[1]};
//         double fobj = 0.;
//         for (size_t i = 0; i< n; i++){
//             fobj += MinimizeLANGMUIR(pexp[i], nexp[i], param);
//         }
//         // printf("%f\n", fobj);
//         return fobj;
//     };
//     double (*MINFUN)(double *) = [](double x[]) { return OptiminFunc(x); };

//     double * new_param = simplex(MINFUN, ini, 4, 1e-6, .5, optconst);

//     // double new_param[] = {5.72399534e+00, 1.64655106e-04};

//     double * ncalc = (double *) malloc (n * sizeof (double));
//     double * pcalc =  (double *) malloc (n * sizeof (double));
    
//     #pragma omp parallel for
//     printf("P \t nexp \t ncalc\n");
//     for( i=0; i<n; i++){
//         pcalc[i] = pexp[i];
//         ncalc[i] = CalculateLangmuir(&pexp[i],new_param[0],new_param[1]);
//         printf("%f \t %f \t %f \n", pexp[i], nexp[i],  ncalc[i]);
//     }
    

//     gplot plt;
//     plt.plot(pcalc,ncalc , &n);
//     plt.show();


//     return 0;
// }
