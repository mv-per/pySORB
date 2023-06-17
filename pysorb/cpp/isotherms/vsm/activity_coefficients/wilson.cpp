#include "wilson.h"
/*
    Calculates the difference from the calculated and experimental pressure
    based on the parameters and n_exp
*/

double GetVSMWilsonPressure(double loading, double maximum_loading, double b, double a_1v, double a_v1)
{
    double theta = loading / maximum_loading;
    double Pressure = (maximum_loading / b * (theta / (1.0 - theta))) *
                      a_1v * ((1.0 - (1.0 - a_v1) * theta) / (a_1v + (1.0 - a_1v) * theta)) *
                      std::exp(-(a_v1 * (1.0 - a_v1) * theta) / (1.0 - (1.0 - a_v1) * theta) - ((1.0 - a_1v) * theta / (a_1v + (1.0 - a_1v) * theta)));
    return Pressure;
}

double GetVSMWilsonPressureEquilibrium(double loading, double Pressure, std::vector<double> parameters)
{
    double CalculatedPressure = GetVSMWilsonPressure(loading, parameters[0], parameters[1], parameters[2], parameters[3]);
    return 1.0e6 * (Pressure - CalculatedPressure) / Pressure;
}

/*
    Optimizes the parameters based on the pressure
*/
double
GetLoadingWILSON(double P, double T, std::vector<double> parameters)
{
    /* Makes a Initial estimate of the parameters using the Langmuir Equation */
    double loading_estimate = langmuir(P, {parameters[0], parameters[1]}) * .8;

    // Set the Equilibrium lambda function
    auto min = [&](double x)
    {
        return GetVSMWilsonPressureEquilibrium(std::fabs(x),
                                               P, parameters);
    };

    /* Solve using the Brent Method */
    double loading = brent_zeroin(min, loading_estimate, 1.0e-16);

    return loading;
}

// /*
//     Returns the absolute difference of n_calc and n_exp
// */
// double *CalculateLogGamma(double *xs, double **Aij, size_t *ncomp)
// {
//     /* Variable declaration */
//     double sum_1, sum_2, sum_3;
//     double *gamma = (double *)malloc(*ncomp + 1 * sizeof(double));
//     size_t i, j, k;

//     /* calculate gamma */
//     for (k = 0; k <= *ncomp; k++)
//     {
//         sum_1 = 0.;
//         for (j = 0; j <= *ncomp; j++)
//         {
//             sum_1 += xs[j] * Aij[k][j];
//         }
//         sum_3 = 0.;
//         for (i = 0; i <= *ncomp; i++)
//         {
//             sum_2 = 0.;
//             for (j = 0; j <= *ncomp; j++)
//             {
//                 sum_2 += xs[j] * Aij[i][j];
//             }

//             sum_3 += (xs[i] * Aij[i][k]) / sum_2;
//         }
//         gamma[k] = exp(1.0 - log(sum_1) - sum_3);
//     }
//     return gamma;
// }

// /*
//     Returns the binary coefficient interaction matrix
// */
// double **GetAijMatrix(vector<vector<double>> param, size_t *ncomp)
// {
//     /* Variable Declaration */
//     double **Aij = (double **)malloc(*ncomp + 1 * sizeof(double *));
//     size_t i, j;
//     for (i = 0; i <= *ncomp; i++)
//     {
//         Aij[i] = (double *)malloc(*ncomp + 1 * sizeof(double));
//     }
//     /* Obtain Aij matrix */

//     for (i = 0; i <= *ncomp; i++)
//     {
//         for (j = 0; j <= *ncomp; j++)
//         {
//             Aij[i][j] = 1.0;
//         }
//     }
//     for (i = 0; i < *ncomp; i++)
//     {
//         Aij[i][*ncomp] = param[i][2];
//         Aij[*ncomp][i] = param[i][3];
//     }
//     return Aij;
// }

// /*
//     Returns [x1,x2,x3,...,x_n, n_max]
// */
// void MinimizeVSMMixture_WILSON(int n, point_t *point, const void *arg)
// {
//     const mix_vsm_params *params = (const mix_vsm_params *)arg;

//     // Function computation
//     /* Variable declaration */
//     size_t i;
//     size_t ncomp = params->y.size();
//     double *x = (double *)malloc(ncomp * sizeof(double));
//     double *x_s = (double *)malloc(ncomp + 1 * sizeof(double));
//     double piA_RT, AdsorbedFugacity, ObjectiveFunction, n_max_mix, theta, nm, logGammaVacancy;
//     double **Aij;
//     double *gammaF;

//     /* Set adsorbed phase composition */
//     for (i = 0; i < ncomp; i++)
//     {
//         x[i] = point->x[i];
//     }
//     /* Calculates maximum adsorbed quantity */
//     nm = point->x[ncomp];

//     /* Calculate the n_max for the mixture */
//     n_max_mix = 0.0;
//     for (i = 0; i < ncomp; i++)
//     {
//         n_max_mix += x[i] * params->param[i][0];
//     }

//     /* Defines the "surface used" theta value */
//     theta = nm / n_max_mix;

//     /* Calculates the vacancy fraction */
//     for (i = 0; i < ncomp; i++)
//     {
//         x_s[i] = x[i] * theta;
//     }
//     x_s[ncomp] = 1.0 - theta;

//     Aij = GetAijMatrix(params->param, &ncomp);
//     gammaF = CalculateLogGamma(x_s, Aij, &ncomp);
//     logGammaVacancy = log(gammaF[ncomp] * x_s[ncomp]);

//     double *BulkFugacity = (double *)malloc(ncomp * sizeof(double));
//     for (i = 0; i < ncomp; i++)
//     {
//         BulkFugacity[i] = params->y[i] * params->P;
//     }

//     ObjectiveFunction = 0.0;
//     double sumx = 0.;

//     for (i = 0; i < ncomp; i++)
//     {
//         sumx += x[i];
//         piA_RT = ((params->param[i][0] - n_max_mix) / nm - 1.0) * logGammaVacancy;
//         AdsorbedFugacity = gammaF[i] * x[i] * (nm / n_max_mix) * (params->param[i][0] / params->param[i][1]) * params->param[i][2] * exp(params->param[i][3] - 1.0) * exp(piA_RT);
//         ObjectiveFunction += fabs(BulkFugacity[i] - AdsorbedFugacity) / BulkFugacity[i];
//     }
//     point->fx = ObjectiveFunction + 100. / (double)n * fabs(sumx - 1.0);
// }

// /*
//     Calculates the mixture adsorbed quantities
//     returns [n_1, n_2, n_3,..., n_n]
// */
// vector<double> CalculateMixtureWilsonVSM(double P, vector<double> y, vector<vector<double>> param)
// {
//     int n = y.size();

//     double *ini = CalculateExtendedLangmuir(param, y, &P, &n);

//     point_t start; // initial point
//     start.x = (double *)malloc(n + 1 * sizeof(double));
//     for (int i = 0; i < n + 1; i++)
//     {
//         start.x[i] = ini[i];
//     }

//     mix_vsm_params params;
//     params.P = P;
//     params.y = y;
//     params.param = param;

//     /* Set optimization settings for the Nelder-mead solver */
//     optimset_t optimset;
//     optimset.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
//     optimset.tolf = 1e-16;    // tolerance on the function value
//     optimset.max_iter = 1500; // maximum number of allowed iterations
//     optimset.max_eval = 1500; // maximum number of allowed function evaluations
//     optimset.verbose = 0;     // toggle verbose output during minimization

//     /* Defines solution pointer */
//     point_t solution;

//     /* Solve function with the nelder-mead solver */
//     nelder_mead(n + 1, &start, &solution, &MinimizeVSMMixture_WILSON, &params, &optimset);

//     /* Calculate the adsorbed quantity based on the maximum adsorbed quantity and
//     adsorbed molar fraction */
//     vector<double> result(n, 0.0);
//     for (int i = 0; i < n; i++)
//     {
//         result[i] = solution.x[i] * solution.x[n];
//     }

//     return result;
// }
