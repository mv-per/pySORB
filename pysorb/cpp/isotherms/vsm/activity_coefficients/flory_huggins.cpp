#include "flory_huggins.h"

typedef struct
{
    double n_max, b, a_1v;
    double P;
} pure_fh_vsm;

double OptimizePressureFloryHuggins(double loading, double P, double n_max, double b, double a_1v)
{
    double theta = loading / n_max;
    double Pcalc = (n_max / b * (theta / (1.0 - theta))) *
                   std::exp(a_1v * a_1v * theta / (1.0 + a_1v * theta));
    return 1.0e6 * (P - Pcalc) / P;
}

double GetLoadingFloryHuggins(double P, double T, std::vector<double> parameters)
{
    // Initial loading estimate using langmuir;
    double loading_estimate = langmuir(P, {parameters[0], parameters[1]}) * .8;

    // Set the Equilibrium lambda function
    auto min = [&](double x)
    {
        return OptimizePressureFloryHuggins(std::fabs(x),
                                            P, parameters[0], parameters[1],
                                            parameters[2]);
    };

    /* Solve using the Brent Method */
    double loading = brent_zeroin(min, loading_estimate, 1.0e-16);

    return loading;
    // point_t start_pure;
    // start_pure.x = (double *)malloc(1 * sizeof(double));
    // start_pure.x[0] = loading_estimate;

    // pure_fh_vsm params_pure;
    // params_pure.n_max = parameters[0];
    // params_pure.b = parameters[1];
    // params_pure.a_1v = parameters[2];
    // params_pure.P = P;

    // // optimisation settings
    // optimset_t optimset_pure;
    // optimset_pure.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
    // optimset_pure.tolf = 1e-16;    // tolerance on the function value
    // optimset_pure.max_iter = 1500; // maximum number of allowed iterations
    // optimset_pure.max_eval = 1500; // maximum number of allowed function evaluations
    // optimset_pure.verbose = 0;     // toggle verbose output during minimization
    // point_t solution_pure;
    // nelder_mead(1, &start_pure, &solution_pure, &FHPcalc, &params_pure, &optimset_pure);

    // return solution_pure.x[0];
}

// double MinimizeFH(double P, double n_exp, double param[])
// {
//     // Minimiza
//     double n_calc = CalculateFH(P, param);
//     return fabs(n_exp - n_calc);
// }

std::vector<double> CalculateLogGammaFH(std::vector<double> xs, std::vector<std::vector<double>> Aij, size_t ncomp)
{
    /* Variable declaration */
    double sum_1;
    std::vector<double> gamma(ncomp + 1, 0.0);
    size_t i, j;

    /* calculate gamma */
    for (i = 0; i <= ncomp; i++)
    {
        sum_1 = 0.;
        for (j = 0; j <= ncomp; j++)
        {
            sum_1 += xs[j] / (Aij[i][j] + 1.0);
        }

        gamma[i] = exp((1.0 - 1.0 / sum_1) - log(sum_1));
    }
    return gamma;
}

std::vector<std::vector<double>> GetAijMatrixFH(std::vector<std::vector<double>> param, std::size_t ncomp)
{
    /* Variable Declaration */
    std::vector<std::vector<double>> Aij(ncomp + 1, std::vector<double>(ncomp + 1, 0.0));
    size_t i, j;

    for (i = 0; i < ncomp; i++)
    {
        for (j = 0; j < ncomp; j++)
        {
            if (i == j)
            {
                Aij[i][j] = 0.0;
            }
            else
            {
                Aij[i][j] = (param[i][2] + 1.0) / (param[j][2] + 1.0) - 1.0;
            }
        }
    }
    for (i = 0; i < ncomp; i++)
    {
        Aij[i][ncomp] = param[i][2];
        Aij[ncomp][i] = param[i][2];
    }
    Aij[ncomp][ncomp] = 0.0;
    return Aij;
}

void MinimizeFHVSMMixture(int n, point_t *point, const void *arg)
{
    const MixtureOptimizationArguments *optimization_arguments = (const MixtureOptimizationArguments *)arg;

    // Function computation
    /* Variable declaration */
    size_t i;
    size_t ncomp = optimization_arguments->BulkComposition.size();
    std::vector<double> x(ncomp, 0.0);
    std::vector<double> x_s(ncomp, 0.0);
    double piA_RT, AdsorbedPartialPressure, ObjectiveFunction, n_max_mix, theta, nm, logGammaVacancy;

    /* Set adsorbed phase composition */
    for (i = 0; i < ncomp; i++)
    {
        x[i] = point->x[i];
    }
    /* Calculates maximum adsorbed quantity */
    nm = point->x[ncomp];

    /* Calculate the n_max for the mixture */
    n_max_mix = 0.0;
    for (i = 0; i < ncomp; i++)
    {
        n_max_mix += x[i] * optimization_arguments->Parameters[i][0];
    }

    /* Defines the "surface used" theta value */
    theta = nm / n_max_mix;

    for (i = 0; i < ncomp; i++)
    {
        x_s[i] = x[i] * theta;
    }
    x_s[ncomp] = 1.0 - theta;

    std::vector<std::vector<double>> Aij = GetAijMatrixFH(optimization_arguments->Parameters, ncomp);

    std::vector<double> gammaF = CalculateLogGammaFH(x_s, Aij, ncomp);

    std::vector<double> PartialPressures(ncomp, 0.0);
    for (i = 0; i < ncomp; i++)
    {
        PartialPressures[i] = optimization_arguments->BulkComposition[i] * optimization_arguments->Pressure;
    }

    ObjectiveFunction = 0.0;
    double sumx = 0.;
    // double CompositionObjectiveFunction = 0.0;
    for (i = 0; i < ncomp; i++)
    {
        sumx += x[i];
        piA_RT = ((optimization_arguments->Parameters[i][0] - n_max_mix) / nm - 1.0) * log(gammaF[ncomp] * x_s[ncomp]);
        AdsorbedPartialPressure = gammaF[i] * x[i] * (nm / n_max_mix) * (optimization_arguments->Parameters[i][0] / optimization_arguments->Parameters[i][1]) * (exp(optimization_arguments->Parameters[i][2]) / (1.0 + optimization_arguments->Parameters[i][2])) * exp(piA_RT);
        ObjectiveFunction += std::fabs(PartialPressures[i] - AdsorbedPartialPressure) / PartialPressures[i];
    }
    // printf("fx = %f", ObjectiveFunction + fabs(sumx-1.0)/sumx);
    point->fx = 1.0e6 * ObjectiveFunction + 100. / (double)n * std::fabs(sumx - 1.0);
}

std::vector<double> GetMixtureLoadingFloryHuggins(double Pressure, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters)
{
    int ncomp = BulkComposition.size();

    point_t start; // initial point
    std::vector<double> initial_estimates = extended_langmuir(Pressure, BulkComposition, Parameters);

    start.x = (double *)malloc(ncomp + 1 * sizeof(double));
    double sum_n = 0;
    for (int i = 0; i < ncomp; i++)
    {
        start.x[i] = initial_estimates[i] * 0.8;
        sum_n += initial_estimates[i];
    }
    start.x[ncomp] = sum_n;

    MixtureOptimizationArguments params;
    params.Pressure = Pressure;
    params.BulkComposition = BulkComposition;
    params.Parameters = Parameters;

    // optimisation settings
    optimset_t optimset;
    optimset.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
    optimset.tolf = 1e-16;    // tolerance on the function value
    optimset.max_iter = 1500; // maximum number of allowed iterations
    optimset.max_eval = 1500; // maximum number of allowed function evaluations
    optimset.verbose = 0;     // toggle verbose output during minimization

    point_t solution;
    nelder_mead(ncomp + 1, &start, &solution, &MinimizeFHVSMMixture, &params, &optimset);

    std::vector<double> result(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++)
    {
        result[i] = solution.x[i] * solution.x[ncomp];
    }

    return result;
}
