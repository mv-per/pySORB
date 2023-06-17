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
GetLoadingWILSON(double Pressure, std::vector<double> Parameters)
{
    /* Makes a Initial estimate of the parameters using the Langmuir Equation */
    double loading_estimate = langmuir(Pressure, {Parameters[0], Parameters[1]}) * .8;

    // Set the Equilibrium lambda function
    auto min = [&](double x)
    {
        return GetVSMWilsonPressureEquilibrium(std::fabs(x),
                                               Pressure, Parameters);
    };

    /* Solve using the Brent Method */
    double loading = brent_zeroin(min, loading_estimate, 1.0e-16);

    return loading;
}

/*
    Returns the absolute difference of n_calc and n_exp
*/
std::vector<double> CalculateLogGamma(std::vector<double> x_s, std::vector<std::vector<double>> Aij, std::size_t ncomp)
{
    /* Variable declaration */
    double sum_1, sum_2, sum_3;
    std::vector<double> gamma(ncomp + 1, 0.0);
    std::size_t i, j, k;

    /* calculate gamma */
    for (k = 0; k <= ncomp; k++)
    {
        sum_1 = 0.;
        for (j = 0; j <= ncomp; j++)
        {
            sum_1 += x_s[j] * Aij[k][j];
        }
        sum_3 = 0.;
        for (i = 0; i <= ncomp; i++)
        {
            sum_2 = 0.;
            for (j = 0; j <= ncomp; j++)
            {
                sum_2 += x_s[j] * Aij[i][j];
            }

            sum_3 += (x_s[i] * Aij[i][k]) / sum_2;
        }
        gamma[k] = std::exp(1.0 - std::log(sum_1) - sum_3);
    }
    return gamma;
}

/*
    Returns the binary coefficient interaction matrix
*/
std::vector<std::vector<double>> GetAijMatrix(std::vector<std::vector<double>> Parameters, std::size_t ncomp)
{
    /* Variable Declaration */
    std::vector<std::vector<double>> Aij(ncomp + 1, std::vector<double>(ncomp, 0.0));
    size_t i, j;

    /* Obtain Aij matrix */
    for (i = 0; i <= ncomp; i++)
    {
        for (j = 0; j <= ncomp; j++)
        {
            Aij[i][j] = 1.0;
        }
    }
    for (i = 0; i < ncomp; i++)
    {
        Aij[i][ncomp] = Parameters[i][2];
        Aij[ncomp][i] = Parameters[i][3];
    }
    return Aij;
}

/*
    Returns [x1,x2,x3,...,x_n, n_max]
*/
void MinimizeVSMMixture_WILSON(int n, point_t *point, const void *arg)
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

    /* Calculates the vacancy fraction */
    for (i = 0; i < ncomp; i++)
    {
        x_s[i] = x[i] * theta;
    }
    x_s[ncomp] = 1.0 - theta;

    std::vector<std::vector<double>> Aij = GetAijMatrix(optimization_arguments->Parameters, ncomp);
    std::vector<double> gammaF = CalculateLogGamma(x_s, Aij, ncomp);
    logGammaVacancy = std::log(gammaF[ncomp] * x_s[ncomp]);

    std::vector<double> BulkPartialPressures(ncomp, 0.0);
    for (i = 0; i < ncomp; i++)
    {
        BulkPartialPressures[i] = optimization_arguments->BulkComposition[i] * optimization_arguments->Pressure;
    }

    ObjectiveFunction = 0.0;
    double adsorbed_loading_ratio_sum = 0.;

    for (i = 0; i < ncomp; i++)
    {
        adsorbed_loading_ratio_sum += x[i];
        piA_RT = ((optimization_arguments->Parameters[i][0] - n_max_mix) / nm - 1.0) * logGammaVacancy;
        AdsorbedPartialPressure = gammaF[i] * x[i] * (nm / n_max_mix) * (optimization_arguments->Parameters[i][0] / optimization_arguments->Parameters[i][1]) * optimization_arguments->Parameters[i][2] * std::exp(optimization_arguments->Parameters[i][3] - 1.0) * exp(piA_RT);
        ObjectiveFunction += 1.0e6 * std::fabs(BulkPartialPressures[i] - AdsorbedPartialPressure) / BulkPartialPressures[i];
    }
    point->fx = ObjectiveFunction + 100. / (double)n * std::fabs(adsorbed_loading_ratio_sum - 1.0);
}

/*
    Calculates the mixture adsorbed quantities
    returns [n_1, n_2, n_3,..., n_n]
*/
std::vector<double> GetMixtureLoadingWILSON(double Pressure, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters)
{
    std::size_t ncomp = BulkComposition.size();

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

    MixtureOptimizationArguments OptimizationArguments;
    OptimizationArguments.Pressure = Pressure;
    OptimizationArguments.BulkComposition = BulkComposition;
    OptimizationArguments.Parameters = Parameters;

    /* Set optimization settings for the Nelder-mead solver */
    optimset_t optimset;
    optimset.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
    optimset.tolf = 1e-16;    // tolerance on the function value
    optimset.max_iter = 1500; // maximum number of allowed iterations
    optimset.max_eval = 1500; // maximum number of allowed function evaluations
    optimset.verbose = 0;     // toggle verbose output during minimization

    /* Defines solution pointer */
    point_t solution;

    /* Solve function with the nelder-mead solver */
    nelder_mead(ncomp + 1, &start, &solution, &MinimizeVSMMixture_WILSON, &OptimizationArguments, &optimset);

    /* Calculate the adsorbed quantity based on the maximum adsorbed quantity and
    adsorbed molar fraction */
    std::vector<double> result(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++)
    {
        result[i] = solution.x[i] * solution.x[ncomp];
    }

    return result;
}
