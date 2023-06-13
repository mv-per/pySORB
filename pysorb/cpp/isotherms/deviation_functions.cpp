#ifndef ERRORS_H
#define ERRORS_H

/*
    Sources:
    1 - 10.1016/j.cej.2009.09.013
    2 - 10.1006/jcis.2002.8664
    3 - 10.1021/la800725s
*/

#include "deviation_functions.h"

double GetSSEDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += std::pow(CalculatedLoadings[i] - ExperimentalLoadings[i], 2.0);
    }

    return deviation;
}

double GetAREDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += (ExperimentalLoadings[i] - CalculatedLoadings[i]) / ExperimentalLoadings[i];
    }

    return 100.0 / double(ExperimentalLoadings.size()) * deviation;
}

double GetEABSDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += std::fabs(ExperimentalLoadings[i] - CalculatedLoadings[i]);
    }

    return deviation;
}

double GetHYBRIDDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += std::pow((ExperimentalLoadings[i] - CalculatedLoadings[i]) / ExperimentalLoadings[i], 2.0);
    }

    return 100. / (double(ExperimentalLoadings.size()) - double(NumberOfParameters)) * deviation;
}

double GetMPSDDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += std::pow(ExperimentalLoadings[i] - CalculatedLoadings[i], 2.0) / ExperimentalLoadings[i];
    }

    return 100. * std::sqrt(1.0 / (double(ExperimentalLoadings.size()) - NumberOfParameters) * deviation);
}

double GetSOREDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += std::pow((ExperimentalLoadings[i] - CalculatedLoadings[i]) / ExperimentalLoadings[i] * GetAREDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters), 2.0);
    }

    return std::sqrt(deviation / (double(ExperimentalLoadings.size()) - 1.0));
}

double GetCHI2Deviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 0;

    for (std::size_t i = 0; i < ExperimentalLoadings.size(); i++)
    {
        deviation += std::pow(ExperimentalLoadings[i] - CalculatedLoadings[i], 2.0) / ExperimentalLoadings[i];
    }

    return deviation;
}

double GetRSDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters)
{
    double deviation = 1.0 - 6.0 * GetSSEDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters) / double(ExperimentalLoadings.size()) / std::pow(double(ExperimentalLoadings.size()) - 1.0, 2.0);

    return deviation;
}

#endif