#include "utils.h"

void CheckCompositionFraction(std::vector<double> &composition)
{
    double composition_sum = std::accumulate(composition.begin(), composition.end(), 0.0);
    // Check the factions consistency
    for (std::size_t i = 0; i < composition.size(); i++)
    {
        composition[i] = composition[i] / composition_sum;
    }
}

std::function<double(std::vector<double>, std::vector<double>)> GetDeviationInvoker(std::string DeviationEquation, double NumberOfParameters)
{

    if (DeviationEquation == "ARE")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetAREDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "SSE")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetSSEDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "EABS")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetEABSDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "RABS")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetRABSDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "HYBRID")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetHYBRIDDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "MPSD")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetMPSDDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "SORE")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetSOREDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "CHI_2")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetCHI2Deviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else if (DeviationEquation == "R_S")
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetRSDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
    else
    {
        return [=](std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings)
        {
            return GetAREDeviation(ExperimentalLoadings, CalculatedLoadings, NumberOfParameters);
        };
    }
};

void printString(std::string val)
{
    std::cout << val << std::endl;
}