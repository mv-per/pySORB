#include "classic_isotherm.h"

#include "deviation_functions.h"
#include "_isotherms.h"

ClassicIsotherm::ClassicIsotherm(std::string isotherm)
{
    this->IsothermInvoker = this->GetIsothermInvoker(isotherm);
}

double ClassicIsotherm::GetLoading(double Pressure, double Temperature, std::vector<double> Parameters)
{
    return this->IsothermInvoker(Pressure, Parameters);
}

std::vector<double> ClassicIsotherm::GetLoadings(std::vector<double> Pressures, double Temperature, std::vector<double> Parameters)
{
    std::vector<double> CalculatedLoadings(Pressures.size(), 0.0);

    for (std::size_t i = 0; i < Pressures.size(); i++)
    {
        CalculatedLoadings[i] = this->IsothermInvoker(Pressures[i], Parameters);
    }
    return CalculatedLoadings;
}

double ClassicIsotherm::GetDeviation(std::vector<double> Pressures, std::vector<double> ExperimentalLoadings, double Temperature, std::vector<double> Parameters, std::string DeviationEquation)
{
    std::vector<double> CalculatedLoadings = this->GetLoadings(Pressures, Temperature, Parameters);
    std::function<double(std::vector<double>, std::vector<double>)> DeviationInvoker = this->GetDeviationInvoker(DeviationEquation, Parameters.size());

    return DeviationInvoker(ExperimentalLoadings, CalculatedLoadings);
}

std::function<double(double, std::vector<double>)> ClassicIsotherm::GetIsothermInvoker(std::string isotherm)
{
    if (isotherm == "langmuir")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return langmuir(pressure, parameters);
        };
    }
    else if (isotherm == "sips")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return sips(pressure, parameters);
        };
    }
    else if (isotherm == "redlich-peterson")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return redlich_peterson(pressure, parameters);
        };
    }
    else if (isotherm == "toth")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return toth(pressure, parameters);
        };
    }
    else if (isotherm == "unilan")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return unilan(pressure, parameters);
        };
    }
    else if (isotherm == "freundlich")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return freundlich(pressure, parameters);
        };
    }
    else if (isotherm == "keller-staudt-toth")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return keller_staudt_toth(pressure, parameters);
        };
    }
    else
    {
        throw std::invalid_argument("Isotherm not found/defined.");
    }
}

std::function<double(std::vector<double>, std::vector<double>)> ClassicIsotherm::GetDeviationInvoker(std::string DeviationEquation, double NumberOfParameters)
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