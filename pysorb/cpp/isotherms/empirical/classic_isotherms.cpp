#include "classic_isotherms.h"

std::function<double(double, double, std::vector<double>)>
ClassicIsotherms::GetPureLoadingInvoker(std::string isotherm)
{
    if (isotherm == "langmuir")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return langmuir(pressure, parameters);
        };
    }
    else if (isotherm == "dual-langmuir")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return dual_langmuir(pressure, parameters);
        };
    }
    else if (isotherm == "sips")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return sips(pressure, parameters);
        };
    }
    else if (isotherm == "redlich-peterson")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return redlich_peterson(pressure, parameters);
        };
    }
    else if (isotherm == "toth")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return toth(pressure, parameters);
        };
    }
    else if (isotherm == "unilan")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return unilan(pressure, parameters);
        };
    }
    else if (isotherm == "freundlich")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return freundlich(pressure, parameters);
        };
    }
    else if (isotherm == "freundlich-2")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return freundlich_2(pressure, temperature, parameters);
        };
    }
    else if (isotherm == "keller-staudt-toth")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return keller_staudt_toth(pressure, parameters);
        };
    }
    else
    {
        throw std::invalid_argument("Isotherm not found/defined in the pure isotherms. If you want to calculate mixture please set pure=false.");
    }
}

std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> ClassicIsotherms::GetMixtureLoadingInvoker(std::string isotherm)
{
    if (isotherm == "extended-langmuir")
    {
        return [=](double pressure, double temperature, std::vector<double> bulk_composition, std::vector<std::vector<double>> parameters)
        {
            return extended_langmuir(pressure, bulk_composition, parameters);
        };
    }
    else
    {
        throw std::invalid_argument("Isotherm not found/defined in the Mixture isotherms.");
    }
}
