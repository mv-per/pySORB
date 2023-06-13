#include "classic_isotherms.h"

std::function<double(double, std::vector<double>)> ClassicIsotherms::GetLoadingInvoker(std::string isotherm)
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
