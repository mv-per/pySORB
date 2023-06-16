
#include "vsm.h"

std::function<double(double, double, std::vector<double>)> VacancySolutionMethod::GetPureLoadingInvoker(std::string model)
{

    if (model == "nrtl")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return GetLoadingNRTL(pressure, temperature, parameters);
        };
    }
    else if (model == "wilson")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return GetLoadingWILSON(pressure, temperature, parameters);
        };
    }
    else if (model == "flory-huggins")
    {
        return [=](double pressure, double temperature, std::vector<double> parameters)
        {
            return GetLoadingFloryHuggins(pressure, temperature, parameters);
        };
    }
    else
    {
        throw std::invalid_argument("model not found/defined.");
    }
}