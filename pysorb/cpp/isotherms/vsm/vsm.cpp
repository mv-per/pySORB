
#include "vsm.h"

// VacancySolutionMethod::VacancySolutionMethod(std::string model)
// {
//     this->LoadingInvoker = this->GetLoadingInvoker(model);
// }

std::function<double(double, std::vector<double>)> VacancySolutionMethod::GetLoadingInvoker(std::string model)
{

    if (model == "nrtl")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return GetLoadingNRTL(pressure, parameters);
        };
    }
    else if (model == "wilson")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return GetLoadingWILSON(pressure, parameters);
        };
    }
    else if (model == "flory-huggins")
    {
        return [=](double pressure, std::vector<double> parameters)
        {
            return GetLoadingFloryHuggins(pressure, parameters);
        };
    }
    else
    {
        throw std::invalid_argument("model not found/defined.");
    }
}