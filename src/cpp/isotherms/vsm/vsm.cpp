
#include "vsm.h"

void VacancySolutionMethod::SetupInvokers()
{
    if ((std::find(this->activity_models.begin(), this->activity_models.end(), this->Model)) != this->activity_models.end())
    {
        this->PureLoadingInvoker = this->GetPureLoadingInvoker(this->Model);
        this->MixtureLoadingInvoker = this->GetMixtureLoadingInvoker(this->Model);
    }
    else
    {
        std::runtime_error("Model not found.");
    }
}

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
            return GetLoadingWILSON(pressure, parameters);
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

std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> VacancySolutionMethod::GetMixtureLoadingInvoker(std::string model)
{
    if (model == "nrtl")
    {
        return [=](double pressure, double temperature, std::vector<double> bulk_composition, std::vector<std::vector<double>> parameters)
        {
            return GetMixtureLoadingNRTL(pressure, bulk_composition, parameters);
        };
    }
    else if (model == "wilson")
    {
        return [=](double pressure, double temperature, std::vector<double> bulk_composition, std::vector<std::vector<double>> parameters)
        {
            return GetMixtureLoadingWILSON(pressure, bulk_composition, parameters);
        };
    }
    else if (model == "flory-huggins")
    {
        return [=](double pressure, double temperature, std::vector<double> bulk_composition, std::vector<std::vector<double>> parameters)
        {
            return GetMixtureLoadingFloryHuggins(pressure, bulk_composition, parameters);
        };
    }
    else
    {
        throw std::invalid_argument("model not found/defined.");
    }
}