#include "_base_isotherm_model.h"

BaseIsothermModel::BaseIsothermModel(std::string isotherm)
{
    if ((std::find(this->PureModels.begin(), this->PureModels.end(), isotherm)) != this->PureModels.end())
    {
        this->PureLoadingInvoker = this->GetPureLoadingInvoker(isotherm);
    }
    else if ((std::find(this->MixtureModels.begin(), this->MixtureModels.end(), isotherm)) != this->MixtureModels.end())
    {
        this->MixtureLoadingInvoker = this->GetMixtureLoadingInvoker(isotherm);
    }
    else
    {
        std::invalid_argument("Isotherm not found.");
    }
}

double BaseIsothermModel::GetPureLoading(double Pressure, double Temperature, std::vector<double> Parameters)
{
    return this->PureLoadingInvoker(Pressure, Temperature, Parameters);
}

std::vector<double> BaseIsothermModel::GetPureLoadings(std::vector<double> Pressures, double Temperature, std::vector<double> Parameters)
{
    std::vector<double> CalculatedLoadings(Pressures.size(), 0.0);

    for (std::size_t i = 0; i < Pressures.size(); i++)
    {
        CalculatedLoadings[i] = this->PureLoadingInvoker(Pressures[i], Temperature, Parameters);
    }
    return CalculatedLoadings;
}

double BaseIsothermModel::GetDeviation(std::vector<double> Pressures, std::vector<double> ExperimentalLoadings, double Temperature, std::vector<double> Parameters, std::string DeviationEquation)
{
    std::vector<double> CalculatedLoadings = this->GetPureLoadings(Pressures, Temperature, Parameters);
    std::function<double(std::vector<double>, std::vector<double>)> DeviationInvoker = GetDeviationInvoker(DeviationEquation, Parameters.size());

    return DeviationInvoker(ExperimentalLoadings, CalculatedLoadings);
}