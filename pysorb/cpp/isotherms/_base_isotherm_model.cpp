#include "_base_isotherm_model.h"

BaseIsothermModel::BaseIsothermModel(std::string isotherm)
{
    this->LoadingInvoker = this->GetLoadingInvoker(isotherm);
}

double BaseIsothermModel::GetLoading(double Pressure, double Temperature, std::vector<double> Parameters)
{
    return this->LoadingInvoker(Pressure, Temperature, Parameters);
}

std::vector<double> BaseIsothermModel::GetLoadings(std::vector<double> Pressures, double Temperature, std::vector<double> Parameters)
{
    std::vector<double> CalculatedLoadings(Pressures.size(), 0.0);

    for (std::size_t i = 0; i < Pressures.size(); i++)
    {
        CalculatedLoadings[i] = this->LoadingInvoker(Pressures[i], Temperature, Parameters);
    }
    return CalculatedLoadings;
}

double BaseIsothermModel::GetDeviation(std::vector<double> Pressures, std::vector<double> ExperimentalLoadings, double Temperature, std::vector<double> Parameters, std::string DeviationEquation)
{
    std::vector<double> CalculatedLoadings = this->GetLoadings(Pressures, Temperature, Parameters);
    std::function<double(std::vector<double>, std::vector<double>)> DeviationInvoker = GetDeviationInvoker(DeviationEquation, Parameters.size());

    return DeviationInvoker(ExperimentalLoadings, CalculatedLoadings);
}