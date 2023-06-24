#include "_base_isotherm_model.h"

double BaseIsothermModel::GetPureLoading(double Pressure, double Temperature, std::vector<double> Parameters)
{
    if (this->PureLoadingInvoker)
    {
        return this->PureLoadingInvoker(Pressure, Temperature, Parameters);
    }
    else
    {
        throw std::runtime_error("PureLoadingInvoker not defined");
    }
}

std::vector<double> BaseIsothermModel::GetPureLoadings(std::vector<double> Pressures, double Temperature, std::vector<double> Parameters)
{

    if (this->PureLoadingInvoker)
    {
        std::vector<double> CalculatedLoadings(Pressures.size(), 0.0);

        for (std::size_t i = 0; i < Pressures.size(); i++)
        {
            CalculatedLoadings[i] = this->PureLoadingInvoker(Pressures[i], Temperature, Parameters);
        }
        return CalculatedLoadings;
    }
    else
    {
        throw std::runtime_error("PureLoadingInvoker not defined");
    }
}

double BaseIsothermModel::GetDeviation(std::vector<double> Pressures, std::vector<double> ExperimentalLoadings, double Temperature, std::vector<double> Parameters, std::string DeviationEquation)
{
    std::vector<double> CalculatedLoadings = this->GetPureLoadings(Pressures, Temperature, Parameters);
    std::function<double(std::vector<double>, std::vector<double>)> DeviationInvoker = GetDeviationInvoker(DeviationEquation, Parameters.size());

    return DeviationInvoker(ExperimentalLoadings, CalculatedLoadings);
}

std::vector<double> BaseIsothermModel::GetMixtureLoading(double Pressure, double Temperature, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters)
{
    if (this->MixtureLoadingInvoker)
    {
        return this->MixtureLoadingInvoker(Pressure, Temperature, BulkComposition, Parameters);
    }
    else
    {
        throw std::runtime_error("MixtureLoadingInvoker not defined");
    }
}