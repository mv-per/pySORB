#include "_isotherms.h"

double freundlich(double Pressure, std::vector<double> Parameters)
{
    return (Parameters[0] * std::pow(Pressure, 1.0 / Parameters[1]));
}

double langmuir(double Pressure, std::vector<double> Parameters)
{
    return Parameters[0] * Parameters[1] * Pressure / (1.0 + Parameters[1] * Pressure);
}

double redlich_peterson(double Pressure, std::vector<double> Parameters)
{
    return (Parameters[0] * Pressure) / (1. + Parameters[1] * std::pow(Pressure, Parameters[2]));
}

double sips(double Pressure, std::vector<double> Parameters)
{

    return (Parameters[0] * std::pow(Parameters[1] * Pressure, 1. / Parameters[2])) /
           (1. + std::pow(Parameters[1] * Pressure, 1. / Parameters[2]));
}

double toth(double Pressure, std::vector<double> Parameters)
{

    return (Parameters[0] * Parameters[1] * Pressure) / std::pow(1. + std::pow(Parameters[1] * Pressure, Parameters[2]), 1. / Parameters[2]);
}

double unilan(double Pressure, std::vector<double> Parameters)
{
    double alpha = (1. + Parameters[2] * Parameters[3]) / (1.0 + Parameters[3] * Pressure);
    return Parameters[0] / 2. / Parameters[2] * std::log((1. + Parameters[1] * std::exp(Parameters[2]) * Pressure) / (1. + Parameters[1] * std::exp(-Parameters[2]) * Pressure));
}

double keller_staudt_toth(double Pressure, std::vector<double> Parameters)
{
    double alpha = (1. + Parameters[2] * Parameters[3]) / (1.0 + Parameters[3] * Pressure);
    return (Parameters[0] * Parameters[2] * Parameters[1] * Pressure / std::pow(1. + std::pow(Parameters[1] * Pressure, alpha), 1. / alpha));
}

double jensen_seaton(double Pressure, std::vector<double> Parameters)
{
    return Parameters[0] * Pressure * pow(1. + pow(Parameters[0] * Pressure / (Parameters[1] * (1. + Parameters[2] * Pressure)), Parameters[3]), -1. / Parameters[3]);
}
