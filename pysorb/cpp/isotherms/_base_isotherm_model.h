

#ifndef BASE_ISOTHERM_MODEL_H
#define BASE_ISOTHERM_MODEL_H

#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include "deviation_functions.h"

#include "utils.h"

class BaseIsothermModel
{
public:
    /**
     * @brief Constructs a BaseIsothermModel object and initializes the IsothermInvoker based on the given isotherm name.
     * @param isotherm The name of the isotherm to be used.
     */
    BaseIsothermModel(){};

    virtual void SetupLoadingInvoker(std::string model) = 0;

    /**
     * @brief Calculates the loading for the given pressure, temperature, and parameters.
     * @param Pressure The pressure value.
     * @param Temperature The temperature value.
     * @param Parameters The parameters for the isotherm.
     * @return The calculated loading.
     */
    double GetPureLoading(double Pressure, double Temperature, std::vector<double> Parameters);

    /**
     * @brief Calculates the loadings for the given pressures, temperature, and parameters.
     * @param Pressures The vector of pressure values.
     * @param Temperature The temperature value.
     * @param Parameters The parameters for the isotherm.
     * @return The vector of calculated loadings.
     */
    std::vector<double> GetPureLoadings(std::vector<double> Pressures, double Temperature, std::vector<double> Parameters);

    /**
     * @brief Calculates the loading for the given pressure, temperature, and parameters.
     * @param Pressure The pressure value.
     * @param Temperature The temperature value.
     * @param Parameters The parameters for the isotherm.
     * @return The calculated loading.
     */
    std::vector<double> GetMixtureLoading(double Pressure, double Temperature, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters);

    /**
     * @brief Calculates the deviation between experimental loadings and calculated loadings using the specified deviation function.
     * @param Pressures The vector of pressure values.
     * @param ExperimentalLoadings The vector of experimental loadings.
     * @param Temperature The temperature value.
     * @param Parameters The parameters for the isotherm.
     * @param DeviationEquation The name of the deviation function to be used.
     * @return The calculated deviation value.
     */
    double GetDeviation(std::vector<double> Pressures, std::vector<double> ExperimentalLoadings, double Temperature, std::vector<double> Parameters, std::string DeviationEquation);

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param isotherm The name of the isotherm.
     * @return The corresponding isotherm invoker function.
     * @throw std::invalid_argument If the isotherm is not found or defined.
     */
    virtual std::function<double(double, double, std::vector<double>)> GetPureLoadingInvoker(std::string isotherm) = 0;

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param isotherm The name of the isotherm.
     * @return The corresponding isotherm invoker function.
     * @throw std::invalid_argument If the isotherm is not found or defined.
     */
    virtual std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> GetMixtureLoadingInvoker(std::string isotherm) = 0;

    std::function<double(double, double, std::vector<double>)> PureLoadingInvoker;
    std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> MixtureLoadingInvoker;
    std::vector<std::string> PureModels = {};
    std::vector<std::string> MixtureModels = {};
};

#endif