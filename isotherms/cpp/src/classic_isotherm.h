#ifndef CLASSIC_ISOTHERM_H
#define CLASSIC_ISOTHERM_H

#include <string>
#include <vector>
#include <functional>

class ClassicIsotherms
{
public:
    /**
     * @brief Constructs a ClassicIsotherms object and initializes the IsothermInvoker based on the given isotherm name.
     * @param isotherm The name of the isotherm to be used.
     */
    ClassicIsotherms(std::string isotherm);

    /**
     * @brief Calculates the loading for the given pressure, temperature, and parameters.
     * @param Pressure The pressure value.
     * @param Temperature The temperature value.
     * @param Parameters The parameters for the isotherm.
     * @return The calculated loading.
     */
    double GetLoading(double Pressure, double Temperature, std::vector<double> Parameters);

    /**
     * @brief Calculates the loadings for the given pressures, temperature, and parameters.
     * @param Pressures The vector of pressure values.
     * @param Temperature The temperature value.
     * @param Parameters The parameters for the isotherm.
     * @return The vector of calculated loadings.
     */
    std::vector<double> GetLoadings(std::vector<double> Pressures, double Temperature, std::vector<double> Parameters);

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

private:
    double Temperature;
    std::function<double(double, std::vector<double>)> IsothermInvoker;
    std::function<double(std::vector<double>, std::vector<double>)> DeviationInvoker;

    /**
     * @brief Returns the appropriate isotherm invoker function based on the given isotherm name.
     * @param isotherm The name of the isotherm.
     * @return The corresponding isotherm invoker function.
     * @throw std::invalid_argument If the isotherm is not found or defined.
     */
    std::function<double(double, std::vector<double>)> GetIsothermInvoker(std::string isotherm);

    /**
     * @brief Returns the appropriate deviation invoker function based on the given deviation function name.
     * @param DeviationEquation The name of the deviation function.
     * @param NumberOfParameters The number of parameters.
     * @return The corresponding deviation invoker function.
     */
    std::function<double(std::vector<double>, std::vector<double>)> GetDeviationInvoker(std::string DeviationEquation, double NumberOfParameters);
};

#endif