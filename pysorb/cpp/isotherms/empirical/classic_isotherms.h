#ifndef CLASSIC_ISOTHERM_H
#define CLASSIC_ISOTHERM_H

#include <string>
#include <vector>
#include <functional>
#include "../deviation_functions.h"
#include "_isotherms.h"
#include "extended_isotherms.h"
#include "../utils.h"
#include "../_base_isotherm_model.h"

class ClassicIsotherms : public BaseIsothermModel
{
public:
    /**
     * @brief Constructs a ClassicIsotherms object and initializes the IsothermInvoker based on the given isotherm name.
     * @param isotherm The name of the isotherm to be used.
     */
    ClassicIsotherms(std::string isotherm);

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param isotherm The name of the isotherm.
     * @return The corresponding isotherm invoker function.
     * @throw std::invalid_argument If the isotherm is not found or defined.
     */
    std::function<double(double, double, std::vector<double>)> GetPureLoadingInvoker(std::string isotherm);

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param isotherm The name of the isotherm.
     * @return The corresponding isotherm invoker function.
     * @throw std::invalid_argument If the isotherm is not found or defined.
     */
    std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> GetMixtureLoadingInvoker(std::string isotherm);

private:
    std::vector<std::string> PureModels =
        {
            "langmuir",
            "dual-langmuir",
            "sips",
            "toth",
            "redlich-peterson",
            "unilan",
            "freundlich-2",
            "freundlich",
            "keller-staudt-toth",
        };

    std::vector<std::string> MixtureModels =
        {
            "extended-langmuir",
        };
};

#endif