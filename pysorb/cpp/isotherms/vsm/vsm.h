#ifndef CLASSIC_ISOTHERM_H
#define CLASSIC_ISOTHERM_H

#include <string>
#include <vector>
#include <functional>
#include "activity_coefficients/nrtl.h"
#include "activity_coefficients/wilson.h"
#include "activity_coefficients/flory_huggins.h"
#include "../utils.h"
#include "../_base_isotherm_model.h"
// #include "deviation_functions.h"
// #include "_isotherms.h"

class VacancySolutionMethod : public BaseIsothermModel
{
public:
    /**
     * @brief Constructs a VacancySolutionMethod object and initializes the PureLoadingInvoker based on the given isotherm name.
     * @param model The name of the activity coefficient model to be used, options: `wilson`, `nrtl`, `flory-huggins`.
     */
    VacancySolutionMethod(std::string model);

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param model The name of the activity coefficient model to be used
     * @return The corresponding loading invoker function.
     * @throw std::invalid_argument If the model is not found or defined.
     */
    std::function<double(double, double, std::vector<double>)> GetPureLoadingInvoker(std::string model);
};

#endif