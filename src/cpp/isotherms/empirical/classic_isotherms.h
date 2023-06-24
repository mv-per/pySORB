#ifndef EMPIRICAL_ISOTHERM_H
#define EMPIRICAL_ISOTHERM_H

#include "../_base_isotherm_model.h"
#include "../deviation_functions.h"
#include "../utils.h"
#include "files/_isotherms.h"
#include "files/extended_isotherms.h"
#include <functional>
#include <string>
#include <vector>

class EmpiricalIsotherms : public BaseIsothermModel {
public:
  /**
   * @brief Constructs a EmpiricalIsotherms object and initializes the
   * IsothermInvoker based on the given isotherm name.
   * @param isotherm The name of the isotherm to be used.
   */
  EmpiricalIsotherms(std::string isotherm) {
    this->Isotherm = isotherm;
    this->SetupInvokers();
  }

  void SetupInvokers() override;

  /**
   * @brief Returns the appropriate loading invoker function based on the given
   * model name.
   * @param isotherm The name of the isotherm.
   * @return The corresponding isotherm invoker function.
   * @throw std::invalid_argument If the isotherm is not found or defined.
   */
  std::function<double(double, double, std::vector<double>)>
  GetPureLoadingInvoker(std::string isotherm) override;

  /**
   * @brief Returns the appropriate loading invoker function based on the given
   * model name.
   * @param isotherm The name of the isotherm.
   * @return The corresponding isotherm invoker function.
   * @throw std::invalid_argument If the isotherm is not found or defined.
   */
  std::function<std::vector<double>(double, double, std::vector<double>,
                                    std::vector<std::vector<double>>)>
  GetMixtureLoadingInvoker(std::string isotherm) override;

private:
  std::string Isotherm;
  std::vector<std::string> PureModels = {
      "langmuir",           "dual-langmuir", "sips",         "toth",
      "redlich-peterson",   "unilan",        "freundlich-2", "freundlich",
      "keller-staudt-toth", "jensen-seaton",
  };

  std::vector<std::string> MixtureModels = {
      "extended-langmuir",
  };
};

#endif