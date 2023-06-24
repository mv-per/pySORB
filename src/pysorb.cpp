#include <pybind11/pybind11.h>

#include "cpp/isotherms/_base_isotherm_module.h"
#include "cpp/isotherms/data_classes_module.h"
#include "cpp/isotherms/empirical/classic_isotherms_module.h"
#include "cpp/isotherms/pta/pta_module.h"
#include "cpp/isotherms/vsm/vsm_module.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(pysorb, m) {
  m.doc() = R"pbdoc(
        PySorb API
    )pbdoc";

  bindDataClasses(m);
  bindBaseIsothermModel(m);
  bindVacancySolutionMethod(m);
  bindPotentialTheoryModels(m);
  bindEmpiricalIsotherms(m);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}