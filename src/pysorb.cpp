#include <pybind11/pybind11.h>

#include "cpp/isotherms/_base_isotherm_model.h"
#include "cpp/isotherms/data_classes.h"
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

  bindVacancySolutionMethod(m);
  bindPotentialTheoryModels(m);
  bindEmpiricalIsotherms(m);

  //   py::class_<BaseIsothermModel>(m, "BaseIsothermModel")
  //       .def("get_loading", &BaseIsothermModel::GetPureLoading,
  //       "TESTEEEEEEEEE",
  //            py::arg("Pressure"), py::arg("Temperature"),
  //            py::arg("Parameters"))
  //       .def("get_loadings", &BaseIsothermModel::GetPureLoadings,
  //            py::arg("Pressures"), py::arg("Temperature"),
  //            py::arg("Parameters"))
  //       .def("get_deviation", &BaseIsothermModel::GetDeviation,
  //            py::arg("Pressures"), py::arg("ExperimentalLoadings"),
  //            py::arg("Temperature"), py::arg("Parameters"),
  //            py::arg("DeviationEquation"))
  //       .def("get_mixture_loading", &BaseIsothermModel::GetMixtureLoading,
  //            py::arg("Pressure"), py::arg("Temperature"),
  //            py::arg("BulkComposition"), py::arg("Parameters"));

  py::class_<Adsorbent>(m, "Adsorbent")
      .def(py::init<>())
      .def(py::init([](std::string name, double solid_diameter,
                       double solid_atomic_density) {
             return new Adsorbent(name, solid_diameter, solid_atomic_density);
           }),
           py::arg("name"), py::arg("diameter"), py::arg("atomic_density"),
           R"pobdoc(
            Returns the adsorbent  with the given name.

           Args:
               name (str): The name of the adsorbent
               diameter (float):  Diameters of the adsorbent wall [unit]
               atomic_density (float):  Atomic density of the adsorbent wall [unit]

            Returns:
               The Adsorbent
        )pobdoc");

  py::class_<Fluid>(m, "Fluid")
      .def_readwrite("name", &Fluid::Name)
      .def_readwrite("critical_pressure", &Fluid::CriticalPressure)
      .def_readwrite("critical_temperature", &Fluid::CriticalTemperature)
      .def_readwrite("accentric_factor", &Fluid::AccentricFactor)
      .def_readwrite("critical_compressibility",
                     &Fluid::CriticalCompressibility)
      .def_readwrite("lennard_jonnes_diameter", &Fluid::LennardJonnesDiameter)
      .def(py::init([](std::string name, double critical_pressure,
                       double critical_temperature, double accentric_factor,
                       double critical_compressibility,
                       double lennard_jonnes_diameter) {
             return new Fluid(name, critical_pressure, critical_temperature,
                              accentric_factor, critical_compressibility,
                              lennard_jonnes_diameter);
           }),
           py::arg("name"), py::arg("critical_pressure"),
           py::arg("critical_temperature"), py::arg("accentric_factor"),
           py::arg("critical_compressibility"),
           py::arg("lennard_jonnes_diameter"), R"(
           Returns the fluid with a specific name.

           Args:
               name (str): The name of the fluid
               critical_pressure (float): The fluid critical Pressure [Pa]
               critical_temperature (float): The fluid critical Temperature [K]
               accentric_factor (float): The fluid accentric factor [-]
               critical_compressibility (float): The fluid critical compressibility [-]
               lennard_jonnes_diameter (float): The fluid Lennard-Jonnes diameter [A]

           Returns:
               The Fluid

           )")
      .def(py::init([](std::string name, double critical_pressure,
                       double critical_temperature, double accentric_factor) {
             return new Fluid(name, critical_pressure, critical_temperature,
                              accentric_factor);
           }),
           py::arg("name"), py::arg("critical_pressure"),
           py::arg("critical_temperature"), py::arg("accentric_factor"), R"(

           Returns the species with the given name.

           Args:
               name (str): The name of the fluid
               critical_pressure (float): The fluid critical Pressure [Pa]
               critical_temperature (float): The fluid critical Temperature [K]
               accentric_factor (float): The fluid accentric factor [-]

           Returns:
               The Fluid

           )")
      .def(py::init<>());

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}