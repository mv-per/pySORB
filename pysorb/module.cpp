#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include "cpp/isotherms/empirical/classic_isotherms.h"
#include "cpp/isotherms/pta/pta.h"
#include "cpp/isotherms/vsm/vsm.h"
#include "cpp/isotherms/data_classes.h"
#include "cpp/isotherms/_base_isotherm_model.h"

namespace py = pybind11;

PYBIND11_MODULE(pysorb, m)
{

    m.doc() = "C/C++ methods to run the PTA model";

    py::class_<BaseIsothermModel>(m, "BaseIsothermModel")
        .def("get_loading", &BaseIsothermModel::GetPureLoading, "Calculate loading", py::arg("Pressure"), py::arg("Temperature"), py::arg("Parameters"))
        .def("get_loadings", &BaseIsothermModel::GetPureLoadings, "Calculate loadings", py::arg("Pressures"), py::arg("Temperature"), py::arg("Parameters"))
        .def("get_deviation", &BaseIsothermModel::GetDeviation, "Calculate deviation", py::arg("Pressures"), py::arg("ExperimentalLoadings"), py::arg("Temperature"), py::arg("Parameters"), py::arg("DeviationEquation"))
        .def("get_mixture_loading", &BaseIsothermModel::GetMixtureLoading, py::arg("Pressure"), py::arg("Temperature"), py::arg("BulkComposition"), py::arg("Parameters"));

    py::class_<Adsorbent>(m, "Adsorbent")
        .def(py::init<>())
        .def_readwrite("name", &Adsorbent::Name)
        .def_readwrite("diameter", &Adsorbent::SolidDiameter)
        .def_readwrite("atomic_density", &Adsorbent::SolidAtomicDensity)
        .def(py::init([](std::string name, double solid_diameter, double solid_atomic_density)
                      { return new Adsorbent(name, solid_diameter, solid_atomic_density); }),
             py::arg("name"),
             py::arg("diameter"),
             py::arg("atomic_density"));

    py::class_<Fluid>(m, "Fluid")
        .def_readwrite("name", &Fluid::Name)
        .def_readwrite("critical_pressure", &Fluid::CriticalPressure)
        .def_readwrite("critical_temperature", &Fluid::CriticalTemperature)
        .def_readwrite("accentric_factor", &Fluid::AccentricFactor)
        .def_readwrite("critical_compressibility", &Fluid::CriticalCompressibility)
        .def_readwrite("lennard_jonnes_diameter", &Fluid::LennardJonnesDiameter)
        .def(py::init([](std::string name,
                         double critical_pressure,
                         double critical_temperature,
                         double accentric_factor,
                         double critical_compressibility,
                         double lennard_jonnes_diameter)
                      { return new Fluid(name, critical_pressure, critical_temperature, accentric_factor, critical_compressibility, lennard_jonnes_diameter); }),
             py::arg("name"),
             py::arg("critical_pressure"),
             py::arg("critical_temperature"),
             py::arg("accentric_factor"),
             py::arg("critical_compressibility"),
             py::arg("lennard_jonnes_diameter"))
        .def(py::init([](std::string name,
                         double critical_pressure,
                         double critical_temperature,
                         double accentric_factor)
                      { return new Fluid(name, critical_pressure, critical_temperature, accentric_factor); }),
             py::arg("name"),
             py::arg("critical_pressure"),
             py::arg("critical_temperature"),
             py::arg("accentric_factor"))
        .def(py::init<>());

    py::class_<EmpiricalIsotherms, BaseIsothermModel>(m, "EmpiricalIsotherms")
        .def(py::init<std::string>(),
             py::arg("model"));

    py::class_<VacancySolutionMethod, BaseIsothermModel>(m, "VacancySolutionMethod")
        .def(py::init<std::string>(),
             py::arg("model"));

    py::class_<PotentialTheoryModels, BaseIsothermModel>(m, "PotentialTheoryModels")
        .def(py::init([](std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, std::vector<Fluid> fluids, Adsorbent adsorbent)
                      { return new PotentialTheoryModels(potential, equation_of_state, isotherm_type, num_of_layers, fluids, adsorbent); }),
             py::arg("potential"),
             py::arg("equation_of_state"),
             py::arg("isotherm_type"),
             py::arg("num_of_layers"),
             py::arg("fluids"),
             py::arg("adsorbent"))
        .def(py::init([](std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, std::vector<Fluid> fluids)
                      { return new PotentialTheoryModels(potential, equation_of_state, isotherm_type, num_of_layers, fluids); }),
             py::arg("potential"),
             py::arg("equation_of_state"),
             py::arg("isotherm_type"),
             py::arg("num_of_layers"),
             py::arg("fluids"))
        .def(py::init([](std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, Fluid fluid, Adsorbent adsorbent)
                      { return new PotentialTheoryModels(potential, equation_of_state, isotherm_type, num_of_layers, fluid, adsorbent); }),
             py::arg("potential"),
             py::arg("equation_of_state"),
             py::arg("isotherm_type"),
             py::arg("num_of_layers"),
             py::arg("fluid"),
             py::arg("adsorbent"))
        .def(py::init([](std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, Fluid fluid)
                      { return new PotentialTheoryModels(potential, equation_of_state, isotherm_type, num_of_layers, fluid); }),
             py::arg("potential"),
             py::arg("equation_of_state"),
             py::arg("isotherm_type"),
             py::arg("num_of_layers"),
             py::arg("fluid"));
}