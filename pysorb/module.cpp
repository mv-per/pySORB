#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "cpp/isotherms/empirical/classic_isotherms.h"

namespace py = pybind11;

PYBIND11_MODULE(pysorb, m)
{
     py::class_<ClassicIsotherms>(m, "ClassicIsotherms")
         .def(py::init<std::string>())
         .def("get_loading", &ClassicIsotherms::GetPureLoading, "Calculate loading",
              py::arg("Pressure"), py::arg("Temperature"), py::arg("Parameters"))
         .def("get_loadings", &ClassicIsotherms::GetPureLoadings, "Calculate loadings",
              py::arg("Pressures"), py::arg("Temperature"), py::arg("Parameters"))
         .def("get_deviation", &ClassicIsotherms::GetDeviation, "Calculate deviation",
              py::arg("Pressures"), py::arg("ExperimentalLoadings"), py::arg("Temperature"),
              py::arg("Parameters"), py::arg("DeviationEquation"));
}