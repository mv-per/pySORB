#include "pybind11/pybind11.h";
#include "pybind11/stl.h";
#include "src/classic_isotherm.h";

namespace py = pybind11;

PYBIND11_MODULE(ClassicIsotherm, m)
{
    py::class_<ClassicIsotherm>(m, "ClassicIsotherm")
        .def(py::init<std::string>())
        .def("GetLoading", &ClassicIsotherm::GetLoading, "Calculate loading",
             py::arg("Pressure"), py::arg("Temperature"), py::arg("Parameters"))
        .def("GetLoadings", &ClassicIsotherm::GetLoadings, "Calculate loadings",
             py::arg("Pressures"), py::arg("Temperature"), py::arg("Parameters"))
        .def("GetDeviation", &ClassicIsotherm::GetDeviation, "Calculate deviation",
             py::arg("Pressures"), py::arg("ExperimentalLoadings"), py::arg("Temperature"),
             py::arg("Parameters"), py::arg("DeviationEquation"));
}