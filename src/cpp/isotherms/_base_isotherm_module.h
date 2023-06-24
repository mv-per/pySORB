#ifndef BASE_ISOTHERM_MODULE_H
#define BASE_ISOTHERM_MODULE_H

#include <pybind11/pybind11.h>

#include "_base_isotherm_model.h"

void bindBaseIsothermModel(pybind11::module &m) {

  pybind11::class_<BaseIsothermModel>(m, "BaseIsothermModel")
      .def("get_loading", &BaseIsothermModel::GetPureLoading, R"(
        Calculates the loading for the given pressure, temperature, and parameters.

        :param pressure: The pressure value [Pa].
        :type pressure: float
        :param temperature: The temperature value [K].
        :type temperature: float
        :param parameters: The parameters for the isotherm.
        :type parameters: List[float]
        :return: The calculated loading [mol/kg].
        :rtype: float
      )",
           pybind11::arg("pressure"), pybind11::arg("temperature"),
           pybind11::arg("parameters"))
      .def("get_loadings", &BaseIsothermModel::GetPureLoadings,
           pybind11::arg("pressures"), pybind11::arg("temperature"),
           pybind11::arg("parameters"), R"(
            Calculates the loadings for the given pressures, temperature, and parameters.

            :param pressures: The vector of pressure values [Pa].
            :type pressures: List[float]
            :param temperature: The temperature value [K].
            :type temperature: float
            :param parameters: The parameters for the isotherm.
            :type parameters: List[float]
            :return: The vector of calculated loadings [mol/kg].
            :rtype: List[float]
           )")
      .def("get_deviation", &BaseIsothermModel::GetDeviation,
           pybind11::arg("pressures"), pybind11::arg("ExperimentalLoadings"),
           pybind11::arg("temperature"), pybind11::arg("parameters"),
           pybind11::arg("DeviationEquation"), R"(
            Calculates the deviation between experimental loadings and calculated loadings using the specified deviation function.

            :param pressures: The vector of pressure values [Pa].
            :type pressures: List[float]
            :param experimental_loadings: The vector of experimental loadings [mol/kg].
            :type experimental_loadings: List[float]
            :param temperature: The temperature value [K].
            :type temperature: float
            :param parameters: The parameters for the isotherm.
            :type parameters: List[float]
            :param deviation_equation: The name of the deviation function to be used.
            :type deviation_equation: str
            :return: The calculated deviation value.
            :rtype: float
           )")
      .def("get_mixture_loading", &BaseIsothermModel::GetMixtureLoading,
           pybind11::arg("pressure"), pybind11::arg("temperature"),
           pybind11::arg("BulkComposition"), pybind11::arg("parameters"), R"(
            Calculates the partial loadings for the given pressure, temperature, bulk_composition, and parameters.

            :param pressure: The pressure value [Pa].
            :type pressure: float
            :param temperature: The temperature value [K].
            :type temperature: float
            :param bulk_composition: The bulk composition.
            :type bulk_composition: List[float]
            :param parameters: The parameters for the isotherm.
            :type parameters: List[float]
            :return: The calculated loading.
            :rtype: float
           )");
}

#endif