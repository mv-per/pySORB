#ifndef VSM_MODULE_H
#define VSM_MODULE_H

#include <pybind11/pybind11.h>

#include "vsm.h"

void bindVacancySolutionMethod(pybind11::module &m) {
  pybind11::class_<VacancySolutionMethod, BaseIsothermModel>(
      m, "VacancySolutionMethod", pybind11::multiple_inheritance())
      .def(pybind11::init<std::string>(), pybind11::arg("model"), R"(
        Initializes as VacancySolutionMethod

        :param model: The Activity coefficient method, Options=('wilson', 'nrtl', 'flory-huggins').
        :type model: str
      )");
}

#endif