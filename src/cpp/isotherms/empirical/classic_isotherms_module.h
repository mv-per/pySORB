#ifndef CLASSIC_ISOTHERMS_MODULE_H
#define CLASSIC_ISOTHERMS_MODULE_H

#include <pybind11/pybind11.h>

#include "classic_isotherms.h"

void bindEmpiricalIsotherms(pybind11::module &m) {

  pybind11::class_<EmpiricalIsotherms, BaseIsothermModel>(
      m, "EmpiricalIsotherms", pybind11::multiple_inheritance())
      .def(pybind11::init<std::string>(), pybind11::arg("model"));
}

#endif