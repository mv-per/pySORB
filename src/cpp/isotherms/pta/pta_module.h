#ifndef PTA_MODULE_H
#define PTA_MODULE_H

#include <pybind11/pybind11.h>

#include "pta.h"

void bindPotentialTheoryModels(pybind11::module &m) {
  pybind11::class_<PotentialTheoryModels, BaseIsothermModel>(
      m, "PotentialTheoryModels", pybind11::multiple_inheritance())
      .def(pybind11::init(
               [](std::string potential, std::string equation_of_state,
                  std::string isotherm_type, std::size_t num_of_layers,
                  std::vector<Fluid> fluids, Adsorbent adsorbent) {
                 return new PotentialTheoryModels(potential, equation_of_state,
                                                  isotherm_type, num_of_layers,
                                                  fluids, adsorbent);
               }),
           pybind11::arg("potential"), pybind11::arg("equation_of_state"),
           pybind11::arg("isotherm_type"), pybind11::arg("num_of_layers"),
           pybind11::arg("fluids"), pybind11::arg("adsorbent"))
      .def(pybind11::init(
               [](std::string potential, std::string equation_of_state,
                  std::string isotherm_type, std::size_t num_of_layers,
                  std::vector<Fluid> fluids) {
                 return new PotentialTheoryModels(potential, equation_of_state,
                                                  isotherm_type, num_of_layers,
                                                  fluids);
               }),
           pybind11::arg("potential"), pybind11::arg("equation_of_state"),
           pybind11::arg("isotherm_type"), pybind11::arg("num_of_layers"),
           pybind11::arg("fluids"))
      .def(pybind11::init(
               [](std::string potential, std::string equation_of_state,
                  std::string isotherm_type, std::size_t num_of_layers,
                  Fluid fluid, Adsorbent adsorbent) {
                 return new PotentialTheoryModels(potential, equation_of_state,
                                                  isotherm_type, num_of_layers,
                                                  fluid, adsorbent);
               }),
           pybind11::arg("potential"), pybind11::arg("equation_of_state"),
           pybind11::arg("isotherm_type"), pybind11::arg("num_of_layers"),
           pybind11::arg("fluid"), pybind11::arg("adsorbent"))
      .def(pybind11::init([](std::string potential,
                             std::string equation_of_state,
                             std::string isotherm_type,
                             std::size_t num_of_layers, Fluid fluid) {
             return new PotentialTheoryModels(potential, equation_of_state,
                                              isotherm_type, num_of_layers,
                                              fluid);
           }),
           pybind11::arg("potential"), pybind11::arg("equation_of_state"),
           pybind11::arg("isotherm_type"), pybind11::arg("num_of_layers"),
           pybind11::arg("fluid"));
}

#endif