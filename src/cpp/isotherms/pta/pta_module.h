#ifndef PTA_MODULE_H
#define PTA_MODULE_H

#include <pybind11/pybind11.h>

#include "pta.h"

void bindPotentialTheoryModels(pybind11::module &m) {
  pybind11::class_<PotentialTheoryModels, BaseIsothermModel>(
      m, "PotentialTheoryModels", pybind11::multiple_inheritance(), R"(
        Models to calculate loadings based on the Multicomponent Potential Theory of Adsorption
      )")
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
           pybind11::arg("fluids"), pybind11::arg("adsorbent"), R"(
            Initializes a model for multicomponent loading calculation with adsorbent data

            :param potential: The adsorption potential (Options: 'DRA', 'LEE', 'STEELE').
            :type potential: str
            :param equation_of_state: The equation of state to be used (Options: 'pr77', 'srk', 'pr77-peneloux', 'srk-peneloux').
            :type equation_of_state: str
            :param isotherm_type: The system's loading type (Options: excess, absolute).
            :type isotherm_type: str
            :param num_of_layers: The number of layers between the solid surface and the Gibb's Interface.
            :type num_of_layers: int
            :param fluids: The list of fluids to simulate multicomponent adsorption.
            :type fluids: List[pysorb.Fluid]
            :param adsorbent: The adsorbent to simulate multicomponent adsorption.
            :type adsorbent: pysorb.Adsorbent
            :returns: The class to calculate the loadings.
            :rtype: PotentialTheoryModels
           )")
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
           pybind11::arg("fluids"), R"(
            Initializes a model for multicomponent loading calculation without adsorbent data

            :param potential: The adsorption potential (Options: 'DRA', 'LEE', 'STEELE')
            :type potential: str
            :param equation_of_state: The equation of state to be used (Options: 'pr77', 'srk', 'pr77-peneloux', 'srk-peneloux')
            :type equation_of_state: str
            :param isotherm_type: The system's loading type (Options: excess, absolute)
            :type isotherm_type: str
            :param num_of_layers: The number of layers between the solid surface and the Gibb's Interface
            :type num_of_layers: int
            :param fluids: The list of fluids to simulate multicomponent adsorption
            :type fluids: List[pysorb.Fluid]
            :returns: The class to calculate the loadings
            :rtype: PotentialTheoryModels
           )")
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
           pybind11::arg("fluid"), pybind11::arg("adsorbent"),
           R"(
            Initializes a model for monocomponent loading calculation with adsorbent data

            :param potential: The adsorption potential (Options: 'DRA', 'LEE', 'STEELE')
            :type potential: str
            :param equation_of_state: The equation of state to be used (Options: 'pr77', 'srk', 'pr77-peneloux', 'srk-peneloux')
            :type equation_of_state: str
            :param isotherm_type: The system's loading type (Options: excess, absolute)
            :type isotherm_type: str
            :param num_of_layers: The number of layers between the solid surface and the Gibb's Interface
            :type num_of_layers: int
            :param fluids: The fluid to calculate monocomponent adsorption
            :type fluid: pysorb.Fluid
            :param adsorbent: The adsorbent to calculate monocomponent adsorption
            :type adsorbent: pysorb.Adsorbent
            :returns: The class to calculate the loadings
            :rtype: PotentialTheoryModels
           )")
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
           pybind11::arg("fluid"),
           R"(
            Initializes a model for monocomponent loading calculation without adsorbent data

            :param potential: The adsorption potential (Options: 'DRA', 'LEE', 'STEELE')
            :type potential: str
            :param equation_of_state: The equation of state to be used (Options: 'pr77', 'srk', 'pr77-peneloux', 'srk-peneloux')
            :type equation_of_state: str
            :param isotherm_type: The system's loading type (Options: excess, absolute)
            :type isotherm_type: str
            :param num_of_layers: The number of layers between the solid surface and the Gibb's Interface
            :type num_of_layers: int
            :param fluids: The fluid to calculate monocomponent adsorption
            :type fluid: pysorb.Fluid
            :returns: The class to calculate the loadings
            :rtype: PotentialTheoryModels
           )");
}

#endif