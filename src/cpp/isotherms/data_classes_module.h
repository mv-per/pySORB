#ifndef DATA_CLASSES_MODULE_H
#define DATA_CLASSES_MODULE_H

#include <pybind11/pybind11.h>

#include "data_classes.h"

void bindDataClasses(pybind11::module &m) {
  pybind11::class_<Adsorbent>(m, "Adsorbent", R"(The system's adsorbent)")
      .def_readwrite("name", &Adsorbent::Name, R"(Name of the adsorbent)")
      .def_readwrite("solid_diameter", &Adsorbent::SolidDiameter,
                     R"(
                     The adsorbent Lennard-Jonnes surface diameter  [Angstrom]
                     )")
      .def_readwrite("atomic_density", &Adsorbent::SolidAtomicDensity,
                     R"(
          The Atomic density of the adsorbent surface [Atoms/sq.Angstrom]
          )")
      .def(pybind11::init<>(), R"(Initialize a base adsobent)")
      .def(pybind11::init([](std::string name, double solid_diameter,
                             double solid_atomic_density) {
             return new Adsorbent(name, solid_diameter, solid_atomic_density);
           }),
           pybind11::arg("name"), pybind11::arg("diameter"),
           pybind11::arg("atomic_density"),
           R"pobdoc(
            Returns the adsorbent with the given name.

            :param name: The name of the adsorbent.
            :type name: str
            :param diameter: Diameters of the adsorbent wall [Angstrom].
            :type diameter: float
            :param atomic_density: Atomic density of the adsorbent wall [unit].
            :type atomic_density: float
            :return: The Adsorbent.
            :rtype: Adsorbent
        )pobdoc");

  pybind11::class_<Fluid>(m, "Fluid", R"(
                     The adsorbate fluid
                     )")
      .def_readwrite("name", &Fluid::Name, R"(Name of the fluid)")
      .def_readwrite("critical_pressure", &Fluid::CriticalPressure,
                     R"(
                     Critical Pressure of the fluid [Pascal]
                     )")
      .def_readwrite("critical_temperature", &Fluid::CriticalTemperature,
                     R"(
          Critical Temperature of the fluid [Kelvin]
          )")
      .def_readwrite("accentric_factor", &Fluid::AccentricFactor,
                     R"(
                     Accentric Factor of the fluid [-]
                     )")
      .def_readwrite("critical_compressibility",
                     &Fluid::CriticalCompressibility,
                     R"(
                     Critical compressibility of the fluid [-]
                     )")
      .def_readwrite("lennard_jonnes_diameter", &Fluid::LennardJonnesDiameter,
                     R"(Lennard Jonnes diameter of the fluid [Angstrom])")
      .def(pybind11::init([](std::string name, double critical_pressure,
                             double critical_temperature,
                             double accentric_factor,
                             double critical_compressibility,
                             double lennard_jonnes_diameter) {
             return new Fluid(name, critical_pressure, critical_temperature,
                              accentric_factor, critical_compressibility,
                              lennard_jonnes_diameter);
           }),
           pybind11::arg("name"), pybind11::arg("critical_pressure"),
           pybind11::arg("critical_temperature"),
           pybind11::arg("accentric_factor"),
           pybind11::arg("critical_compressibility"),
           pybind11::arg("lennard_jonnes_diameter"), R"(
              Initialize a fluid with Lennard-Jonnes and Critical compressibility.

              :param name: The name of the fluid.
              :type name: str
              :param critical_pressure: The fluid critical pressure [Pa].
              :type critical_pressure: float
              :param critical_temperature: The fluid critical temperature [K].
              :type critical_temperature: float
              :param accentric_factor: The fluid accentric factor [-].
              :type accentric_factor: float
              :param critical_compressibility: The fluid critical compressibility [-].
              :type critical_compressibility: float
              :param lennard_jonnes_diameter: The fluid Lennard-Jonnes diameter [Angstrom].
              :type lennard_jonnes_diameter: float
              :return: The Fluid.
              :rtype: Fluid
           )")
      .def(pybind11::init([](std::string name, double critical_pressure,
                             double critical_temperature,
                             double accentric_factor) {
             return new Fluid(name, critical_pressure, critical_temperature,
                              accentric_factor);
           }),
           pybind11::arg("name"), pybind11::arg("critical_pressure"),
           pybind11::arg("critical_temperature"),
           pybind11::arg("accentric_factor"), R"(
              Initialize a fluid.

              :param name: The name of the fluid.
              :type name: str
              :param critical_pressure: The fluid critical pressure [Pa].
              :type critical_pressure: float
              :param critical_temperature: The fluid critical temperature [K].
              :type critical_temperature: float
              :param accentric_factor: The fluid accentric factor [-].
              :type accentric_factor: float
              :return: The Fluid.
              :rtype: Fluid
           )")
      .def(pybind11::init<>(), R"(
        Initialize a base fluid
      )");
}

#endif