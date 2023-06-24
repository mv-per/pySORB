#ifndef PTA_MIXTURE_H
#define PTA_MIXTURE_H

#include "pta_helper.h"
#include "helpers.h"
#include "../equations_of_state/eos.h"
#include "../optimization_algorithms/brent.h"
#include "adsorption_potentials.h"
#include "pta_solver.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <cfloat>

std::function<mix_eos(std::vector<double>, double, double)> GetMixtureEquationOfStateInvoker(std::string equation_of_state, Fluid fluid_);
std::function<double(double, Fluid, std::vector<double>)> GetMixtureAdsorptionPotentialInvoker(std::string potential, Adsorbent adsorbent);
/**
 * Get the potential function based on its energy interation and temperature
 *
 */
std::vector<double> GetDRAMixtureLoading(double BulkPressure,
										 double Temperature,
										 std::vector<double> BulkComposition,
										 std::vector<std::vector<double>> Parameters,
										 std::vector<Fluid> fluids,
										 std::size_t NumberOfLayers,
										 std::string IsothermType,
										 std::string Potential,
										 std::string EquationOfState,
										 Adsorbent adsorbent);
std::vector<double> GetLJMixtureLoading(double BulkPressure,
										double Temperature,
										std::vector<double> BulkComposition,
										std::vector<std::vector<double>> Parameters,
										std::vector<Fluid> fluids,
										std::size_t NumberOfLayers,
										std::string IsothermType,
										std::string Potential,
										std::string EquationOfState,
										Adsorbent adsorbent);

#endif