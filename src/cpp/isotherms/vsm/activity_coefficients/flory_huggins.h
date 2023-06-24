#ifndef FLORY_HUGGINS_H
#define FLORY_HUGGINS_H

#include <cmath>
#include <vector>

#include "../../data_classes.h"
#include "../../empirical/files/_isotherms.h"
#include "../../empirical/files/extended_isotherms.h"
#include "../../optimization_algorithms/brent.h"
#include "../../optimization_algorithms/nelder_mead.h"

double GetLoadingFloryHuggins(double P, double T,
                              std::vector<double> parameters);

std::vector<double>
GetMixtureLoadingFloryHuggins(double Pressure,
                              std::vector<double> BulkComposition,
                              std::vector<std::vector<double>> Parameters);

#endif