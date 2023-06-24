#ifndef NRTL_H
#define NRTL_H

#include <cmath>
#include <vector>

#include "../../data_classes.h"
#include "../../empirical/files/_isotherms.h"
#include "../../empirical/files/extended_isotherms.h"
#include "../../optimization_algorithms/nelder_mead.h"

double GetLoadingNRTL(double P, double T, std::vector<double> parameters);
std::vector<double>
GetMixtureLoadingNRTL(double Pressure, std::vector<double> BulkComposition,
                      std::vector<std::vector<double>> Parameters);

#endif