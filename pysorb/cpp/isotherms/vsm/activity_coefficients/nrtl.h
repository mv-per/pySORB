#ifndef NRTL_H
#define NRTL_H

#include <vector>
#include <cmath>

#include "../../empirical/_isotherms.h"
#include "../../empirical/extended_isotherms.h"
#include "../../optimization_algorithms/nelder_mead.h"
#include "../../data_classes.h"

double GetLoadingNRTL(double P, double T, std::vector<double> parameters);
std::vector<double> GetMixtureLoadingNRTL(double Pressure, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters);

#endif