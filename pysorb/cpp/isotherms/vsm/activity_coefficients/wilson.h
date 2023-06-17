#ifndef WILSON_H
#define WILSON_H

#include <vector>
#include <cmath>

#include "../../empirical/_isotherms.h"
#include "../../empirical/extended_isotherms.h"
#include "../../optimization_algorithms/brent.h"
#include "../../optimization_algorithms/nelder_mead.h"
#include "../../data_classes.h"

double GetLoadingWILSON(double Pressure, std::vector<double> Parameters);

std::vector<double> GetMixtureLoadingWILSON(double Pressure, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters);

#endif