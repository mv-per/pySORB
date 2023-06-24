#ifndef WILSON_H
#define WILSON_H

#include <cmath>
#include <vector>

#include "../../data_classes.h"
#include "../../empirical/files/_isotherms.h"
#include "../../empirical/files/extended_isotherms.h"
#include "../../optimization_algorithms/brent.h"
#include "../../optimization_algorithms/nelder_mead.h"

double GetLoadingWILSON(double Pressure, std::vector<double> Parameters);

std::vector<double>
GetMixtureLoadingWILSON(double Pressure, std::vector<double> BulkComposition,
                        std::vector<std::vector<double>> Parameters);

#endif