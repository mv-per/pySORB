#ifndef WILSON_H
#define WILSON_H

#include <vector>
#include <cmath>

#include "../../empirical/_isotherms.h"
#include "../../optimization_algorithms/brent.h"

double GetLoadingWILSON(double P, double T, std::vector<double> parameters);

#endif