#ifndef FLORY_HUGGINS_H
#define FLORY_HUGGINS_H

#include <vector>
#include <cmath>

#include "../../empirical/_isotherms.h"
#include "../../optimization_algorithms/brent.h"

double GetLoadingFloryHuggins(double P, double T, std::vector<double> parameters);

#endif