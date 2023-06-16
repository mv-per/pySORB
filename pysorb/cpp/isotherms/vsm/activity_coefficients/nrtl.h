#ifndef NRTL_H
#define NRTL_H

#include <vector>
#include <cmath>

#include "../../empirical/_isotherms.h"
#include "../../optimization_algorithms/nelder_mead.h"

double GetLoadingNRTL(double P, double T, std::vector<double> parameters);

#endif