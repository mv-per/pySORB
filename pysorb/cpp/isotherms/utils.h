#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <vector>
#include <numeric>
#include "deviation_functions.h"

const double GAS_CONSTANT = 8.314462618;

std::function<double(std::vector<double>, std::vector<double>)> GetDeviationInvoker(std::string DeviationEquation, double NumberOfParameters);

void CheckCompositionFraction(std::vector<double> &composition);

#endif