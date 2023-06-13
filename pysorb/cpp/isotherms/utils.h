#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <vector>
#include "deviation_functions.h"

std::function<double(std::vector<double>, std::vector<double>)> GetDeviationInvoker(std::string DeviationEquation, double NumberOfParameters);

#endif