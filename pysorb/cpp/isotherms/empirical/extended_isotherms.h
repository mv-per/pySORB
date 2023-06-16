#ifndef EXTENDED_ISOTHERMS_H
#define EXTENDED_ISOTHERMS_H

#include <vector>
#include <cassert>
#include "../utils.h"

std::vector<double> extended_langmuir(double Pressure, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters);

#endif