#ifndef BRENT_H
#define BRENT_H

#include <functional>
#include <cmath>

/**
 * Returns the first value in the same unit as de second one
 *
 * @param value1 Value to change sign.
 * @param value2 Value to find sign.
 * @return value1 in the same unit as value2
 */
double sign(double value1, double value2);
/**
 * Perform the Richard Brent's improvements to Dekker's zeroin algorithm
 *
 * @param fun Pointer to Function that receives a double and returns a double (function to be minimized).
 * @param x Initial estimate.
 * @param tol Minimum `TOL` .
 * @return Value of the value with minimal :fun: value.
 */
double brent_zeroin(std::function<double(double)> fun, double x, double tol);
double brent_zeroin2(std::function<double(double)> fun, double x, double tol);

#endif