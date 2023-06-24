/**
 * @file deviation_functions.h
 * @brief Contains functions for calculating different types of deviations.
 */

#ifndef DEVIATION_FUNCTIONS_H
#define DEVIATION_FUNCTIONS_H

#include <vector>
#include <cmath>

/**
 * @brief Calculates the sum of squared errors (SSE) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The SSE deviation.
 */
double GetSSEDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the average relative error (ARE) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The ARE deviation.
 */
double GetAREDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the absolute error (EABS) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The EABS deviation.
 */
double GetEABSDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the hybrid deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The hybrid deviation.
 */
double GetHYBRIDDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the mean percent squared deviation (MPSD).
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The MPSD deviation.
 */
double GetMPSDDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the sum of relative errors (SORE) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The SORE deviation.
 */
double GetSOREDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the chi-squared (CHI2) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The CHI2 deviation.
 */
double GetCHI2Deviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the R-squared (RS) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The RS deviation.
 */
double GetRSDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);

/**
 * @brief Calculates the R-squared (RS) deviation.
 * @param ExperimentalLoadings The vector of experimental loadings.
 * @param CalculatedLoadings The vector of calculated loadings.
 * @param NumberOfParameters The number of parameters in the model.
 * @return The RS deviation.
 */
double GetRABSDeviation(std::vector<double> ExperimentalLoadings, std::vector<double> CalculatedLoadings, double NumberOfParameters);
#endif
