
/**
 * @file _isotherms.h
 * @brief Contains implementations of various isotherm functions.
 */

#ifndef _ISOTHERMS_H
#define _ISOTHERMS_H
#include "../utils.h"
#include <cmath>
#include <vector>

/**
 * @brief Calculates the loading using the Freundlich isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [K, n]
 * @return The calculated loading.
 */
double freundlich(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the Freundlich isotherm equation for gas adsorption.
 *
 * This function calculates the Freundlich isotherm equation for gas adsorption,
 * which describes the relationship between the amount of gas adsorbed by a
 * solid and the gas pressure and temperature. The equation is given by:
 *
 *    K * Pressure^inverse_n
 *
 * Where:
 * - `Pressure` is the pressure of the gas.
 * - `Temperature` is the temperature of the system.
 * - `Parameters` is a vector containing three parameters:
 *   - `Parameters[0]` represents the constant K in the equation.
 *   - `Parameters[1]` represents the parameter used in the exponential term.
 *   - `Parameters[2]` represents a constant used in the calculation of
 * inverse_n.
 *
 * @param Pressure The pressure of the gas.
 * @param Temperature The temperature of the system.
 * @param Parameters A vector containing three parameters: K, exponential term
 * parameter, and constant.
 *
 * @return The calculated value of the Freundlich equation.
 */
double freundlich_2(double Pressure, double Temperature,
                    std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Langmuir isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, b]
 * @return The calculated loading.
 */
double langmuir(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Langmuir isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, b_inf, Q]
 * @return The calculated loading.
 */
double langmuir_2(double Pressure, double Temperature,
                  std::vector<double> Parameters);

/**
 * @brief Calculates the dual Langmuir equation for gas adsorption.
 *
 * This function calculates the dual Langmuir equation for gas adsorption,
 * which models the adsorption of a gas on a surface with two types of active
 * sites. The equation is given by:
 *
 *    (Parameters[0] * Parameters[1] * Pressure) / (1 + Parameters[1] *
 * Pressure)
 *  + (Parameters[2] * Parameters[3] * Pressure) / (1 + Parameters[3] *
 * Pressure)
 *
 * Where:
 * - `Pressure` is the pressure of the gas.
 * - `Parameters` is a vector containing four parameters:
 *   - `Parameters[0]` represents the constant for the first active site.
 *   - `Parameters[1]` represents the Langmuir isotherm parameter for the first
 * active site.
 *   - `Parameters[2]` represents the constant for the second active site.
 *   - `Parameters[3]` represents the Langmuir isotherm parameter for the second
 * active site.
 *
 * @param Pressure The pressure of the gas.
 * @param Parameters A vector containing four parameters: constant1, isotherm
 * parameter1, constant2, and isotherm parameter2.
 *
 * @return The calculated value of the dual Langmuir equation.
 */
double dual_langmuir(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Redlich-Peterson isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [K_RP, a_RP, g_RP]
 * @return The calculated loading.
 */
double redlich_peterson(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Sips isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, b, n]
 * @return The calculated loading.
 */
double sips(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the temperature dependent Sips isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, chi, t_0, alpha,
 * b_0, Q]
 * @return The calculated loading.
 */
double sips_2(double Pressure, double Temperature,
              std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Toth isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_max, b, t]
 * @return The calculated loading.
 */
double toth(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the temperature dependent Toth isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, chi, t_0, alpha,
 * b_0, Q]
 * @return The calculated loading.
 */
double toth_2(double Pressure, double Temperature,
              std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Unilan isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, b_med, s]
 * @return The calculated loading.
 */
double unilan(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Keller-Staudt-Toth isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [N_max, b, alpha_m, beta]
 * @return The calculated loading.
 */
double keller_staudt_toth(double Pressure, std::vector<double> Parameters);

/**
 * @brief Calculates the loading using the Jensen-Seaton isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [K , a, kappa, c]
 * @return The calculated loading.
 */
double jensen_seaton(double Pressure, std::vector<double> Parameters);

#endif