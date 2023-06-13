
/**
 * @file _isotherms.h
 * @brief Contains implementations of various isotherm functions.
 */

#ifndef _ISOTHERMS_H
#define _ISOTHERMS_H
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
 * @brief Calculates the loading using the Langmuir isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_i_max, b]
 * @return The calculated loading.
 */
double langmuir(double Pressure, std::vector<double> Parameters);

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
 * @brief Calculates the loading using the Toth isotherm.
 * @param Pressure The pressure value.
 * @param Parameters The parameters for the isotherm. [n_max, b, t]
 * @return The calculated loading.
 */
double toth(double Pressure, std::vector<double> Parameters);

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