#ifndef PTA_SOLVER_H
#define PTA_SOLVER_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include "../equations_of_state/eos.h"
#include "../optimization_algorithms/brent.h"
#include "../optimization_algorithms/fmin.h"
#include "../data_classes.h"

class PTASolver
{
private:
    double BulkFugacity, PotentialEnergy, Temperature;
    call_mono_eos EquationOfState;
    call_mix_eos MixEquationOfState;
    double OptimizedPressure;

    /// Default tolerance
    double DEFAULT_TOL = 1.48e-6;

    double equilibrium(double p_,
                       double f_eps);

public:
    PTASolver(double bulk_fugacity, double temperature, call_mono_eos eos);
    PTASolver(double bulk_fugacity, double temperature, call_mix_eos eos);
    double findOptimizedPressure(double InitialEstimate, double f_eps);
};
#endif