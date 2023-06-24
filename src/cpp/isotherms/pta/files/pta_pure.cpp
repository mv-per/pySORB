#include "pta_pure.h"

// PurePTA::PurePTA(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers)
// {
//     this->Potential = potential;
//     this->NumberOfLayers = num_of_layers;
//     this->EquationOfState = equation_of_state;
//     this->IsothermType = isotherm_type;
//     this->Loader = GetLoadingFunction();
// };

// /**
//  * Get the calculated loading in a specific pressure
//  *
//  * @param P Pressure of the Fluid
//  * @param T Temperature of the Fluid
//  * @param potential_params Params of the Adsorption Potential for this Fluid
//  * @param Fluid Fluid properties.
//  * @return Calculated adsorbed loading
//  */
// double PurePTA::GetLoading(double P, double T, std::vector<double> potential_params, Fluid Fluid)
// {
//     this->Fluid = Fluid;
//     call_potential PotentialInvoker = GetAdsorptionPotentialInvoker(potential_params, Fluid, this->adsorbent);
//     call_mono_eos eos_caller = GetEquationOfStateInvoker(Fluid);
//     double Loading = this->Loader(P, T, potential_params, eos_caller, PotentialInvoker);
//     return Loading;
// }

// /**
//  * Get the calculated loading in a specific pressure
//  *
//  * @param P Pressure of the Fluid
//  * @param T Temperature of the Fluid
//  * @param potential_params Params of the Adsorption Potential for this Fluid
//  * @param Fluid Fluid properties.
//  * @return Calculated adsorbed loading
//  */
// double PurePTA::GetPressure(double n, double T, std::vector<double> potential_params, Fluid Fluid, double P_estimate_ = 1e6)
// {
//     this->Fluid = Fluid;
//     call_potential PotentialInvoker = GetAdsorptionPotentialInvoker(potential_params, Fluid, this->adsorbent);
//     call_mono_eos eos_caller = GetEquationOfStateInvoker(Fluid);

//     auto equilibrium = [&](double x)
//     {
//         double diff = (this->Loader(x, T, potential_params, eos_caller, PotentialInvoker) - n);
//         return diff * 1000;
//     };

//     return brent_zeroin(equilibrium, 1e6, 1e-5);
// }

// /**
//  * Get the calculated loading for different pressures
//  *
//  * @param P List of pressure points to get loading
//  * @param T Temperature of the Fluid
//  * @param potential_params Params of the Adsorption Potential for this Fluid
//  * @param Fluid Fluid properties.
//  * @return List of calculated loadings
//  */
// std::vector<double> PurePTA::GetLoadings(std::vector<double> P, double T, std::vector<double> potential_params, Fluid Fluid)
// {
//     std::vector<double>
//         loadings(P.size());

//     this->Fluid = Fluid;
//     call_potential PotentialInvoker = GetAdsorptionPotentialInvoker(potential_params, Fluid, this->adsorbent);
//     call_mono_eos eos_caller = GetEquationOfStateInvoker(Fluid);

//     for (std::size_t i = 0; i < P.size(); i++)
//     {
//         loadings[i] = this->Loader(P[i], T, potential_params, eos_caller, PotentialInvoker);
//     }

//     return loadings;
// }

// /**
//  * Get the deviation between calculated experimental loadings
//  *
//  * @param P List of pressure points to get loading
//  * @param T Temperature of the Fluid
//  * @param potential_params Params of the Adsorption Potential for this Fluid
//  * @param Fluid Fluid properties.
//  * @return List of calculated loadings
//  */
// double PurePTA::GetDeviationRange(std::string deviation_type,
//                                   std::vector<double> loading_exp, std::vector<double> P, double T, std::vector<double> potential_params, Fluid Fluid)
// {
//     deviation_caller deviation = get_deviation_function(deviation_type);

//     double difference = 0.;

//     std::vector<double> calc_loading = this->GetLoadings(P, T, potential_params, Fluid);
//     for (std::size_t i = 0; i < P.size(); i++)
//     {
//         difference += deviation(loading_exp[i], calc_loading[i]);
//     }

//     if (deviation_type.find("relative") != std::string::npos)
//     {
//         return 100. / loading_exp.size() * difference;
//     }

//     return difference;
// }

std::function<mono_eos(double, double)> GetPureEquationOfStateInvoker(std::string equation_of_state, Fluid fluid_)
{

    if (equation_of_state == "pr77")
    {
        return [=](double P, double T)
        {
            return pr77().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else if (equation_of_state == "pr77-peneloux")
    {
        fluid_ = CheckForFluidCriticalCompressibility(fluid_);
        return [=](double P, double T)
        {
            return pr77_peneloux().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else if (equation_of_state == "srk")
    {
        return [=](double P, double T)
        {
            return srk().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else if (equation_of_state == "srk-peneloux")
    {
        fluid_ = CheckForFluidCriticalCompressibility(fluid_);
        return [=](double P, double T)
        {
            return srk_peneloux().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else
    {
        throw std::invalid_argument("Equation of State not found/defined.");
    }
}

std::function<double(double, std::vector<double>)> GetPureAdsorptionPotentialInvoker(std::string potential, Fluid Fluid, Adsorbent adsorbent)
{
    if (potential == DRA_POTENTIAL)
    {
        return [=](double z, std::vector<double> potential_params)
        {
            return DRA(z, potential_params[0], potential_params[1], potential_params[2]);
        };
    }
    else if (potential == STEELE_POTENTIAL)
    {
        return [=](double z, std::vector<double> potential_params)
        {
            return STEELE(z, potential_params[0], Fluid.LennardJonnesDiameter, adsorbent);
        };
    }
    else if (potential == LEE_POTENTIAL)
    {
        return [=](double z, std::vector<double> potential_params)
        {
            return LEE(z, potential_params[0], Fluid.LennardJonnesDiameter, adsorbent);
        };
    }
    else
    {
        throw std::invalid_argument("Isotherm potential not found/defined.");
    }
}

double GetDRAPureLoading(
    double BulkPressure,
    double Temperature,
    std::vector<double> Parameters,
    call_mono_eos EquationOfStateInvoker,
    call_potential PotentialInvoker,
    Fluid Fluid,
    std::size_t NumberOfLayers,
    std::string IsothermType,
    std::string Potential)
{

    double Pz, f_eps, delta, eps;

    struct mono_eos bulk = EquationOfStateInvoker(BulkPressure, Temperature);
    PTASolver solver = PTASolver(bulk.fug, Temperature, EquationOfStateInvoker);

    // Initial estimates for the adsorbed region
    Pz = BulkPressure;

    delta = Parameters[1] / (double)NumberOfLayers;

    std::vector<double>
        z = linespace(Parameters[1], delta, NumberOfLayers);

    std::vector<double> integral(NumberOfLayers, 0.0);

    for (std::size_t i = 0; i < z.size(); i++)
    {
        eps = PotentialInvoker(z[i], Parameters);
        f_eps = GetAdsorptionPotentialEnergy(eps, Temperature, Potential);
        Pz = solver.findOptimizedPressure(Pz, f_eps);

        struct mono_eos adsorbed = EquationOfStateInvoker(Pz, Temperature);
        integral[i] = GetCalculatedAdsorbedDensity(adsorbed.dens, bulk.dens, IsothermType) * 1e-3;
    }
    return integrate(integral, delta);
}

double GetLJPureLoading(
    double BulkPressure,
    double Temperature,
    std::vector<double> Parameters,
    call_mono_eos EquationOfStateInvoker,
    call_potential PotentialInvoker,
    Fluid Fluid,
    std::size_t NumberOfLayers,
    std::string IsothermType,
    std::string Potential)
{
    double Pz, f_eps, delta, eps;

    struct mono_eos bulk = EquationOfStateInvoker(BulkPressure, Temperature);
    PTASolver solver = PTASolver(bulk.fug, Temperature, EquationOfStateInvoker);

    // Initial estimates for the adsorbed region
    Pz = BulkPressure;
    if (!Fluid.LennardJonnesDiameter)
    {
        throw std::invalid_argument("Missing the Lennard-Jonnes diameter for the Fluid " + Fluid.Name);
    }
    double minimal_space = 0.7 * Fluid.LennardJonnesDiameter;

    std::vector<double>
        z = linespace(Parameters[1] / 2.0, minimal_space, NumberOfLayers);

    delta = z[1] - z[0];

    std::vector<double> integral(NumberOfLayers, 0.0);

    for (std::size_t i = 0; i < z.size(); i++)
    {
        eps = PotentialInvoker(z[i], Parameters) + PotentialInvoker(Parameters[1] - z[i], Parameters);
        f_eps = GetAdsorptionPotentialEnergy(eps, Temperature, Potential);
        Pz = solver.findOptimizedPressure(Pz, f_eps);

        struct mono_eos adsorbed = EquationOfStateInvoker(Pz, Temperature);
        integral[i] = GetCalculatedAdsorbedDensity(adsorbed.dens, bulk.dens, IsothermType);
    }
    return -Parameters[2] * integrate(integral, delta) * 1e-7;
}
