#include "pta_mixture.h"

double get_equilibrium_difference(double Pz,
								  std::vector<double> adsorbed_composition,
								  std::vector<double> bulk_fugacity,
								  std::vector<double> f_eps,
								  double T,
								  call_mix_eos eos)
{

	if (std::isinf(Pz) == 1 || std::isnan(Pz) == 1 || Pz <= 0)
	{
		return 1.0;
	}

	double difference = 0.0;
	struct mix_eos adsorbed = eos(adsorbed_composition, Pz, T);

	// Find the Equilibrium difference
	for (size_t i = 0; i < adsorbed_composition.size(); i++)
	{
		difference += bulk_fugacity[i] * f_eps[i] - adsorbed.fug[i];
	}
	return difference;
}

struct Equilibrium find_equilibrium_properties(double Pz_old, std::vector<double> xold, std::vector<double> fb,
											   std::vector<double> f_eps, double T, call_mix_eos eos)
{
	double xsum, composition_balancer;
	double Pz = Pz_old;
	std::size_t i, k;
	std::vector<double> new_adsorbed_composition = xold;
	struct mix_eos adsorbed = eos(new_adsorbed_composition, Pz, T);
	new_adsorbed_composition = get_new_adsorbed_fractions(fb, f_eps, adsorbed.phi, Pz);
	k = 0;

	do
	{

		// Set the Equilibrium lambda function
		auto min = [&](double x)
		{
			return get_equilibrium_difference(x, new_adsorbed_composition, fb, f_eps, T, eos);
		};

		// Find the Equilibrium pressure for the current composition
		Pz = brent_zeroin(min, Pz, 1.0e-8);

		// Get the fluid properties
		struct mix_eos adsorbed = eos(new_adsorbed_composition, Pz, T);

		// Finds the difference from the current and previous composition
		xsum = 0;
		for (i = 0; i < fb.size(); i++)
		{
			xsum += std::fabs(new_adsorbed_composition[i] / xold[i]) - 1.0;
		}

		xsum *= xsum;

		// Limits to 25 iterations
		if (k >= 25)
		{
			break;
		}

		// Update composition
		xold = new_adsorbed_composition;
		k++;
		new_adsorbed_composition = get_new_adsorbed_fractions(fb, f_eps, adsorbed.phi, Pz);

	} while (xsum >= 1.0e-12);
	// Keep updating the composition until it doesn't changes anymore

	struct Equilibrium data = {Pz, new_adsorbed_composition};
	return data;
}

std::function<mix_eos(std::vector<double>, double, double)> GetMixtureEquationOfStateInvoker(std::string equation_of_state, std::vector<Fluid> fluids)
{
	if (equation_of_state == "pr77")
	{
		return [=](std::vector<double> composition, double P, double T)
		{
			return pr77().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else if (equation_of_state == "pr77-peneloux")
	{
		for (std::size_t i = 0; i < fluids.size(); i++)
		{
			fluids[i] = CheckForFluidCriticalCompressibility(fluids[i]);
		}
		return [=](std::vector<double> composition, double P, double T)
		{
			return pr77_peneloux().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else if (equation_of_state == "srk")
	{
		return [=](std::vector<double> composition, double P, double T)
		{
			return srk().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else if (equation_of_state == "srk-peneloux")
	{
		for (std::size_t i = 0; i < fluids.size(); i++)
		{
			fluids[i] = CheckForFluidCriticalCompressibility(fluids[i]);
		}
		return [=](std::vector<double> composition, double P, double T)
		{
			return srk_peneloux().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else
	{
		throw std::invalid_argument("Equation of State not found/defined.");
	}
}

call_potential_mix GetMixtureAdsorptionPotentialInvoker(std::string potential, Adsorbent adsorbent)
{
	if (potential == DRA_POTENTIAL)
	{
		return [=](double z, Fluid fluid, std::vector<double> parameters)
		{
			return DRA(z, parameters[0], parameters[1], parameters[2]);
		};
	}
	else if (potential == STEELE_POTENTIAL)
	{
		return [=](double z, Fluid fluid, std::vector<double> parameters)
		{
			return STEELE(z, parameters[0], fluid.LennardJonnesDiameter, adsorbent);
		};
	}
	else if (potential == LEE_POTENTIAL)
	{
		return [=](double z, Fluid fluid, std::vector<double> parameters)
		{
			return LEE(z, parameters[0], fluid.LennardJonnesDiameter, adsorbent);
		};
	}
	else
	{
		throw std::invalid_argument("Isotherm potential not found/defined.");
	}
}

/**
 * Get the potential function based on its energy interation and temperature
 *
 */
std::vector<double> GetDRAMixtureLoading(double BulkPressure,
										 double Temperature,
										 std::vector<double> BulkComposition,
										 std::vector<std::vector<double>> Parameters,
										 std::vector<Fluid> fluids,
										 std::size_t NumberOfLayers,
										 std::string IsothermType,
										 std::string Potential,
										 std::string EquationOfState,
										 Adsorbent adsorbent)
{
	std::size_t i, j;
	std::size_t ncomp = BulkComposition.size();
	double eps;

	call_potential_mix get_potential = GetMixtureAdsorptionPotentialInvoker(Potential, adsorbent);

	call_mix_eos eos = GetMixtureEquationOfStateInvoker(EquationOfState, fluids);

	struct mix_eos bulk = eos(BulkComposition, BulkPressure, Temperature);

	// Initial Conditions
	double Pz = BulkPressure;
	std::vector<double> adsorbed_composition = BulkComposition;

	std::vector<std::vector<double>> z(ncomp, std::vector<double>(NumberOfLayers, 0.0));
	std::vector<double> deltas(ncomp, 0.0);

	for (j = 0; j < ncomp; j++)
	{
		deltas[j] = Parameters[j][1] / (double)NumberOfLayers;
		z[j] = linespace(Parameters[j][1], deltas[j], NumberOfLayers);
	}
	// Instantiate blank vectors
	std::vector<double> f_eps(ncomp, 0.0);
	std::vector<double> adsorbed_loading(ncomp, 0.0);
	std::vector<std::vector<double>> mix_integral(ncomp, std::vector<double>(NumberOfLayers, 0.0));

	for (i = 0; i < NumberOfLayers; i++)
	{
		for (j = 0; j < ncomp; j++)
		{
			eps = get_potential(z[j][i], fluids[j], Parameters[j]);
			f_eps[j] = GetAdsorptionPotentialEnergy(eps, Temperature, Potential);
		}

		// Find pressure and composition of equilibrium
		struct Equilibrium Props = find_equilibrium_properties(Pz, adsorbed_composition, bulk.fug, f_eps, Temperature, eos);
		Pz = Props.Pz;
		adsorbed_composition = Props.x;

		struct mix_eos adsorbed = eos(adsorbed_composition, Pz, Temperature);

		for (j = 0; j < ncomp; j++)
		{
			mix_integral[j][i] = GetCalculatedAdsorbedDensity(adsorbed.dens * adsorbed_composition[j], bulk.dens * BulkComposition[j], IsothermType);
		}
	}

	for (j = 0; j < ncomp; j++)
	{
		adsorbed_loading[j] = integrate(mix_integral[j], deltas[j]) / 1000.0;
	}

	return adsorbed_loading;
}
std::vector<double> GetLJMixtureLoading(double BulkPressure,
										double Temperature,
										std::vector<double> BulkComposition,
										std::vector<std::vector<double>> Parameters,
										std::vector<Fluid> fluids,
										std::size_t NumberOfLayers,
										std::string IsothermType,
										std::string Potential,
										std::string EquationOfState,
										Adsorbent adsorbent)
{
	std::size_t i;
	std::size_t j;
	double minimal_space, eps;
	std::size_t ncomp = BulkComposition.size();

	call_potential_mix get_potential = GetMixtureAdsorptionPotentialInvoker(Potential, adsorbent);

	call_mix_eos eos = GetMixtureEquationOfStateInvoker(EquationOfState, fluids);

	struct mix_eos bulk = eos(BulkComposition, BulkPressure, Temperature);

	// Initial Conditions
	double Pz = BulkPressure;
	std::vector<double> adsorbed_composition = BulkComposition;

	std::vector<std::vector<double>> z(ncomp, std::vector<double>(NumberOfLayers, 0.0));
	std::vector<double> deltas(ncomp, 0.0);
	std::vector<double> f_eps(ncomp, 0.0);
	std::vector<double> adsorbed_loading(ncomp, 0.0);
	std::vector<std::vector<double>> mix_integral(ncomp, std::vector<double>(NumberOfLayers, 0.0));

	for (j = 0; j < ncomp; j++)
	{
		if (!fluids[j].LennardJonnesDiameter)
		{
			throw std::invalid_argument("Missing the Lennard-Jonnes diameter for the fluid " + fluids[j].Name);
		}

		minimal_space = 0.7 * fluids[j].LennardJonnesDiameter;

		z[j] = linespace(Parameters[j][1] / 2.0, minimal_space, NumberOfLayers);
		deltas[j] = z[j][1] - z[j][0];
	}

	for (i = 0; i < NumberOfLayers; i++)
	{
		for (j = 0; j < ncomp; j++)
		{
			eps = get_potential(z[j][i], fluids[j], Parameters[j]) + get_potential(Parameters[j][1] - z[j][i], fluids[j], Parameters[j]);
			f_eps[j] = GetAdsorptionPotentialEnergy(eps, Temperature, Potential);

			// TODO: study the proper rate with parameter L and the adsorbent sizes
		}

		// Find pressure and composition of equilibrium
		struct Equilibrium Props = find_equilibrium_properties(Pz, adsorbed_composition, bulk.fug, f_eps, Temperature, eos);
		Pz = Props.Pz;
		adsorbed_composition = Props.x;

		struct mix_eos adsorbed = eos(adsorbed_composition, Pz, Temperature);

		for (j = 0; j < ncomp; j++)
		{
			mix_integral[j][i] = GetCalculatedAdsorbedDensity(adsorbed.dens * adsorbed_composition[j], bulk.dens * BulkComposition[j], IsothermType);
		}
	}

	for (j = 0; j < ncomp; j++)
	{
		adsorbed_loading[j] = -Parameters[j][2] * integrate(mix_integral[j], deltas[j]) * 1.0e-7;
	}

	return adsorbed_loading;
}