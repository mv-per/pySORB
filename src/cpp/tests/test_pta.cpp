#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <gtest/gtest.h>
#include <vector>
#include <cassert>
#include "../isotherms/pta/pta.h"

TEST(test_pta_pure, TestSinglePTALoading)
{

    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    PotentialTheoryModels pta_model = PotentialTheoryModels(potential, "pr77", "excess", 100, co2);

    std::vector<double> DRA_params = {7880.19, 0.29, 2.};

    double adsorbed_loading = pta_model.GetPureLoading(1.e6, 303., DRA_params);

    EXPECT_NEAR(adsorbed_loading, 5.8862, 1e-4);
};

TEST(test_pta_pure, TestSinglePTALoading2)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    auto pta_model = PotentialTheoryModels(potential, "pr77", "excess", 100, co2);

    std::vector<double> DRA_params = {7880.19, 0.29, 2.};

    double adsorbed_loading = pta_model.GetPureLoading(10.e6, 323., DRA_params);

    EXPECT_NEAR(adsorbed_loading, 5.62658, 1e-4);
};

TEST(test_pta_pure, TestGetMultiplePurePTA)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;
    PotentialTheoryModels pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, co2);

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> ExpectedLoadings = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    EXPECT_NEAR(0.0001, 0.0001, 0.2);
    std::vector<double>
        CalculatedLoadings = pta_model.GetPureLoadings(pressures, 318.2, CO2DRAParams);

    ASSERT_EQ(pressures.size(), ExpectedLoadings.size());
    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(ExpectedLoadings[i], CalculatedLoadings[i], 1e-2);
        ;
    }
};

TEST(test_pta_pure, TestGetMultiplePureDeviationDRA)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    PotentialTheoryModels pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, co2);

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> experimental = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviation(pressures, experimental, 318.2, CO2DRAParams, "EABS");

    ASSERT_EQ(pressures.size(), experimental.size());

    EXPECT_NEAR(CalculatedDeviation, 0.0, 1e-4);
};

TEST(test_pta_pure, TestGetMultiplePureDeviationLEE)
{
    std::string potential = LEE_POTENTIAL;

    Adsorbent Z01x = Adsorbent("Z01x", 3.35, 0.382);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    PotentialTheoryModels pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, co2, Z01x);

    std::vector<double> CO2LEEParams = {125.63, 12.26, 765.70}; // eps/K, L, A

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> experimental = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviation(pressures, experimental, 318.2, CO2LEEParams, "EABS");

    ASSERT_EQ(pressures.size(), experimental.size());

    EXPECT_NEAR(CalculatedDeviation, 2.5183, 1e-4);
};

TEST(test_pta_pure, TestGetMultiplePureDeviationSTEELE)
{
    std::string potential = STEELE_POTENTIAL;

    Fluid fluid;
    fluid.CriticalPressure = 73.773e5;
    fluid.CriticalTemperature = 304.13;
    fluid.AccentricFactor = 0.22394;
    fluid.LennardJonnesDiameter = 3.941;

    Adsorbent adsorbent = Adsorbent("Z01x", 3.35, 0.382);

    PotentialTheoryModels pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, fluid, adsorbent);

    std::vector<double> params = {109.32, 13.34, 611.88}; // eps/K, L, A

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> experimental = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviation(pressures, experimental, 318.2, params, "EABS");

    ASSERT_EQ(pressures.size(), experimental.size());

    EXPECT_NEAR(CalculatedDeviation, 41.7384, 1e-4);
};

// // equations of state

TEST(test_pta_pure, TestGetMultiplePureDeviationDRAWithSRK)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    PotentialTheoryModels pta_model = PotentialTheoryModels(potential, "srk", "excess", 155, co2);

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> expected = {0.076485, 0.38631, 0.546913, 0.795782, 1.013569, 3.157997, 4.30760, 5.46686, 5.83229, 6.119404, 6.16664, 6.12293};

    std::vector<double> calculated = pta_model.GetPureLoadings(pressures, 318.2, CO2DRAParams);

    ASSERT_EQ(pressures.size(), calculated.size());

    for (std::size_t i = 0; i < pressures.size(); i++)
    {
        EXPECT_NEAR(calculated[i], expected[i], 1e-4);
    }
};

TEST(test_pta_mixture, testMixtureCalculation)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    auto pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, {co2, ch4});

    std::vector<double> CO2_DRA_params = {7880.19, 0.29, 2.};
    std::vector<double> CH4_DRA_params = {5600, 0.36, 3.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(pressures[7], 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationDRASRK)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    auto pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, {co2, ch4});

    std::vector<double> CO2_DRA_params = {7880.19, 0.29, 2.};
    std::vector<double> CH4_DRA_params = {5600, 0.36, 3.};

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(101000.0, 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationDRASRKSameAsFront)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    auto pta_model = PotentialTheoryModels(potential, "srk", "excess", 555, {co2, ch4});

    std::vector<double> CO2_DRA_params = {7880.19, 0.29, 2.};
    std::vector<double> CH4_DRA_params = {5600, 0.36, 3.};

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(101000.0, 310.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationLeeSRK)
{
    std::string potential = LEE_POTENTIAL;

    Adsorbent Z01x = Adsorbent("Z01x", 3.35, 0.382);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;
    ch4.LennardJonnesDiameter = 3.789;

    auto pta_model = PotentialTheoryModels(potential, "srk", "excess", 155, {co2, ch4}, Z01x);

    std::vector<double> CO2_DRA_params = {109.32, 13.34, 611.8};
    std::vector<double> CH4_DRA_params = {92.32, 13.34, 456.8}; // eps_k, L, A

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(101000.0, 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationSteeleSRK)
{
    std::string potential = STEELE_POTENTIAL;

    Adsorbent Z01x = Adsorbent("Z01x", 3.35, 0.382);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;
    ch4.LennardJonnesDiameter = 3.789;

    std::vector<double> CO2_DRA_params = {109.32, 14.90, 611.8};
    std::vector<double> CH4_DRA_params = {92.32, 14.90, 456.8}; // eps_k, L, A

    auto pta_model = PotentialTheoryModels(potential, "srk", "excess", 155, {co2, ch4}, Z01x);

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(101000.0, 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationLee)
{
    std::string potential = LEE_POTENTIAL;

    Adsorbent Z01x = Adsorbent("Z01x", 3.35, 0.382);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;
    ch4.LennardJonnesDiameter = 3.789;

    auto pta_model = PotentialTheoryModels(potential, "pr77", "excess", 555, {co2, ch4}, Z01x);

    std::vector<double> CO2_DRA_params = {109.32, 13.34, 611.8};
    std::vector<double> CH4_DRA_params = {92.32, 13.34, 456.8}; // eps_k, L, A

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(pressures[7], 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationSteele)
{
    std::string potential = STEELE_POTENTIAL;

    Adsorbent Z01x = Adsorbent("Z01x", 3.35, 0.382);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;
    ch4.LennardJonnesDiameter = 3.798;

    std::vector<double> CO2_DRA_params = {109.32, 13.34, 611.8};
    std::vector<double> CH4_DRA_params = {92.32, 13.34, 456.8}; // eps_k, L, A

    auto pta_model = PotentialTheoryModels(potential, "pr77", "excess", 555, {co2, ch4}, Z01x);

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(pressures[7], 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationThreeComponents)
{
    std::string potential = DRA_POTENTIAL;

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    Fluid n2;
    n2.CriticalPressure = 33.958e5;
    n2.CriticalTemperature = 126.192;
    n2.AccentricFactor = 0.0372;

    std::vector<double> CO2_DRA_params = {7880.19, 0.29, 2.};
    std::vector<double> CH4_DRA_params = {7000.19, 0.2, 3.};
    std::vector<double> N2_DRA_params = {5000.5, 0.57, 2.};

    auto pta_model = PotentialTheoryModels(potential, "pr77", "excess", 155, {co2, ch4, n2});

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double> bulk_composition_fraction = {0.25, 0.25, 0.50};

    std::vector<double>
        ns = pta_model.GetMixtureLoading(pressures[11], 318.2, bulk_composition_fraction, {CO2_DRA_params, CH4_DRA_params, N2_DRA_params});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}