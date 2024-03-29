#include <gtest/gtest.h>
#include "../isotherms/vsm/vsm.h"

double getVSMLoading(std::string isotherm, std::vector<double> parameters)
{
    VacancySolutionMethod VSM = VacancySolutionMethod(isotherm);
    return VSM.GetPureLoading(100, 298.15, parameters);
}

std::vector<double> getVSMLoadings(std::string isotherm, std::vector<double> parameters)
{
    VacancySolutionMethod VSM = VacancySolutionMethod(isotherm);
    return VSM.GetPureLoadings({10, 50, 100}, 298.15, parameters);
}

TEST(test_vsm, GetLoading)
{
    EXPECT_NEAR(getVSMLoading("wilson", {10, 0.05, 0.3, 0.5}), 8.90614, 1e-5);
    EXPECT_NEAR(getVSMLoading("nrtl", {10, 0.05, 0.1, 0.2, 0.5}), 9.30121, 1e-5);
    EXPECT_NEAR(getVSMLoading("flory-huggins", {10, 5e-2, 0.9}), 2.929822, 1e-5);
}

TEST(test_vsm, GetLoadingsWilson)
{
    std::vector<double> ExpectedLoadings = {0.64428, 7.23952, 8.906141};
    std::vector<double> CalculatedLoadings = getVSMLoadings("wilson", {10, 0.05, 0.3, 0.5});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    }
}

TEST(test_vsm, GetLoadingsNRTL)
{
    std::vector<double> ExpectedLoadings = {5.8685772, 8.7056045, 9.301218};
    std::vector<double> CalculatedLoadings = getVSMLoadings("nrtl", {10, 0.05, 0.1, 0.2, 0.5});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    }
}

TEST(test_vsm, GetLoadingsFloryHuggins)
{
    std::vector<double> ExpectedLoadings = {0.460217, 1.806217, 2.92982};
    std::vector<double> CalculatedLoadings = getVSMLoadings("flory-huggins", {10, 5e-2, 0.9});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    }
}

TEST(test_vsm, GetDeviation)
{
    std::vector<double> ExperimetalLoadings = {1.64428, 8.23952, 9.906141};

    VacancySolutionMethod VSM = VacancySolutionMethod("wilson");

    double Deviation = VSM.GetDeviation({10, 50, 100}, ExperimetalLoadings, 298.15, {10, 0.05, 0.3, 0.5}, "ARE");

    EXPECT_NEAR(Deviation, 27.682714, 1e-3);
}

TEST(test_vsm_mixture, GetLoadingWilson)
{
    std::vector<double> ExpectedLoadings = {1.971666, 1.877777};

    VacancySolutionMethod VSM = VacancySolutionMethod("wilson");

    std::vector<double> CalculatedLoadings = VSM.GetMixtureLoading(10, 298.15, {0.5, 0.5}, {{10, 0.05, 0.3, 0.5}, {3, 0.05, 0.3, 0.5}});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    };
}

TEST(test_vsm_mixture, GetLoadingNRTL)
{
    std::vector<double> ExpectedLoadings = {0.177909, 0.215217};

    VacancySolutionMethod VSM = VacancySolutionMethod("nrtl");

    std::vector<double> CalculatedLoadings = VSM.GetMixtureLoading(10, 298.15, {0.5, 0.5}, {{2, 0.05, 0.1, 0.2, 0.5}, {3, 0.05, 0.1, 0.2, 0.5}});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    };
}

TEST(test_vsm_mixture, GetLoadingFloryHuggins)
{
    std::vector<double> ExpectedLoadings = {0.1191434, 0.4873899};

    VacancySolutionMethod VSM = VacancySolutionMethod("flory-huggins");

    std::vector<double> CalculatedLoadings = VSM.GetMixtureLoading(10, 298.15, {0.5, 0.5}, {{2, 0.05, 0.1}, {3, 0.05, 0.1}});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    };
}