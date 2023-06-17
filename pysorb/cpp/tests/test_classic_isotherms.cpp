#include <gtest/gtest.h>
#include "../isotherms/empirical/classic_isotherms.h"

#include <iostream>

double getLoading(std::string isotherm, std::vector<double> parameters)
{
    ClassicIsotherms IsothermModel = ClassicIsotherms(isotherm);
    return IsothermModel.GetPureLoading(10, 298.15, parameters);
}

std::vector<double> getLoadings(std::string isotherm, std::vector<double> parameters)
{
    ClassicIsotherms IsothermModel = ClassicIsotherms(isotherm);
    return IsothermModel.GetPureLoadings({10, 20, 30, 40}, 298.15, parameters);
}

TEST(test_classics, GetLoading)
{
    EXPECT_NEAR(getLoading("jensen-seaton", {10, 0.05, 0.3, 0.5}), 0.1832436, 1e-5);
    EXPECT_NEAR(getLoading("langmuir", {10, 0.05}), 3.33333, 1e-5);
    EXPECT_NEAR(getLoading("dual-langmuir", {10, 0.05, 3, 0.2}), 5.33333, 1e-5);
    EXPECT_NEAR(getLoading("unilan", {10, 0.05, 3, 0.2}), 3.96197, 1e-5);
    EXPECT_NEAR(getLoading("keller-staudt-toth", {10, 0.05, 3, 0.2}), 5.60199, 1e-5);
    EXPECT_NEAR(getLoading("freundlich", {2, 1.2}), 13.62584, 1e-5);
    EXPECT_NEAR(getLoading("sips", {10, 0.05, 0.3}), 0.902579, 1e-5);
    EXPECT_NEAR(getLoading("toth", {10, 0.05, 0.3}), 0.689034, 1e-5);
    EXPECT_NEAR(getLoading("freundlich-2", {2, 0.5, -15000}), 1.4847480, 1e-5);
    EXPECT_NEAR(getLoading("redlich-peterson", {2, 1.2, 2}), 0.165289, 1e-5);
}

TEST(test_classics, IsothermNotFound)
{

    EXPECT_THROW({
        try
        {
            getLoading("unknown-isotherm", {2, 1.2, 2});
        }
        catch (const std::runtime_error &e)
        {
            // and this tests that it has the correct message
            EXPECT_STREQ("PureLoadingInvoker not defined", e.what());
            throw;
        }
    },
                 std::runtime_error);
}

TEST(test_classics, GetLoadings)
{
    std::vector<double> ExpectedLoadings = {3.3333, 5, 6, 6.6666};
    std::vector<double> CalculatedLoadings = getLoadings("langmuir", {10, 0.05});

    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(CalculatedLoadings[i], ExpectedLoadings[i], 1e-4);
    }
}

TEST(test_classics, GetDeviation)
{
    ClassicIsotherms IsothermModel = ClassicIsotherms("langmuir");

    std::vector<std::string> deviationFunctions = {"ARE", "SSE", "EABS", "HYBRID", "MPSD", "SORE", "CHI_2", "R_S"};
    std::vector<double> expectedDeviations = {-77.20797, 41.56132, 11.11965, 189.48791, 296.04621, 75.15154, 8.76433, -19.78066};

    for (std::size_t i = 0; i < deviationFunctions.size(); i++)
    {
        EXPECT_NEAR(IsothermModel.GetDeviation({10, 20, 30}, {4, 5, 6}, 298.15, {10, 0.4}, deviationFunctions[i]), expectedDeviations[i], 1e-4);
    }
}

TEST(test_classics, GetMixtureLoading)
{
    ClassicIsotherms IsothermModel = ClassicIsotherms("extended-langmuir");

    std::vector<double> loadings = IsothermModel.GetMixtureLoading(10, 298.15, {0.3, 0.7}, {{10, 0.05}, {6, 0.08}});
    std::vector<double> expectedLoadings = {0.631578, 2.3578};
    for (std::size_t i = 0; i < loadings.size(); i++)
    {
        EXPECT_NEAR(loadings[i], expectedLoadings[i], 1e-4);
    }
}
