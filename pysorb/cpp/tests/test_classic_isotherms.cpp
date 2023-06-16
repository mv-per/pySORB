#include <gtest/gtest.h>
#include "../isotherms/empirical/classic_isotherms.h"

TEST(test_loadings, DemostrateDeviations)
{
    auto IsothermModel = ClassicIsotherms('langmuir');

    double loading = IsothermModel.GetLoading(1000, {10, 0.05});

    ASSERT_TRUE(loading == 0.9);
}