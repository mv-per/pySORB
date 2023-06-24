#include "extended_isotherms.h"

std::vector<double> extended_langmuir(double Pressure, std::vector<double> BulkComposition, std::vector<std::vector<double>> Parameters)
{

    std::size_t i;
    std::size_t ncomp = BulkComposition.size();

    assert(ncomp == Parameters.size());

    CheckCompositionFraction(BulkComposition);

    std::vector<double> Loadings(ncomp, 0.0);

    double b_mix = 0.;
    double n_mix = 0.;
    for (i = 0; i < ncomp; i++)
    {
        b_mix += Parameters[i][1] * Pressure * BulkComposition[i];
        n_mix += Parameters[i][0] * BulkComposition[i];
    }

    for (i = 0; i < ncomp; i++)
    {
        Loadings[i] = n_mix * Parameters[i][1] * Pressure * BulkComposition[i] / (1.0 + b_mix);
    }

    return Loadings;
}
