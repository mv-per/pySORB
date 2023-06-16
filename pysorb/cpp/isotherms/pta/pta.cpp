#include "pta.h"

std::function<double(double, double, std::vector<double>)> PotentialTheoryModels::GetPureLoadingInvoker()
{
    if (this->Potential == DRA_POTENTIAL)
    {
        return [this](double P, double T, std::vector<double> params)
        { return GetDRAPureLoading(P, T, params, this->MonoEosInvoker, this->MonoPotentialInvoker, this->fluid, this->NumberOfLayers, this->IsothermType, this->Potential); };
    }
    else if (this->Potential == STEELE_POTENTIAL || this->Potential == LEE_POTENTIAL)
    {
        return [this](double P, double T, std::vector<double> params)
        { 
            if (!this->AdsorbentConfigured) {
                throw std::invalid_argument("Adsorbent properties are needed for LJ-based potentials and is not defined.");
            }
            return GetLJPureLoading(P, T, params, this->MonoEosInvoker, this->MonoPotentialInvoker, this->fluid, this->NumberOfLayers, this->IsothermType, this->Potential); };
    }
    else
    {
        throw std::invalid_argument("Isotherm potential not found/defined.");
    }
}

void PotentialTheoryModels::SetAdsorbent(Adsorbent adsorbent)
{
    this->adsorbent = adsorbent;
    this->AdsorbentConfigured = true;
}

void PotentialTheoryModels::SetupInvokers(bool pure)
{

    if (pure)
    {
        this->MonoEosInvoker = GetPureEquationOfStateInvoker(this->EquationOfState, this->fluid);
        this->MonoPotentialInvoker = GetPureAdsorptionPotentialInvoker(this->Potential, this->fluid, this->adsorbent);
        this->PureLoadingInvoker = this->GetPureLoadingInvoker();
    }
    else
    {
        assert(this->Fluids.size() > 1);
        // this->MixEosInvoker = '22';
    }
}