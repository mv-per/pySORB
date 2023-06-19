#include "pta.h"

std::function<double(double, double, std::vector<double>)> PotentialTheoryModels::GetPureLoadingInvoker(std::string potential)
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

std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> PotentialTheoryModels::GetMixtureLoadingInvoker(std::string potential)
{
    if (this->Potential == DRA_POTENTIAL)
    {
        return [this](double P, double T, std::vector<double> BulkComposition, std::vector<std::vector<double>> params)
        { return GetDRAMixtureLoading(P, T, BulkComposition, params, this->fluids, this->NumberOfLayers, this->IsothermType, this->Potential, this->EquationOfState, this->adsorbent); };
    }
    else if (this->Potential == STEELE_POTENTIAL || this->Potential == LEE_POTENTIAL)
    {
        return [this](double P, double T, std::vector<double> BulkComposition, std::vector<std::vector<double>> params)
        {
            if (!this->AdsorbentConfigured) {
                throw std::invalid_argument("Adsorbent properties are needed for LJ-based potentials and is not defined.");
            }
            return GetLJMixtureLoading(P, T, BulkComposition, params, this->fluids, this->NumberOfLayers, this->IsothermType, this->Potential, this->EquationOfState, this->adsorbent); };
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

void PotentialTheoryModels::SetupInvokers()
{

    if (this->IsPure)
    {
        this->MonoEosInvoker = GetPureEquationOfStateInvoker(this->EquationOfState, this->fluid);
        this->MonoPotentialInvoker = GetPureAdsorptionPotentialInvoker(this->Potential, this->fluid, this->adsorbent);
        this->PureLoadingInvoker = this->GetPureLoadingInvoker(this->Potential);
    }
    else
    {
        assert(this->fluids.size() > 1);
        this->MixtureLoadingInvoker = GetMixtureLoadingInvoker(this->Potential);
    }
}