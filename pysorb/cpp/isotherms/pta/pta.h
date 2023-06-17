#ifndef PTA_H
#define PTA_H

#include <string>
#include <vector>
#include <functional>
#include <cassert>
#include "../utils.h"
#include "../_base_isotherm_model.h"
#include "../equations_of_state/eos.h"
#include "../data_classes.h"
#include "pta_pure.h"

class PotentialTheoryModels : public BaseIsothermModel
{
public:
    /**
     * @brief Constructs a PotentialTheoryModels object and initializes the PureLoadingInvoker based on the given isotherm name.
     * @param model The name of the activity coefficient model to be used, options: `wilson`, `nrtl`, `flory-huggins`.
     */
    PotentialTheoryModels(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, std::vector<Fluid> fluids, Adsorbent adsorbent) : BaseIsothermModel()
    {
        this->Potential = potential;
        this->IsothermType = isotherm_type;
        this->NumberOfLayers = num_of_layers;
        this->EquationOfState = equation_of_state;
        this->fluids = fluids;
        this->SetAdsorbent(adsorbent);
        this->IsPure = false;
        this->SetupInvokers();
    }

    PotentialTheoryModels(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, std::vector<Fluid> fluids) : BaseIsothermModel()
    {
        this->Potential = potential;
        this->IsothermType = isotherm_type;
        this->NumberOfLayers = num_of_layers;
        this->EquationOfState = equation_of_state;
        this->fluids = fluids;
        this->IsPure = false;

        this->SetupInvokers();
    }

    PotentialTheoryModels(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, Fluid fluid, Adsorbent adsorbent) : BaseIsothermModel()
    {
        this->Potential = potential;
        this->IsothermType = isotherm_type;
        this->NumberOfLayers = num_of_layers;
        this->EquationOfState = equation_of_state;
        this->fluid = fluid;
        this->SetAdsorbent(adsorbent);
        this->IsPure = true;
        this->SetupInvokers();
    }

    PotentialTheoryModels(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers, Fluid fluid) : BaseIsothermModel()
    {
        this->Potential = potential;
        this->IsothermType = isotherm_type;
        this->NumberOfLayers = num_of_layers;
        this->EquationOfState = equation_of_state;
        this->fluid = fluid;
        this->IsPure = true;
        this->SetupInvokers();
    }

    void SetupInvokers() override;

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param model The name of the activity coefficient model to be used
     * @return The corresponding loading invoker function.
     * @throw std::invalid_argument If the model is not found or defined.
     */
    std::function<double(double, double, std::vector<double>)> GetPureLoadingInvoker(std::string potential) override;

    /**
     * @brief Returns the appropriate loading invoker function based on the given model name.
     * @param isotherm The name of the isotherm.
     * @return The corresponding isotherm invoker function.
     * @throw std::invalid_argument If the isotherm is not found or defined.
     */
    std::function<std::vector<double>(double, double, std::vector<double>, std::vector<std::vector<double>>)> GetMixtureLoadingInvoker(std::string potential) override;

private:
    bool IsPure = false;
    std::function<mono_eos(double, double)> MonoEosInvoker;
    std::function<mix_eos(std::vector<double>, double, double)> MixEosInvoker;
    std::function<double(double, std::vector<double>)> MonoPotentialInvoker;
    std::function<double(double, std::vector<double>)> MixPotentialInvoker;
    std::vector<Fluid> fluids;
    Fluid fluid;
    Adsorbent adsorbent;
    std::string Potential;
    bool AdsorbentConfigured = false;
    std::size_t NumberOfLayers;
    std::string IsothermType;
    std::string EquationOfState;

    void SetAdsorbent(Adsorbent adsorbent);
};

#endif