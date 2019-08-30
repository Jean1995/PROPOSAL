
#include <functional>

#include "PROPOSAL/crossection/BhabhaMollerIntegral.h"
#include "PROPOSAL/crossection/parametrization/BhabhaMoller.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Logging.h"
#include <cmath>

using namespace PROPOSAL;

BhabhaMollerIntegral::BhabhaMollerIntegral(const BhabhaMoller& param)
        : CrossSectionIntegral(GetType(param), param)
{
}

BhabhaMollerIntegral::BhabhaMollerIntegral(const BhabhaMollerIntegral& bmi)
        : CrossSectionIntegral(bmi)
{
}

BhabhaMollerIntegral::~BhabhaMollerIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //
double BhabhaMollerIntegral::CalculatedEdxWithoutMultiplier(double energy){
    double sum = 0;

    for (int i = 0; i < (parametrization_->GetMedium().GetNumComponents()); i++)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        auto integrand_substitution = [&](double energy, double t){
            return std::exp(t) * parametrization_->FunctionToDEdxIntegral(energy, 1 - std::exp(t));
        };

        double t_min = std::log(1. - limits.vUp);
        double t_max = std::log(1. - limits.vMin);

        sum += dedx_integral_.Integrate(
                limits.vMin,
                limits.vUp,
                std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
                4);

/*
        sum += dedx_integral_.Integrate(
                t_min,
                t_max,
                std::bind(integrand_substitution, energy, std::placeholders::_1),
                2);
*/
    }

    return energy * sum;
}

double BhabhaMollerIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * BhabhaMollerIntegral::CalculatedEdxWithoutMultiplier(energy);
}

DynamicData::Type BhabhaMollerIntegral::GetType(const BhabhaMoller& param){
    if(param.GetParticleDef().charge > 0.){
        return DynamicData::Bhabha;
    }
    else if(param.GetParticleDef().charge < 0.){
        return DynamicData::Moller;
    }
    else{
        log_fatal("Bhabha/Moller scattering not possible for uncharged particles");
        return DynamicData::None;
    }
}
