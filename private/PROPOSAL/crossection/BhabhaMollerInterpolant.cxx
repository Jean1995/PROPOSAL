
#include <functional>

#include "PROPOSAL/crossection/BhabhaMollerIntegral.h"
#include "PROPOSAL/crossection/BhabhaMollerInterpolant.h"
#include "PROPOSAL/crossection/parametrization/BhabhaMoller.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

BhabhaMollerInterpolant::BhabhaMollerInterpolant(const BhabhaMoller& param, InterpolationDef def)
        : CrossSectionInterpolant(GetType(param), param)
{
    // Use parent CrossSecition dNdx interpolation
    InitdNdxInterpolation(def);

    // --------------------------------------------------------------------- //
    // Builder for DEdx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    // Needed for CalculatedEdx integration
    BhabhaMollerIntegral bmi(param);

    builder1d.SetMax(def.nodes_cross_section)
            .SetXMin(param.GetParticleDef().mass)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(true)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(true)
            .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedEdxWithoutMultiplier, &bmi, std::placeholders::_1));

    builder_container.push_back(std::make_pair(&builder1d, &dedx_interpolant_));

    // --------------------------------------------------------------------- //
    // Builder for DE2dx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder_de2dx;
    Helper::InterpolantBuilderContainer builder_container_de2dx;

    builder_de2dx.SetMax(def.nodes_continous_randomization)
            .SetXMin(param.GetParticleDef().mass)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedE2dxWithoutMultiplier, &bmi, std::placeholders::_1));

    builder_container_de2dx.push_back(std::make_pair(&builder_de2dx, &de2dx_interpolant_));

    Helper::InitializeInterpolation("dEdx", builder_container, std::vector<Parametrization*>(1, parametrization_), def);
    Helper::InitializeInterpolation(
            "dE2dx", builder_container_de2dx, std::vector<Parametrization*>(1, parametrization_), def);
}

BhabhaMollerInterpolant::BhabhaMollerInterpolant(const BhabhaMollerInterpolant& bmi)
        : CrossSectionInterpolant(bmi)
{
}

BhabhaMollerInterpolant::~BhabhaMollerInterpolant() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double BhabhaMollerInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}

DynamicData::Type BhabhaMollerInterpolant::GetType(const BhabhaMoller& param){
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
