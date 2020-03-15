
#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/ComptonInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

ComptonInterpolant::ComptonInterpolant(const Compton& param, std::shared_ptr<const EnergyCutSettings> cuts, const InterpolationDef& def)
        : CrossSectionInterpolant(param, cuts)
{
    // Use own CrossSecition dNdx interpolation
    ComptonInterpolant::InitdNdxInterpolation(def);

    // --------------------------------------------------------------------- //
    // Builder for DEdx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    // Needed for CalculatedEdx integration
    ComptonIntegral compton(param, cuts);

    builder1d.SetMax(def.nodes_cross_section)
            .SetXMin(ME)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(true)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(true)
            .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedEdxWithoutMultiplier, &compton, std::placeholders::_1));

    builder_container.push_back(&builder1d);

    // --------------------------------------------------------------------- //
    // Builder for DE2dx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder_de2dx;
    Helper::InterpolantBuilderContainer builder_container_de2dx;

    builder_de2dx.SetMax(def.nodes_continous_randomization)
            .SetXMin(ME)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedE2dxWithoutMultiplier, &compton, std::placeholders::_1));

    builder_container_de2dx.push_back(&builder_de2dx);

    auto dedx_interpolant_vec = Helper::InitializeInterpolation("dEdx", builder_container, std::vector<Parametrization*>(1, parametrization_), def);
    dedx_interpolant_ = std::move(dedx_interpolant_vec.at(0));

    auto de2dx_interpolant_vec = Helper::InitializeInterpolation(
            "dE2dx", builder_container_de2dx, std::vector<Parametrization*>(1, parametrization_), def);
    de2dx_interpolant_ = std::move(de2dx_interpolant_vec.at(0));
}

/*ComptonInterpolant::ComptonInterpolant(const ComptonInterpolant& compton)
        : CrossSectionInterpolant(compton)
{
}*/

ComptonInterpolant::~ComptonInterpolant() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double ComptonInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}

//----------------------------------------------------------------------------//
double ComptonInterpolant::FunctionToBuildDNdxInterpolant2D(double energy,
                                                                 double v,
                                                                 Integral& integral,
                                                                 int component)
{
    parametrization_->SetCurrentComponent(component);
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

    double vUp = GetEnergyCut(energy);

    if (vUp == limits.vMax)
    {
        return 0;
    }

    v = vUp + (limits.vMax - vUp) * v;

    // Integrate with the substitution t = ln(1-v) to avoid numerical problems
    auto integrand_substitution = [&](double energy, double t){
        return std::exp(t) * parametrization_->FunctionToDNdxIntegral(energy, 1 - std::exp(t));
    };

    double t_min = std::log(1. - v);
    double t_max = std::log(1. - vUp);


    return integral.Integrate(
            t_min,
            t_max,
            std::bind(integrand_substitution, energy, std::placeholders::_1),
            2);
}

double ComptonInterpolant::CalculateCumulativeCrossSection(double energy, int component, double v)
{
    parametrization_->SetCurrentComponent(component);
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

    double vUp = GetEnergyCut(energy);

    v = (v - vUp) / (limits.vMax - vUp);

    return dndx_interpolant_2d_.at(component)->Interpolate(energy, v);
}

std::pair<double, double> ComptonInterpolant::StochasticDeflection(double energy, double energy_loss) {
    double theta_deflect = RandomGenerator::Get().RandomDouble() * 2 * PI; // random azimuth
    double cosphi = 1. - (ME * (1. / (energy - energy_loss) - 1. / energy));

    return std::make_pair(cosphi, theta_deflect);
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double ComptonInterpolant::CalculateStochasticLoss(double energy, double rnd1)
{

    double rnd;
    double rsum;
    double vUp;

    rnd  = rnd1 * sum_of_rates_;
    rsum = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        rsum += prob_for_component_[i];

        if (rsum > rnd)
        {
            parametrization_->SetCurrentComponent(i);
            Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

            vUp = GetEnergyCut(energy);

            if (vUp == limits.vMax)
            {
                return energy * vUp;
            }

            // Linear interpolation in v
            return energy * ( vUp + (limits.vMax - limits.vMin) * dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i]) );
        }
    }

    // sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero = true;
    for (size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

        vUp = GetEnergyCut(energy);

        if (vUp != limits.vMax)
            prob_for_all_comp_is_zero = false;
    }

    if (prob_for_all_comp_is_zero)
        return 0;

    log_fatal("sum was not initialized correctly");
    return 0; // just to prevent warnings
}

// ------------------------------------------------------------------------- //
void ComptonInterpolant::InitdNdxInterpolation(const InterpolationDef& def)
{
    // --------------------------------------------------------------------- //
    // Builder for dNdx
    // --------------------------------------------------------------------- //

    std::vector<Interpolant1DBuilder> builder1d(components_.size());
    std::vector<Interpolant2DBuilder> builder2d(components_.size());

    Helper::InterpolantBuilderContainer builder_container1d(components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(components_.size());

    Integral integral(IROMB, IMAXS, IPREC);

    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        // !!! IMPORTANT !!!
        // Order of builder matter because the functions needed for 1d interpolation
        // needs the already intitialized 2d interpolants.
        builder2d[i]
                .SetMax1(def.nodes_cross_section)
                .SetX1Min(ME)
                .SetX1Max(def.max_node_energy)
                .SetMax2(def.nodes_cross_section)
                .SetX2Min(1. / (2. * (1. - def.nodes_cross_section)))
                .SetX2Max((1. - 2. * def.nodes_cross_section) / (2. * (1. - def.nodes_cross_section)))
                .SetRomberg1(def.order_of_interpolation)
                .SetRational1(false)
                .SetRelative1(false)
                .SetIsLog1(true)
                .SetRomberg2(def.order_of_interpolation)
                .SetRational2(true)
                .SetRelative2(false)
                .SetIsLog2(false)
                .SetRombergY(def.order_of_interpolation)
                .SetRationalY(true)
                .SetRelativeY(false)
                .SetLogSubst(false)
                .SetFunction2D(std::bind(
                        &CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D,
                        this,
                        std::placeholders::_1,
                        std::placeholders::_2,
                        std::ref(integral),
                        i));

        builder_container2d[i] = &builder2d[i];

        builder1d[i]
                .SetMax(def.nodes_cross_section)
                .SetXMin(ME)
                .SetXMax(def.max_node_energy)
                .SetRomberg(def.order_of_interpolation)
                .SetRational(false)
                .SetRelative(false)
                .SetIsLog(true)
                .SetRombergY(def.order_of_interpolation)
                .SetRationalY(true)
                .SetRelativeY(false)
                .SetLogSubst(false)
                .SetFunction1D(std::bind(&CrossSectionInterpolant::FunctionToBuildDNdxInterpolant, this, std::placeholders::_1, i));

        builder_container1d[i]  = &builder1d[i];
    }


    dndx_interpolant_2d_ = Helper::InitializeInterpolation("dNdx_diff", builder_container2d, std::vector<Parametrization*>(1, parametrization_), def);
    dndx_interpolant_1d_ = Helper::InitializeInterpolation("dNdx", builder_container1d, std::vector<Parametrization*>(1, parametrization_), def);
}
