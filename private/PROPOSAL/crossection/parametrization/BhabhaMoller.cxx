
#include <functional>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/BhabhaMoller.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

/******************************************************************************
 *                               BhabhaMoller                                 *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

BhabhaMoller::BhabhaMoller(const ParticleDef& particle_def,
                               const Medium& medium,
                               const EnergyCutSettings& cuts,
                               double multiplier,
                               double threshold)
        : Parametrization(particle_def, medium, cuts, multiplier)
        , threshold_(threshold)
{
}

BhabhaMoller::BhabhaMoller(const BhabhaMoller& bbm)
        : Parametrization(bbm)
        , threshold_(bbm.threshold_)
{
}

BhabhaMoller::~BhabhaMoller() {}

bool BhabhaMoller::compare(const Parametrization& parametrization) const
{
    const BhabhaMoller* bbm = static_cast<const BhabhaMoller*>(&parametrization);
    if (threshold_ != bbm->threshold_)
        return false;
    else
        return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

size_t BhabhaMoller::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed, threshold_);

    return seed;
}

/******************************************************************************
 *                          Specifc Parametrizations                           *
 ******************************************************************************/

// ------------------------------------------------------------------------- //
// Bhabha Scattering
// ------------------------------------------------------------------------- //


BhabhaScattering::BhabhaScattering(const ParticleDef& particle_def,
                                               const Medium& medium,
                                               const EnergyCutSettings& cuts,
                                               double multiplier,
                                               double threshold)
        : BhabhaMoller(particle_def, medium, cuts, multiplier, threshold)
{
}

BhabhaScattering::BhabhaScattering(const BhabhaScattering& bs)
        : BhabhaMoller(bs)
{
}

BhabhaScattering::~BhabhaScattering()
{
}

const std::string BhabhaScattering::name_ = "BhabhaScattering";

Parametrization::IntegralLimits BhabhaScattering::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = std::min(threshold_/energy, 1.0);


    limits.vMax = 1. - ME / energy;

    if (limits.vMax < 0) {
        limits.vMax = 0;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    if (limits.vMax < limits.vUp)
    {
        limits.vMax = limits.vUp;
    }

    return limits;
}

double BhabhaScattering::DifferentialCrossSection(double energy, double v)
{
    double aux    = 0;

    double gamma = energy/ME;
    double epsilon = (v * energy)/(energy - ME);
    double betasquared = 1. - 1. / (gamma * gamma);
    double y = 1. / (gamma + 1.);
    double B1 = 2. - y * y;
    double B2 = (1. - 2. * y) * (3. + y * y);
    double B3 = std::pow(1. - 2. * y, 2.) + std::pow(1. - 2. * y, 3.);
    double B4 = std::pow(1. - 2. * y, 3.);

    aux = 1. / (betasquared * epsilon * epsilon) - B1 / epsilon + B2 - B3 * epsilon + B4 * epsilon * epsilon;

    aux *= 1. / (gamma - 1.);
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI * std::pow(RE, 2.) * components_[component_index_]->GetNucCharge();

    return aux;

}


// ------------------------------------------------------------------------- //
// Moller Scattering
// ------------------------------------------------------------------------- //


MollerScattering::MollerScattering(const ParticleDef& particle_def,
                                   const Medium& medium,
                                   const EnergyCutSettings& cuts,
                                   double multiplier,
                                   double threshold)
        : BhabhaMoller(particle_def, medium, cuts, multiplier, threshold)
{
}

MollerScattering::MollerScattering(const MollerScattering& bs)
        : BhabhaMoller(bs)
{
}

MollerScattering::~MollerScattering()
{
}

const std::string MollerScattering::name_ = "MollerScattering";

Parametrization::IntegralLimits MollerScattering::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    limits.vMin = std::min(threshold_/energy, 1.0);

    // ambiguity of initial and final electron: for v = v_max, both electrons have the same energy
    limits.vMax = 0.5 * (1. - ME / energy);

    if (limits.vMax < 0)
    {
        limits.vMax = 0;
    }

    limits.vUp = std::min(limits.vMax, cut_settings_.GetCut(energy));

    if (limits.vUp < limits.vMin)
    {
        limits.vUp = limits.vMin;
    }

    if (limits.vMax < limits.vUp)
    {
        limits.vMax = limits.vUp;
    }

    return limits;
}

double MollerScattering::DifferentialCrossSection(double energy, double v)
{
    double aux    = 0;

    double gamma = energy/ME;
    double epsilon = (v * energy)/(energy - ME);
    double epsilonprime = 1 - epsilon;
    double C1 = std::pow((gamma - 1.)/gamma, 2.);
    double C2 = (2. * gamma - 1.) / (gamma * gamma);
    double betasquared = 1. - 1. / (gamma * gamma);

    aux = C1 + 1. / epsilon * (1. / epsilon - C2) + 1. / epsilonprime * (1. / epsilonprime - C2);

    aux *= 1. / (betasquared * (gamma - 1.));
    aux *= 1. / (1. - 1. / gamma); // conversion from epsilon to v
    aux *= 2. * PI * std::pow(RE, 2.) * components_[component_index_]->GetNucCharge();

    return aux;

}

