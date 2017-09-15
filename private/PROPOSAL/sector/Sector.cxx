#include <boost/bind.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
// #include "PROPOSAL/continous_randomization/ContinuousRandomization.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/scattering/ScatteringDefault.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/sector/Sector.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                                 Sector                                 *
******************************************************************************/

Sector::Definition::Definition()
    : PropagationUtility::Definition()
    , do_interpolation(true)
    , do_weighting(false)
    , weighting_order(0)
    , do_continuous_randomization(false)
    , location(Sector::ParticleLocation::InsideDetector)
{
}

Sector::Definition::~Definition()
{
}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

// Standard constructor
Sector::Sector(PROPOSALParticle& particle)
    : ini_(0)
    , sector_def_()
    , weighting_starts_at_(0)
      //TODO(mario): init different Fri 2017/09/01
    , particle_(particle)
    , geometry_(new Sphere(Vector3D(), 1e18, 0))
    // , randomizer_(NULL)
    , utility(new PropagationUtilityInterpolant(particle.GetParticleDef(), Water(), EnergyCutSettings(), sector_def_))
    , scattering_(new ScatteringDefault())
{
}

Sector::Sector(PROPOSALParticle& particle, const Medium& medium,
                       const Geometry& geometry,
                       const EnergyCutSettings& cut_settings,
                       const Scattering& scattering,
                       const Definition& def)
    : ini_(0)
    , sector_def_(def)
    , weighting_starts_at_(0)
    , particle_(particle)
    , geometry_(geometry.clone())
    // , randomizer_(NULL)
    , utility(NULL)
    , scattering_(scattering.clone())
{
    if (def.do_interpolation)
    {
        utility = new PropagationUtilityInterpolant(particle_.GetParticleDef(), medium, cut_settings, def);
    }
    else
    {
        utility = new PropagationUtilityIntegral(particle_.GetParticleDef(), medium, cut_settings, def);
    }

    //TODO(mario): Polymorphic initilaization in collections childs  Sun 2017/08/27
    // if (sector_def_.do_continuous_randomization)
    // {
    //     randomizer_ = new ContinuousRandomization(particle);
    // }

    // if (sector_def_.do_scattering)
    // {
    //     scattering_ = ScatteringFactory::Get().CreateScattering(sector_def_.scattering_model);
    // }
}

Sector::Sector(const Sector& collection)
    :ini_(collection.ini_)
    ,sector_def_(collection.sector_def_)
    ,weighting_starts_at_(collection.weighting_starts_at_)
    ,particle_(collection.particle_)
    ,geometry_(collection.geometry_->clone())
    // ,randomizer_(collection.randomizer_) //TODO(mario): ranomizer clone Sat 2017/08/26
    ,utility(collection.utility->clone())
    ,scattering_(collection.scattering_->clone())
{
}

Sector::~Sector()
{
    delete geometry_;
    delete utility;
    delete scattering_;

    // if (randomizer_)
    // {
    //     delete randomizer_;
    // }

    //TODO(mario): delete scatter Sat 2017/08/26
}

// ------------------------------------------------------------------------- //
double Sector::Propagate(double distance)
{
    bool flag;
    double displacement;

    double propagated_distance = 0;

    double initial_energy = particle_.GetEnergy();
    double final_energy   = particle_.GetEnergy();

    bool particle_interaction = false;

    pair<double, ParticleType::Enum> decay;
    std::vector<PROPOSALParticle*> decay_products;

    pair<double, DynamicData::Type> energy_loss;

    // TODO(mario): check Fri 2017/08/25
    // int secondary_id    =   0;

    // first: final energy befor first interaction second: energy at which the
    // particle decay
    // first and second are compared to decide if interaction happens or decay
    pair<double, double> energy_till_stochastic_;

    if (distance < 0)
    {
        distance = 0;
    }

    if (initial_energy <= particle_.GetLow() || distance == 0)
    {
        flag = false;
    } else
    {
        flag = true;
    }

    while (flag)
    {
        energy_till_stochastic_ = CalculateEnergyTillStochastic( initial_energy);
        if (energy_till_stochastic_.first > energy_till_stochastic_.second)
        {
            particle_interaction = true;
            final_energy         = energy_till_stochastic_.first;
        } else
        {
            particle_interaction = false;
            final_energy         = energy_till_stochastic_.second;
        }

        // Calculate the displacement according to initial energy initial_energy and final_energy
        displacement = utility->CalculateDisplacement(initial_energy,
                                                      final_energy,
                                                      utility->GetMedium().GetDensityCorrection() *
                                                          (distance - propagated_distance)) /
                       utility->GetMedium().GetDensityCorrection();

        // The first interaction or decay happens behind the distance we want to propagate
        // So we calculate the final energy using only continuous losses
        if (displacement > distance - propagated_distance)
        {
            displacement = distance - propagated_distance;

            final_energy = utility->CalculateFinalEnergy( initial_energy, utility->GetMedium().GetDensityCorrection() * displacement);
        }
        // Advance the Particle according to the displacement
        // Initial energy and final energy are needed if Molier Scattering is enabled
        AdvanceParticle( displacement, initial_energy, final_energy);

        propagated_distance += displacement;

        if (abs(distance - propagated_distance) < abs(distance) * COMPUTER_PRECISION)
        {
            propagated_distance = distance; // computer precision control
        }

        //TODO(mario): Revert randomizer Fri 2017/08/25
        // Randomize the continuous energy loss if this option is enabled
        // if (sector_def_.do_continuous_randomization)
        // {
        //     if (final_energy != particle_.GetLow())
        //     {
        //         double rnd = RandomGenerator::Get().RandomDouble();
        //         final_energy = randomizer_->Randomize(crosssections_, initial_energy, final_energy, rnd);
        //     }
        // }

        // Lower limit of particle energy is reached or
        // or complete particle is propagated the whole distance
        if (final_energy == particle_.GetLow() || propagated_distance == distance)
        {
            break;
        }

        // Set the particle energy to the current energy before making
        // stochatic losses or decay
        particle_.SetEnergy(final_energy);

        if (particle_interaction)
        {
            energy_loss = MakeStochasticLoss();
            if (energy_loss.second == DynamicData::None)
            {
                // in this case, no cross section is chosen, so there is no interaction
                // due to the parameterization of the cross section cutoffs
                log_debug("no interaction due to the parameterization of the cross section cutoffs. final energy: %f\n",
                          final_energy);
                initial_energy = final_energy;
                continue;
            }
            final_energy -= energy_loss.first;
            // log_debug("Energyloss: %f\t%s", energy_loss.first,
            // PROPOSALParticle::GetName(energy_loss.second).c_str());
            // //TODO(mario): hack Thu 2017/08/24
            Output::getInstance().FillSecondaryVector(particle_, energy_loss.second, energy_loss.first, 0);
            // secondary_id    =   particle_.GetParticleId() + 1;
            // Output::getInstance().FillSecondaryVector(& secondary_id, energy_loss, 0);
        } else
        {
            energy_loss = MakeStochasticLoss();
            DecayChannel* mode = &particle_.GetDecayTable().SelectChannel();
            decay_products     = mode->Decay(&particle_);
            Output::getInstance().FillSecondaryVector(decay_products);

            // TODO(mario): Delete decay products Tue 2017/08/22

            // decay           =   current_collection_->MakeDecay();
            // final_energy    =   0;
            // log_debug("Decay of particle: %s", particle_->GetName().c_str());
            // secondary_id    = particle_->GetParticleId()  +   1;
            // Output::getInstance().FillSecondaryVector(particle_, secondary_id, decay ,0);
        }

        // break if the lower limit of particle energy is reached
        if (final_energy <= particle_.GetLow())
        {
            break;
        }

        // Next round: update the inital energy
        initial_energy = final_energy;
    }

    // if(stopping_decay_)
    // {
    //     if(propagated_distance!=distance && final_energy!=0 && particle_->GetLifetime()>=0)
    //     {
    //         particle_->SetEnergy(particle_->GetMass());
    //
    //         double t    =   particle_->GetTime() -particle_->GetLifetime()*log(RandomDouble());
    //         double product_energy   =   0;
    //
    //         pair<double, ParticleType::Enum> decay_to_store;
    //         secondary_id    =   particle_->GetParticleId() + 1;
    //
    //         particle_->SetTime( t );
    //
    //         if(particle_->GetType()==2)
    //         {
    //             // --------------------------------------------------------------------- //
    //             // Calculate random numbers before passing to a fuction, because
    //             // the order of argument evaluation is unspecified in c++ standards and
    //             // therfor depend on the compiler.
    //             // --------------------------------------------------------------------- //
    //
    //             double rnd1 = RandomDouble();
    //             double rnd2 = RandomDouble();
    //
    //             product_energy  =   current_collection_->GetDecay()->CalculateProductEnergy(rnd1, 0.5, rnd2);
    //         }
    //         else
    //         {
    //             product_energy  =   current_collection_->GetDecay()->CalculateProductEnergy(RandomDouble(), 0.5,
    //             0.5);
    //         }
    //
    //         decay_to_store.first    =   product_energy;
    //         decay_to_store.second   =   current_collection_->GetDecay()->GetOut();
    //
    //         final_energy  =   0;
    //
    //         Output::getInstance().FillSecondaryVector(particle_,secondary_id, decay_to_store, final_energy);
    //     }
    // }

    particle_.SetEnergy(final_energy);

    // Particle reached the border, final energy is returned
    if (propagated_distance == distance)
    {
        return final_energy;
    }
    // The particle stopped/decayed, the propageted distance is return with a minus sign
    else
    {
        return -propagated_distance;
    }
    // Should never be here
    return 0;
}

std::pair<double, double> Sector::CalculateEnergyTillStochastic(
                                                                    double initial_energy)
{
    double rndd = -log(RandomGenerator::Get().RandomDouble());
    double rndi = -log(RandomGenerator::Get().RandomDouble());

    double rndiMin = 0;
    double rnddMin = 0;

    pair<double, double> final;

    // solving the tracking integral
    if (particle_.GetLifetime() < 0)
    {
        rnddMin = 0;
    } else
    {
        rnddMin = utility->CalculateTrackingIntegal( initial_energy, rndd, false) / utility->GetMedium().GetDensityCorrection();
    }

    rndiMin = utility->CalculateTrackingIntegal( initial_energy, rndi, true);

    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0)
    {
        final.second = particle_.GetLow();
    } else
    {
        final.second = utility->CalculateFinalEnergy( initial_energy, rndd * utility->GetMedium().GetDensityCorrection(), false);
    }

    if (rndi >= rndiMin || rndiMin <= 0)
    {
        final.first = particle_.GetLow();
    } else
    {
        final.first = utility->CalculateFinalEnergy( initial_energy, rndi, true);
    }

    return final;
}

void Sector::AdvanceParticle(double dr, double ei, double ef)
{

    double dist       = particle_.GetPropagatedDistance();
    double time       = particle_.GetTime();
    Vector3D position = particle_.GetPosition();

    dist += dr;

    if (sector_def_.do_exact_time_calculation)
    {
        time += utility->CalculateParticleTime( ei, ef) / utility->GetMedium().GetDensityCorrection();
    } else
    {
        time += dr / SPEED;
    }

    // TODO(mario): Adjucst the whole scatteing class Thu 2017/08/24
    scattering_->Scatter(particle_, utility->GetCrosssections(), dr, ei, ef);

    // if(scattering_model_!=-1)
    // {
    //     switch(scattering_model_)
    //     {
    //         case 0:
    //             current_collection_->GetScattering()->Scatter(dr,ei,ef);
    //             break;
    //
    //         case 1:
    //             scatteringFirstOrder_->Scatter(dr, particle_, current_collection_->GetMedium());
    //             break;
    //
    //         case 2:
    //             scatteringFirstOrderMoliere_->Scatter(dr, particle_, current_collection_->GetMedium());
    //             break;
    //         default:
    //             log_error("Never should be here! scattering_model = %i !",scattering_model_);
    //     }
    //
    // }
    // else
    // {
    //     position = position + dr*particle_.GetDirection();
    //     particle_.SetPosition(position);
    // }

    particle_.SetPropagatedDistance(dist);
    particle_.SetTime(time);
}

pair<double, DynamicData::Type> Sector::MakeStochasticLoss()
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    double total_rate          = 0;
    double total_rate_weighted = 0;
    double rates_sum           = 0;

    std::vector<CrossSection*> cross_sections = utility->GetCrosssections();

    // return 0 and unknown, if there is no interaction
    pair<double, DynamicData::Type> energy_loss;
    energy_loss.first  = 0.;
    energy_loss.second = DynamicData::None;

    std::vector<double> rates;

    rates.resize(cross_sections.size());

    if (sector_def_.do_weighting)
    {
        if (particle_.GetPropagatedDistance() > weighting_starts_at_)
        {
            double exp   = abs(sector_def_.weighting_order);
            double power = pow(rnd2, exp);

            if (sector_def_.weighting_order > 0)
            {
                rnd2 = 1 - power * rnd2;
            } else
            {
                rnd2 = power * rnd2;
            }

            sector_def_.weighting_order     = (1 + exp) * power;
            weighting_starts_at_ = particle_.GetPropagatedDistance();
            sector_def_.do_weighting        = false;
        }
    }
    // if (particle_->GetEnergy() < 650) printf("energy: %f\n", particle_->GetEnergy());
    for (unsigned int i = 0; i < cross_sections.size(); i++)
    {
        rates[i] = cross_sections[i]->CalculatedNdx(particle_.GetEnergy(), rnd2);
        total_rate += rates[i];
        // if (rates.at(i) == 0) printf("%i = 0, energy: %f\n", i, particle_->GetEnergy());
        // log_debug("Rate for %s = %f", crosssections_.at(i)->GetName().c_str(), rates.at(i));
    }

    total_rate_weighted = total_rate * rnd1;

    log_debug("Total rate = %f, total rate weighted = %f", total_rate, total_rate_weighted);

    for (unsigned int i = 0; i < rates.size(); i++)
    {
        rates_sum += rates[i];

        if (rates_sum > total_rate_weighted)
        {
            energy_loss.first  = cross_sections[i]->CalculateStochasticLoss(particle_.GetEnergy(), rnd2, rnd3);
            energy_loss.second = cross_sections[i]->GetTypeId();
            break;
        }
    }

    return energy_loss;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// DecayChannel::DecayProducts Sector::MakeDecay()
// {
//     const DecayChannel* mode = particle_->GetDecayTable().SelectChannel();
//     return mode->Decay(particle_);
// }

// pair<double, ParticleType::Enum> Sector::MakeDecay()
// {
//     // --------------------------------------------------------------------- //
//     // Calculate random numbers before passing to a fuction, because
//     // the order of argument evaluation is unspecified in c++ standards and
//     // therfor depend on the compiler.
//     // --------------------------------------------------------------------- //
//
//     pair<double, ParticleType::Enum> decay_pair;
//     const DecayChannel* mode = particle_->GetDecayTable().SelectChannel();
//     decay_pair.first = mode->Decay(particle_);
//
//     // double rnd1 = MathModel::RandomDouble();
//     // double rnd2 = MathModel::RandomDouble();
//     // double rnd3 = MathModel::RandomDouble();
//     //
//     // return MakeDecay(rnd1, rnd2, rnd3);
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// pair<double, ParticleType::Enum> Sector::MakeDecay(double rnd1,double rnd2, double rnd3)
// {
//     pair<double, ParticleType::Enum> decay;
//
//     if(particle_->GetType() == ParticleType::TauPlus || particle_->GetType() == ParticleType::TauMinus)
//     {
//         decay.first     =   decay_->CalculateProductEnergy(rnd1, rnd2, rnd3);
//     }
//     else
//     {
//         decay.first     =   decay_->CalculateProductEnergy(rnd2, rnd3, 0.5);
//     }
//
//     decay.second    =   decay_->GetOut();
//
//     return decay;
// }

// double Sector::Randomize(double initial_energy, double final_energy)
// {
//     double rnd = RandomGenerator::Get().RandomDouble();
//     return randomizer_->Randomize(initial_energy, final_energy, rnd);
// }

// ------------------------------------------------------------------------- //
// Lpm effect & randomization
// ------------------------------------------------------------------------- //


// void Sector::EnableContinuousRandomization()
// {
//     randomizer_                  = new ContinuousRandomization(particle_, medium_, crosssections_);
//     do_continuous_randomization_ = true;
// }
//
// void Sector::DisableContinuousRandomization()
// {
//     delete randomizer_;
//     randomizer_                  = NULL;
//     do_continuous_randomization_ = false;
// }