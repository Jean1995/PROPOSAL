
#pragma once

// #include <string>
// #include <vector>

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityFactory.h"

namespace PROPOSAL {

// class ContinuousRandomization;
class CrossSection;
class Medium;
class EnergyCutSettings;
class Geometry;

/*! \class ProcessSector ProcessSector.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */
class Sector
{
    public:
    struct ParticleLocation
    {
        enum Enum
        {
            InfrontDetector = 0,
            InsideDetector,
            BehindDetector
        };
    };

    struct Definition : PropagationUtility::Definition
    {
        bool do_interpolation;

        bool do_weighting;      //!< Do weigthing? Set to false in constructor
        double weighting_order; //!< Re-weighting order. Set to 0 in constructor

        // int location;              //!< 0 = infront of the detector, 1 = inside the detector, 2 = behind the detector
        Sector::ParticleLocation::Enum
            location; //!< 0 = infront of the detector, 1 = inside the detector, 2 = behind the detector

        Definition();
        ~Definition();
    };

    public:
    Sector(PROPOSALParticle&);
    Sector(PROPOSALParticle&, const Medium&, const Geometry&, const EnergyCutSettings&, const Scattering&, const Definition& def = Definition());
    Sector(const Sector&);
    virtual ~Sector();

    // Sector& operator=(const Sector& collection);
    // bool operator==(const Sector& collection) const;
    // bool operator!=(const Sector& collection) const;
    // friend std::ostream& operator<<(std::ostream& os, Sector const& collection);

    // --------------------------------------------------------------------- //
    // Member functions
    // --------------------------------------------------------------------- //

    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to the
     * point of disappearance with a minus sign otherwise.
     *
     *  \param  distance   maximum track length
     *  \param  energy   initial energy
     *  \return energy at distance OR -(track length)
     */

    virtual double Propagate(double distance);

    /**
     * Calculates the contiuous loss till the first stochastic loss happend
     * and subtract it from initial energy
     * Also caluclate the energy at which the particle decay
     * These to energys can be compared to decide if a decay or particle interaction
     * happens
     *
     *  \param  initial_energy   initial energy
     *  \return pair.first final energy befor first interaction pair.second decay energy at which the
     *          particle decay
     */
    virtual std::pair<double, double> CalculateEnergyTillStochastic(double initial_energy);

    /*!
    * advances the particle by the given distance
    * Sets the x,y and z coordinates of particle_
    * and its time and propagation distance
    *
    * \param    dr  flight distance
    * \param    ei  initial energy
    * \param    ef  final energy
    */

    void AdvanceParticle(double dr, double ei, double ef);

    /**
     *  Makes Stochastic Energyloss
     *
     *  \return pair of energy loss [MeV] and kind of interaction
     */
    virtual std::pair<double, DynamicData::Type> MakeStochasticLoss();

    // --------------------------------------------------------------------- //
    // Enable options & Setter
    // --------------------------------------------------------------------- //

    void SetLocation(ParticleLocation::Enum location) { sector_def_.location = location; }

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation() const { return sector_def_.location; }
    double GetIni() const { return ini_; }

    bool GetLpmEffectEnabled() const { return sector_def_.lpm_effect_enabled; }
    // bool GetDoRandomization() const { return do_continuous_randomization_; }
    // bool GetEnableRandomization() const { return enable_randomization_; }

    Scattering* GetScattering() const { return scattering_; }
    PROPOSALParticle& GetParticle() const { return particle_; }
    Geometry* GetGeometry() const { return geometry_; }
    const Medium* GetMedium() const { return &utility->GetMedium(); }
    // ContinuousRandomization* GetContinuousRandomization() const { return randomizer_; }

    protected:
    Sector& operator=(const Sector&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Protected members
    // --------------------------------------------------------------------- //

    // Just a temporary to store -> bad design
    double ini_;

    Definition sector_def_;

    // TODO(mario): Do better weight enabling Fri 2017/08/25
    double weighting_starts_at_; //!< Distance at which re-weighting starts. Set to 0 in constructor

    PROPOSALParticle& particle_;
    Geometry* geometry_;
    PropagationUtility* utility;

    // ContinuousRandomization* randomizer_;
    Scattering* scattering_;

};
}
