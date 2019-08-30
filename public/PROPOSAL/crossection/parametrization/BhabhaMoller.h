
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace PROPOSAL {
    class BhabhaMoller : public Parametrization
    {
    public:
        BhabhaMoller(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, double threshold);
        BhabhaMoller(const BhabhaMoller&);
        virtual ~BhabhaMoller();

        virtual Parametrization* clone() const = 0;

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        virtual double DifferentialCrossSection(double energy, double v) = 0;

        virtual IntegralLimits GetIntegralLimits(double energy) = 0;

        // ----------------------------------------------------------------- //
        // Getter
        // ----------------------------------------------------------------- //

        virtual size_t GetHash() const;

    protected:
        virtual bool compare(const Parametrization&) const;
        double threshold_;

    };

// ------------------------------------------------------------------------- //
// Declare Bhabha parametrization
// ------------------------------------------------------------------------- //

    class BhabhaScattering : public BhabhaMoller
    {
    public:
        BhabhaScattering(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, double threshold);
        BhabhaScattering(const BhabhaScattering&);
        ~BhabhaScattering();

        Parametrization* clone() const { return new BhabhaScattering(*this); }
        static BhabhaMoller* create(const ParticleDef& particle_def,
                                      const Medium& medium,
                                      const EnergyCutSettings& cuts,
                                      double multiplier,
                                      double threshold)
        {
            return new BhabhaScattering(particle_def, medium, cuts, multiplier, threshold);
        }

        double DifferentialCrossSection(double energy, double v);
        IntegralLimits GetIntegralLimits(double energy);

        const std::string& GetName() const { return name_; }

    private:
        static const std::string name_;
    };

// ------------------------------------------------------------------------- //
// Declare Moller parametrization
// ------------------------------------------------------------------------- //

    class MollerScattering : public BhabhaMoller
    {
    public:
        MollerScattering(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, double threshold);
        MollerScattering(const MollerScattering&);
        ~MollerScattering();

        Parametrization* clone() const { return new MollerScattering(*this); }
        static BhabhaMoller* create(const ParticleDef& particle_def,
                                    const Medium& medium,
                                    const EnergyCutSettings& cuts,
                                    double multiplier,
                                    double threshold)
        {
            return new MollerScattering(particle_def, medium, cuts, multiplier, threshold);
        }

        double DifferentialCrossSection(double energy, double v);
        IntegralLimits GetIntegralLimits(double energy);

        const std::string& GetName() const { return name_; }

    private:
        static const std::string name_;
    };

} // namespace PROPOSAL
