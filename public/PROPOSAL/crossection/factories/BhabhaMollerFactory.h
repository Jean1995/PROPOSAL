
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

#include <boost/bimap.hpp>
#include <functional>

#include <map>
#include <string>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

    class CrossSection;
    class BhabhaMoller;
    struct ParticleDef;
    class EnergyCutSettings;
    class Medium;

    class BhabhaMollerFactory
    {
    public:
        // --------------------------------------------------------------------- //
        // Declare usable enums
        // --------------------------------------------------------------------- //

        enum Enum
        {
            Fail = 0,
            None,
            BhabhaScattering,
            MollerScattering,
        };

        struct Definition
        {
            Definition()
                    : parametrization(None)
                    , multiplier(1.0)
                    , threshold(1e-3)
            {
            }

            bool operator==(const BhabhaMollerFactory::Definition& def) const
            {
                if (parametrization != def.parametrization)
                    return false;
                else if (multiplier != def.multiplier)
                    return false;
                else if (threshold != def.threshold)
                    return false;

                return true;
            }

            bool operator!=(const BhabhaMollerFactory::Definition& def) const
            {
                return !(*this == def);
            }

            Enum parametrization;
            double multiplier;
            double threshold;
        };

        // --------------------------------------------------------------------- //
        // Typedefs for readablitiy
        // --------------------------------------------------------------------- //

        typedef std::function<
                BhabhaMoller*(const ParticleDef&, const Medium&, const EnergyCutSettings&,
                        double multiplier, double threshold)>
                RegisterFunction;

        typedef std::map<std::string, RegisterFunction> BhabhaMollerMapString;
        typedef std::map<Enum, RegisterFunction> BhabhaMollerMapEnum;
        typedef boost::bimap<std::string, Enum> BimapStringEnum;

        // --------------------------------------------------------------------- //
        // Most general creation
        // --------------------------------------------------------------------- //

        CrossSection* CreateBhabhaMoller(const ParticleDef&,
                                           const Medium&,
                                           const EnergyCutSettings&,
                                           const Definition&) const;

        CrossSection* CreateBhabhaMoller(const ParticleDef&,
                                           const Medium&,
                                           const EnergyCutSettings&,
                                           const Definition&,
                                           InterpolationDef) const;

        // ----------------------------------------------------------------------------
        /// @brief string to enum conversation for bhabhamoller parametrizations
        // ----------------------------------------------------------------------------
        Enum GetEnumFromString(const std::string&);

        // ----------------------------------------------------------------------------
        /// @brief enum to string conversation for bhabhamoller parametrizations
        // ----------------------------------------------------------------------------
        std::string GetStringFromEnum(const Enum&);

        // --------------------------------------------------------------------- //
        // Singleton pattern
        // --------------------------------------------------------------------- //

        static BhabhaMollerFactory& Get()
        {
            static BhabhaMollerFactory instance;
            return instance;
        }

    private:
        BhabhaMollerFactory();
        ~BhabhaMollerFactory();

        // ----------------------------------------------------------------------------
        /// @brief Register BhabhaMoller parametrizations
        ///
        /// @param name
        /// @param Enum
        /// @param RegisterFunction
        // ----------------------------------------------------------------------------
        void Register(const std::string& name, Enum, RegisterFunction);

        BhabhaMollerMapString bhabhamoller_map_str_;
        BhabhaMollerMapEnum bhabhamoller_map_enum_;
        BimapStringEnum string_enum_;
    };

} // namespace PROPOSAL
