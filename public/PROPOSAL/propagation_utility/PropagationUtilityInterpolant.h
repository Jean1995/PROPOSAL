
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

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

namespace PROPOSAL {
class Integral;
class UtilityIntegral;
} // namespace PROPOSAL

#define UTILITY_INTERPOLANT_DEC(cls)                                                                                   \
    class UtilityInterpolant##cls : public UtilityInterpolant                                                          \
    {                                                                                                                  \
    public:                                                                                                            \
        UtilityInterpolant##cls(CrossSectionList, const ParticleDef&, std::shared_ptr<InterpolationDef>);              \
                                                                                                                       \
        double Calculate(double ei, double ef, double rnd);                                                            \
        double GetUpperLimit(double ei, double rnd);                                                                   \
                                                                                                                       \
    private:                                                                                                           \
        double BuildInterpolant(double, UtilityIntegral&, Integral&);                                                  \
        void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);                   \
    };

namespace PROPOSAL {

class Interpolant;

class UtilityInterpolant : public UtilityDecorator
{
public:
    UtilityInterpolant(CrossSectionList, const ParticleDef&, std::shared_ptr<InterpolationDef>);

    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd);

protected:
    virtual double BuildInterpolant(double, UtilityIntegral&, Integral&)                                = 0;
    virtual void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points) = 0;

    double stored_result_;
    std::unique_ptr<Interpolant> interpolant_;
    std::unique_ptr<Interpolant> interpolant_diff_;
    double low_;

    std::shared_ptr<InterpolationDef> interpolation_def_;
};

class UtilityInterpolantInteraction : public UtilityInterpolant
{
public:
    UtilityInterpolantInteraction(CrossSectionList, const ParticleDef&, std::shared_ptr<InterpolationDef>);

    double Calculate(double ei, double ef, double rnd);
    double GetUpperLimit(double ei, double rnd);

private:
    double BuildInterpolant(double, UtilityIntegral&, Integral&);
    void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);

    double big_low_;
    double up_;
};

class UtilityInterpolantDecay : public UtilityInterpolant
{
public:
    UtilityInterpolantDecay(CrossSectionList, const ParticleDef&, std::shared_ptr<InterpolationDef>);

    double Calculate(double ei, double ef, double rnd);
    double GetUpperLimit(double ei, double rnd);

private:
    double BuildInterpolant(double, UtilityIntegral&, Integral&);
    void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);

    double big_low_;
    double up_;
};

class UtilityInterpolantContRand : public UtilityInterpolant
    {
    public:
        UtilityInterpolantContRand(CrossSectionList, const ParticleDef&, std::shared_ptr<InterpolationDef>);

        double Calculate(double ei, double ef, double rnd);
        double GetUpperLimit(double ei, double rnd);

    private:
        std::unique_ptr<UtilityIntegralContRand> integral_;

        double BuildInterpolant(double, UtilityIntegral&, Integral&);
        void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);
    };

UTILITY_INTERPOLANT_DEC(Displacement)
UTILITY_INTERPOLANT_DEC(Time)
UTILITY_INTERPOLANT_DEC(Scattering)

} // namespace PROPOSAL
