
#include <algorithm>

#include "PROPOSAL/crossection/BhabhaMollerIntegral.h"
#include "PROPOSAL/crossection/BhabhaMollerInterpolant.h"
#include "PROPOSAL/crossection/factories/BhabhaMollerFactory.h"
#include "PROPOSAL/crossection/parametrization/BhabhaMoller.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

BhabhaMollerFactory::BhabhaMollerFactory()
        : bhabhamoller_map_str_()
        , bhabhamoller_map_enum_()
        , string_enum_()
{
    // Register all bhabhamoller parametrizations in lower case!

    Register("bhabhascattering", BhabhaScattering, &BhabhaScattering::create);
    Register("mollerscattering", MollerScattering, &MollerScattering::create);
    Register("none", None, nullptr); //empty parametrization
}

// ------------------------------------------------------------------------- //
BhabhaMollerFactory::~BhabhaMollerFactory()
{
    bhabhamoller_map_str_.clear();
    bhabhamoller_map_enum_.clear();
    string_enum_.clear();
}

// ------------------------------------------------------------------------- //
void BhabhaMollerFactory::Register(const std::string& name, Enum enum_t, RegisterFunction create)
{
    bhabhamoller_map_str_[name]    = create;
    bhabhamoller_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* BhabhaMollerFactory::CreateBhabhaMoller(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          const Definition& def) const
{
    if(def.parametrization == BhabhaMollerFactory::Enum::None){
        log_fatal("Can't return BhabhaMoller Crosssection if parametrization is None");
        return NULL;
    }

    BhabhaMollerMapEnum::const_iterator it = bhabhamoller_map_enum_.find(def.parametrization);

    if (it != bhabhamoller_map_enum_.end())
    {
        return new BhabhaMollerIntegral(*it->second(particle_def, medium, cuts, def.multiplier, def.threshold));
    } else
    {
        log_fatal("BhabhaMoller %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* BhabhaMollerFactory::CreateBhabhaMoller(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          const Definition& def,
                                                          InterpolationDef interpolation_def) const
{
    if(def.parametrization == BhabhaMollerFactory::Enum::None){
        log_fatal("Can't return BhabhaMoller Crosssection if parametrization is None");
        return NULL;
    }

    BhabhaMollerMapEnum::const_iterator it = bhabhamoller_map_enum_.find(def.parametrization);

    if (it != bhabhamoller_map_enum_.end())
    {
        return new BhabhaMollerInterpolant(*it->second(particle_def, medium, cuts, def.multiplier, def.threshold),
                                    interpolation_def);
    } else
    {
        log_fatal("BhabhaMoller %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
BhabhaMollerFactory::Enum BhabhaMollerFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("BhabhaMoller %s not registered!", name.c_str());
        return BhabhaMollerFactory::Fail; // Just to prevent warinngs

    }
}

// ------------------------------------------------------------------------- //
std::string BhabhaMollerFactory::GetStringFromEnum(const BhabhaMollerFactory::Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("BhabhaMoller %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
