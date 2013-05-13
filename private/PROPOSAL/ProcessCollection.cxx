/*
 * ProcessCollection.cxx
 *
 *  Created on: 29.04.2013
 *      Author: koehne
 */

#include "PROPOSAL/ProcessCollection.h"
#include <cmath>
#include "boost/function.hpp"
#include "boost/bind.hpp"
#include "PROPOSAL/MathModel.h"


using namespace std;


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//




double ProcessCollection::CalculateDisplacement(double ei, double ef, double dist)
{
    if(do_interpolation_)
    {
        if(fabs(ei-ef) > fabs(ei)*HALF_PRECISION)
        {
            double aux;

            ini_    =   interpolant_->Interpolate(ei);
            aux     =   ini_ - interpolant_->Interpolate(ef);

            if(fabs(aux) > fabs(ini_)*HALF_PRECISION)
            {
                return max(aux, 0.0);
            }

        }

        ini_    =   0;
        return max((interpolant_diff_->Interpolate((ei + ef)/2))*(ef - ei), 0.0);

    }
    else
    {
        return integral_->IntegrateWithLog(ei, ef, boost::bind(&ProcessCollection::FunctionToIntegral, this, _1), -dist);
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::CalculateTrackingIntegal(double initial_energy, double rnd, bool particle_interaction)
{
    if(do_interpolation_)
    {
        if(particle_interaction)
        {
            storeDif_.at(1) =   interpol_prop_interaction_->Interpolate(initial_energy);
        }
        else
        {
            storeDif_.at(0) =   interpol_prop_decay_->Interpolate(initial_energy);
        }

        if(up_&&particle_interaction)
        {
            if(particle_interaction)
            {
                return max(storeDif_.at(1), 0.0);
            }
            else
            {
                return max(storeDif_.at(0), 0.0);
            }
        }
        else
        {
            if(particle_interaction)
            {
                return max(bigLow_.at(1)-storeDif_.at(1), 0.0);
            }
            else
            {
                return max(bigLow_.at(0)-storeDif_.at(0), 0.0);
            }
        }
    }
    else
    {

        if(particle_interaction)
        {
            return prop_interaction_->IntegrateWithLog(initial_energy, particle_->GetLow(), boost::bind(&ProcessCollection::FunctionToPropIntegralInteraction, this, _1), -rnd);
        }
        else
        {
            return prop_decay_->IntegrateWithLog(initial_energy, particle_->GetLow(), boost::bind(&ProcessCollection::FunctionToPropIntegralDecay, this, _1), -rnd);
        }
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::CalculateFinalEnergy(double ei, double dist)
{

    if(do_interpolation_)
    {
        if(ini_ != 0)
        {
            double aux;
            aux     =   interpolant_->FindLimit(ini_-dist);

            if(fabs(aux) > fabs(ei)*HALF_PRECISION)
            {
                return min( max(aux, particle_->GetLow()), ei);
            }
        }

        return min( max(ei+dist/interpolant_diff_->Interpolate(ei + dist/(2*interpolant_diff_->Interpolate(ei))), particle_->GetLow()), ei);
    }
    else
    {
        return integral_->GetUpperLimit();
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::CalculateFinalEnergy(double ei, double rnd, bool particle_interaction)
{
    if( do_interpolation_ )
    {
        if( particle_interaction )
        {
            if( abs(rnd) > abs(storeDif_.at(1))*HALF_PRECISION)
            {
                double aux;

                if(up_&&particle_interaction)
                {
                    if(particle_interaction)
                    {
                        aux =   interpol_prop_interaction_->FindLimit(storeDif_.at(1)-rnd);
                    }
                    else
                    {
                        aux =   interpol_prop_decay_->FindLimit(storeDif_.at(1)-rnd);
                    }
                }
                else
                {
                    if(particle_interaction)
                    {
                        aux =   interpol_prop_interaction_->FindLimit(storeDif_.at(1)+rnd);
                    }
                    else
                    {
                        aux =   interpol_prop_decay_->FindLimit(storeDif_.at(0)+rnd);
                    }
                }

                if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
                {
                    return min(max(aux, particle_->GetLow()), ei);
                }
            }
        }
        else
        {
            if(abs(rnd) > abs(storeDif_.at(0))*HALF_PRECISION)
            {
                double aux;

                if(up_&&particle_interaction)
                {
                    if(particle_interaction)
                    {
                        aux =   interpol_prop_interaction_->FindLimit(storeDif_.at(1)-rnd);
                    }
                    else
                    {
                        aux =   interpol_prop_decay_->FindLimit(storeDif_.at(0)-rnd);
                    }
                }
                else
                {
                    if(particle_interaction)
                    {
                        aux =   interpol_prop_interaction_->FindLimit(storeDif_.at(1)+rnd);
                    }
                    else
                    {
                        aux =   interpol_prop_decay_->FindLimit(storeDif_.at(0)+rnd);
                    }
                }

                if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
                {
                    return min(max(aux, particle_->GetLow()), ei);
                }
            }
        }

        if(particle_interaction)
        {
            return min(max(ei + rnd/interpol_prop_interaction_diff_->Interpolate(ei + rnd/(2*interpol_prop_interaction_diff_->Interpolate(ei))), particle_->GetLow()), ei);
        }
        else
        {
            return min(max(ei + rnd/interpol_prop_decay_diff_->Interpolate(ei + rnd/(2*interpol_prop_decay_diff_->Interpolate(ei))), particle_->GetLow()), ei);
        }
    }
    else
    {
        if(particle_interaction)
        {
            return prop_interaction_->GetUpperLimit();
        }
        else
        {
            return prop_decay_->GetUpperLimit();
        }

    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


pair<double,string> ProcessCollection::MakeStochasticLoss()
{
    double rnd1                =    MathModel::RandomDouble();
    double rnd2                =    MathModel::RandomDouble();
    double rnd3                =    MathModel::RandomDouble();
    double total_rate          =    0;
    double total_rate_weighted =    0;
    double decayS              =    0;
    double rates_sum           =    0;

    pair<double,string> energy_loss;

    std::vector<double> rates;

    rates.resize( crosssections_.size() );

    if(do_weighting_)
    {
        if(particle_->GetPropagatedDistance() > weighting_starts_at_)
        {
            double exp      =   abs(weighting_order_);
            double power    =   pow(rnd2, exp);

            if(weighting_order_>0)
            {
                rnd2    =   1 - power*rnd2;
            }
            else
            {
                rnd2    =   power*rnd2;
            }

            weighting_order_        =   (1 + exp)*power;
            weighting_starts_at_    =   particle_->GetPropagatedDistance();
            do_weighting_           =   false;
        }
    }

    for(unsigned int i = 0 ; i < GetCrosssections().size(); i++)
    {
        rates.at(i) =  crosssections_.at(i)->CalculatedNdx( rnd2 );
        total_rate  +=  rates.at(i);
        if(debug_)
        {
            cerr<<"\t"<<crosssections_.at(i)->GetName()<<" = "<<rates.at(i);
        }

    }

    total_rate_weighted = total_rate*rnd1;

    if(debug_)
    {
        cerr<<" . rnd1 = "<<rnd1<<" rnd2 = "<<rnd2<<
            " rnd3 = "<<rnd3<<" decay = "<<decayS<<endl;
    }


    for(unsigned int i = 0 ; i < rates.size(); i++)
    {
        rates_sum   += rates.at(i);

        if(rates_sum > total_rate_weighted)
        {
            energy_loss.first   =   crosssections_.at(i)->CalculateStochasticLoss(rnd2,rnd3);
            energy_loss.second  =   crosssections_.at(i)->GetName();
            break;
        }
    }

//        else  // due to the parameterization of the cross section cutoffs
//        {
//            ei  =   ef;
//            continue;
//        }

    return energy_loss;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


pair<double,string> ProcessCollection::MakeDecay()
{
    double rnd1    =   MathModel::RandomDouble();
    double rnd2    =   MathModel::RandomDouble();
    double rnd3    =   MathModel::RandomDouble();

    pair<double,string> decay;

    if(particle_->GetType() ==2)
    {
        decay.first     =   decay_->CalculateProductEnergy(rnd1, rnd2, rnd3);
        decay.second    =   decay_->GetOut();
    }
    else
    {
        decay.second    =   decay_->CalculateProductEnergy(rnd2, rnd3, 0.5);
        decay.second    =   decay_->GetOut();

    }

    return decay;
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------Enable and Disable interpolation------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::EnableInterpolation()
{
    if(do_interpolation_)return;

    EnableDEdxInterpolation();
    EnableDNdxInterpolation();

    double energy = particle_->GetEnergy();

    interpolant_        =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::FunctionToBuildInterpolant, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpolant_diff_   =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::FunctionToBuildInterpolantDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    particle_->SetEnergy(energy);

    if(abs(-prop_interaction_->Integrate(particle_->GetLow(), particle_->GetLow()*10, boost::bind(&ProcessCollection::FunctionToPropIntegralInteraction, this, _1),4))
            < abs(-prop_interaction_->Integrate(BIGENERGY, BIGENERGY/10, boost::bind(&ProcessCollection::FunctionToPropIntegralInteraction, this, _1),4)))
    {
        up_  =   true;
    }
    else
    {
        up_  =   false;
    }

    interpol_prop_decay_            =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::InterpolPropDecay, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpol_prop_decay_diff_       =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::InterpolPropDecayDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpol_prop_interaction_      =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::InterpolPropInteraction, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    interpol_prop_interaction_diff_ =   new Interpolant(NUM3, particle_->GetLow(), BIGENERGY, boost::bind(&ProcessCollection::InterpolPropInteractionDiff, this, _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    bigLow_.at(0)=interpol_prop_decay_->Interpolate(particle_->GetLow());
    bigLow_.at(1)=interpol_prop_interaction_->Interpolate(particle_->GetLow());

    particle_->SetEnergy(energy);

    do_interpolation_=true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



void ProcessCollection::EnableDEdxInterpolation()
{
    for(unsigned int i =0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->EnableDEdxInterpolation();
        cout<<"dEdx for "<<crosssections_.at(i)->GetName()<<" interpolated"<<endl;

    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



void ProcessCollection::EnableDNdxInterpolation()
{
    for(unsigned int i =0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->EnableDNdxInterpolation();
        cout<<"dNdx for "<<crosssections_.at(i)->GetName()<<" interpolated"<<endl;

    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::DisableDEdxInterpolation()
{
    for(unsigned int i =0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->DisableDEdxInterpolation();
        cout<<"Interpolation for dEdx for "<<crosssections_.at(i)->GetName()<<" disabled"<<endl;

    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::DisableDNdxInterpolation()
{
    for(unsigned int i =0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->DisableDNdxInterpolation();
        cout<<"Interpolation fordNdx for "<<crosssections_.at(i)->GetName()<<" disbaled"<<endl;

    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::DisableInterpolation()
{
    do_interpolation_  =   false;
    DisableDEdxInterpolation();
    DisableDNdxInterpolation();
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------------lpm effect------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



void ProcessCollection::EnableLpmEffect()
{
    lpm_effect_enabled_ =true;
    for(unsigned int i=0;i<crosssections_.size();i++){
        crosssections_.at(i)->EnableLpmEffect(lpm_effect_enabled_);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::DisableLpmEffect()
{
    lpm_effect_enabled_ =false;
    for(unsigned int i=0;i<crosssections_.size();i++){
        crosssections_.at(i)->EnableLpmEffect(lpm_effect_enabled_);
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//Standard constructor
ProcessCollection::ProcessCollection()
    :order_of_interpolation_    ( 5 )
    ,do_interpolation_          ( false )
    ,lpm_effect_enabled_        ( false )
    ,ini_                       ( 0 )
    ,debug_                     ( false )
    ,do_weighting_              ( false )
    ,weighting_order_           ( 0 )
    ,weighting_starts_at_       ( 0 )
    ,up_                        ( false )
    ,bigLow_                    ( 2,0 )
    ,storeDif_                  ( 2,0 )
    ,crosssections_             ( )
{

    interpolant_                    = NULL;
    interpolant_diff_               = NULL;
    particle_                       = new Particle();
    medium_                         = new Medium();
    integral_                       = new Integral();
    cut_settings_                   = new EnergyCutSettings();
    prop_decay_                     = new Integral();
    prop_interaction_               = new Integral();
    decay_                          = new Decay();

    interpol_prop_decay_            = NULL;
    interpol_prop_decay_diff_       = NULL;
    interpol_prop_interaction_      = NULL;
    interpol_prop_interaction_diff_ = NULL;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


//Copyconstructor
ProcessCollection::ProcessCollection(const ProcessCollection &collection)
    :order_of_interpolation_    ( collection.order_of_interpolation_ )
    ,do_interpolation_          ( collection.do_interpolation_ )
    ,lpm_effect_enabled_        ( collection.lpm_effect_enabled_ )
    ,ini_                       ( collection.ini_ )
    ,debug_                     ( collection.debug_ )
    ,do_weighting_              ( collection.do_weighting_ )
    ,weighting_order_           ( collection.weighting_order_ )
    ,weighting_starts_at_       ( collection.weighting_starts_at_ )
    ,up_                        ( collection.up_)
    ,bigLow_                    ( collection.bigLow_ )
    ,storeDif_                  ( collection.storeDif_ )
    ,particle_                  ( new Particle(*collection.particle_) )
    ,medium_                    ( new Medium( *collection.medium_) )
    ,integral_                  ( new Integral(*collection.integral_) )
    ,cut_settings_              ( new EnergyCutSettings(*collection.cut_settings_) )
    ,decay_                     ( new Decay(*collection.decay_) )
    ,prop_decay_                ( new Integral(*collection.prop_decay_) )
    ,prop_interaction_          ( new Integral(*collection.prop_interaction_) )
{
    crosssections_.resize(collection.crosssections_.size());

    for(unsigned int i =0; i<collection.crosssections_.size(); i++)
    {
        if(collection.crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            crosssections_.at(i) = new Bremsstrahlung( *(Bremsstrahlung*)collection.crosssections_.at(i) );
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Ionization")==0)
        {
            crosssections_.at(i) = new Ionization( *(Ionization*)collection.crosssections_.at(i) );
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            crosssections_.at(i) = new Epairproduction( *(Epairproduction*)collection.crosssections_.at(i) );
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            crosssections_.at(i) = new Photonuclear( *(Photonuclear*)collection.crosssections_.at(i) );
        }
        else
        {
            cout<<"In copy constructor of ProcessCollection: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    if(collection.interpolant_ != NULL)
    {
        interpolant_ = new Interpolant(*collection.interpolant_) ;
    }
    else
    {
        interpolant_ = NULL;
    }

    if(collection.interpolant_diff_ != NULL)
    {
        interpolant_diff_ = new Interpolant(*collection.interpolant_diff_) ;
    }
    else
    {
        interpolant_diff_ = NULL;
    }

    if(collection.interpol_prop_decay_ != NULL)
    {
        interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_) ;
    }
    else
    {
        interpol_prop_decay_ = NULL;
    }

    if(collection.interpol_prop_decay_diff_ != NULL)
    {
        interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_) ;
    }
    else
    {
        interpol_prop_decay_diff_ = NULL;
    }

    if(collection.interpol_prop_interaction_ != NULL)
    {
        interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_) ;
    }
    else
    {
        interpol_prop_interaction_ = NULL;
    }

    if(collection.interpol_prop_interaction_diff_ != NULL)
    {
        interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_) ;
    }
    else
    {
        interpol_prop_interaction_diff_ = NULL;
    }

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ProcessCollection::ProcessCollection(Particle *particle, Medium *medium, EnergyCutSettings* cut_settings)
    :order_of_interpolation_    ( 5 )
    ,do_interpolation_          ( false )
    ,lpm_effect_enabled_        ( false )
    ,ini_                       ( 0 )
    ,debug_                     ( false )
    ,do_weighting_              ( false )
    ,weighting_order_           ( 0 )
    ,weighting_starts_at_       ( 0 )
    ,up_                        ( false )
    ,bigLow_                    ( 2,0 )
    ,storeDif_                  ( 2,0 )
{

    particle_           =   particle;
    medium_             =   medium;
    cut_settings_       =   cut_settings;

    integral_           =   new Integral(IROMB, IMAXS, IPREC2);
    prop_decay_         =   new Integral(IROMB, IMAXS, IPREC2);
    prop_interaction_   =   new Integral(IROMB, IMAXS, IPREC2);

    crosssections_.resize(4);
    crosssections_.at(0) = new Ionization(particle_, medium_, cut_settings_);
    crosssections_.at(1) = new Bremsstrahlung(particle_, medium_, cut_settings_);
    crosssections_.at(2) = new Photonuclear(particle_, medium_, cut_settings_);
    crosssections_.at(3) = new Epairproduction(particle_, medium_, cut_settings_);

    decay_               = new Decay(particle_);

    interpolant_                    = NULL;
    interpolant_diff_               = NULL;
    interpol_prop_decay_            = NULL;
    interpol_prop_decay_diff_       = NULL;
    interpol_prop_interaction_      = NULL;
    interpol_prop_interaction_diff_ = NULL;



}




//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ProcessCollection& ProcessCollection::operator=(const ProcessCollection &collection){
    if (this != &collection)
    {
      ProcessCollection tmp(collection);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ProcessCollection::operator==(const ProcessCollection &collection) const
{
    if( order_of_interpolation_    != collection.order_of_interpolation_ )  return false;
    if( do_interpolation_          != collection.do_interpolation_ )        return false;
    if( lpm_effect_enabled_        != collection.lpm_effect_enabled_ )      return false;
    if( ini_                       != collection.ini_ )                     return false;
    if( debug_                     != collection.debug_ )                   return false;
    if( *particle_                 != *collection.particle_ )               return false;
    if( *medium_                   != *collection.medium_ )                 return false;
    if( *integral_                 != *collection.integral_ )               return false;
    if( *cut_settings_             != *collection.cut_settings_ )           return false;
    if( *prop_decay_               != *collection.prop_decay_ )             return false;
    if( *prop_interaction_         != *collection.prop_interaction_ )       return false;
    if( up_                        != collection.up_)                       return false;
    if( do_weighting_              != collection.do_weighting_ )            return false;
    if( weighting_order_           != collection.weighting_order_ )         return false;
    if( weighting_starts_at_       != collection.weighting_starts_at_ )     return false;
    if( *decay_                    != *collection.decay_ )                  return false;

    if( crosssections_.size()      != collection.crosssections_.size() )    return false;
    if( bigLow_.size()             != collection.bigLow_.size() )           return false;
    if( storeDif_.size()           != collection.storeDif_.size() )         return false;

    for(unsigned int i =0; i<collection.bigLow_.size(); i++)
    {
        if( bigLow_.at(i) !=  collection.bigLow_.at(i) )        return false;
    }

    for(unsigned int i =0; i<collection.storeDif_.size(); i++)
    {
        if( storeDif_.at(i) !=  collection.storeDif_.at(i) )    return false;
    }

    for(unsigned int i =0; i<collection.crosssections_.size(); i++)
    {
        if(collection.crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            if( *(Bremsstrahlung*)crosssections_.at(i) !=  *(Bremsstrahlung*)collection.crosssections_.at(i) ) return false;
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Ionization")==0)
        {
            if( *(Ionization*)crosssections_.at(i) != *(Ionization*)collection.crosssections_.at(i) ) return false;
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            if( *(Epairproduction*)crosssections_.at(i) !=  *(Epairproduction*)collection.crosssections_.at(i) ) return false;
        }
        else if(collection.crosssections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            if( *(Photonuclear*)crosssections_.at(i) !=  *(Photonuclear*)collection.crosssections_.at(i) )  return false;
        }
        else
        {
            cout<<"In copy constructor of ProcessCollection: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    if( interpolant_ != NULL && collection.interpolant_ != NULL)
    {
        if( *interpolant_   != *collection.interpolant_)                                        return false;
    }
    else if( interpolant_ != collection.interpolant_)                                           return false;

    if( interpolant_diff_ != NULL && collection.interpolant_diff_ != NULL)
    {
        if( *interpolant_diff_   != *collection.interpolant_diff_)                              return false;
    }
    else if( interpolant_diff_ != collection.interpolant_diff_)                                 return false;

    if( interpol_prop_decay_ != NULL && collection.interpol_prop_decay_ != NULL)
    {
        if( *interpol_prop_decay_   != *collection.interpol_prop_decay_)                        return false;
    }
    else if( interpol_prop_decay_ != collection.interpol_prop_decay_)                           return false;

    if( interpol_prop_decay_diff_ != NULL && collection.interpol_prop_decay_diff_ != NULL)
    {
        if( *interpol_prop_decay_diff_   != *collection.interpol_prop_decay_diff_)              return false;
    }
    else if( interpol_prop_decay_diff_ != collection.interpol_prop_decay_diff_)                 return false;

    if( interpol_prop_interaction_ != NULL && collection.interpol_prop_interaction_ != NULL)
    {
        if( *interpol_prop_interaction_   != *collection.interpol_prop_interaction_)            return false;
    }
    else if( interpol_prop_interaction_ != collection.interpol_prop_interaction_)               return false;

    if( interpol_prop_interaction_diff_ != NULL && collection.interpol_prop_interaction_diff_ != NULL)
    {
        if( *interpol_prop_interaction_diff_   != *collection.interpol_prop_interaction_diff_)  return false;
    }
    else if( interpol_prop_interaction_diff_ != collection.interpol_prop_interaction_diff_)     return false;


    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ProcessCollection::operator!=(const ProcessCollection &collection) const
{
    return !(*this == collection);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::swap(ProcessCollection &collection)
{
    using std::swap;

    Particle tmp_particle1(*collection.particle_);
    Particle tmp_particle2(*particle_);

    EnergyCutSettings tmp_cuts1(*collection.cut_settings_);
    EnergyCutSettings tmp_cuts2(*cut_settings_);

    Medium tmp_medium1(*collection.medium_);
    Medium tmp_medium2(*medium_);

    swap( order_of_interpolation_    , collection.order_of_interpolation_ );
    swap( do_interpolation_          , collection.do_interpolation_ );
    swap( lpm_effect_enabled_        , collection.lpm_effect_enabled_ );
    swap( ini_                       , collection.ini_ );
    swap( debug_                     , collection.debug_ );
    swap( up_                        , collection.up_ );
    swap( do_weighting_              , collection.do_weighting_ );
    swap( weighting_order_           , collection.weighting_order_ );
    swap( weighting_starts_at_       , collection.weighting_starts_at_ );

    particle_->swap( *collection.particle_ );       //particle pointer swap
    medium_->swap( *collection.medium_ );
    integral_->swap( *collection.integral_ );
    cut_settings_->swap( *collection.cut_settings_ );
    crosssections_.swap(collection.crosssections_); //particle pointer swap
    prop_decay_->swap( *collection.prop_decay_ );
    prop_interaction_->swap( *collection.prop_interaction_ );
    decay_->swap(*collection.decay_);               //particle pointer swap

    storeDif_.swap(collection.storeDif_);
    bigLow_.swap(collection.bigLow_);

    // Set pointers again (to many swapping above....)
    SetParticle( new Particle(tmp_particle1) );
    collection.SetParticle( new Particle(tmp_particle2) );

    SetMedium( new Medium(tmp_medium1) );
    collection.SetMedium( new Medium(tmp_medium2) );


    SetCutSettings(  new EnergyCutSettings(tmp_cuts1) );
    collection.SetCutSettings( new EnergyCutSettings(tmp_cuts2) );



    if( interpolant_ != NULL && collection.interpolant_ != NULL)
    {
        interpolant_->swap(*collection.interpolant_);
    }
    else if( interpolant_ == NULL && collection.interpolant_ != NULL)
    {
        interpolant_ = new Interpolant(*collection.interpolant_);
        collection.interpolant_ = NULL;
    }
    else if( interpolant_ != NULL && collection.interpolant_ == NULL)
    {
        collection.interpolant_ = new Interpolant(*interpolant_);
        interpolant_ = NULL;
    }

    if( interpolant_diff_ != NULL && collection.interpolant_diff_ != NULL)
    {
        interpolant_diff_->swap(*collection.interpolant_diff_);
    }
    else if( interpolant_diff_ == NULL && collection.interpolant_diff_ != NULL)
    {
        interpolant_diff_ = new Interpolant(*collection.interpolant_diff_);
        collection.interpolant_diff_ = NULL;
    }
    else if( interpolant_diff_ != NULL && collection.interpolant_diff_ == NULL)
    {
        collection.interpolant_diff_ = new Interpolant(*interpolant_diff_);
        interpolant_diff_ = NULL;
    }

    if( interpol_prop_decay_ != NULL && collection.interpol_prop_decay_ != NULL)
    {
        interpol_prop_decay_->swap(*collection.interpol_prop_decay_);
    }
    else if( interpol_prop_decay_ == NULL && collection.interpol_prop_decay_ != NULL)
    {
        interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_);
        collection.interpol_prop_decay_ = NULL;
    }
    else if( interpol_prop_decay_ != NULL && collection.interpol_prop_decay_ == NULL)
    {
        collection.interpol_prop_decay_ = new Interpolant(*interpol_prop_decay_);
        interpol_prop_decay_ = NULL;
    }

    if( interpol_prop_decay_diff_ != NULL && collection.interpol_prop_decay_diff_ != NULL)
    {
        interpol_prop_decay_diff_->swap(*collection.interpol_prop_decay_diff_);
    }
    else if( interpol_prop_decay_diff_ == NULL && collection.interpol_prop_decay_diff_ != NULL)
    {
        interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_);
        collection.interpol_prop_decay_diff_ = NULL;
    }
    else if( interpol_prop_decay_diff_ != NULL && collection.interpol_prop_decay_diff_ == NULL)
    {
        collection.interpol_prop_decay_diff_ = new Interpolant(*interpol_prop_decay_diff_);
        interpol_prop_decay_diff_ = NULL;
    }

    if( interpol_prop_interaction_ != NULL && collection.interpol_prop_interaction_ != NULL)
    {
        interpol_prop_interaction_->swap(*collection.interpol_prop_interaction_);
    }
    else if( interpol_prop_interaction_ == NULL && collection.interpol_prop_interaction_ != NULL)
    {
        interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_);
        collection.interpol_prop_interaction_ = NULL;
    }
    else if( interpol_prop_interaction_ != NULL && collection.interpol_prop_interaction_ == NULL)
    {
        collection.interpol_prop_interaction_ = new Interpolant(*interpol_prop_interaction_);
        interpol_prop_interaction_ = NULL;
    }

    if( interpol_prop_interaction_diff_ != NULL && collection.interpol_prop_interaction_diff_ != NULL)
    {
        interpol_prop_interaction_diff_->swap(*collection.interpol_prop_interaction_diff_);
    }
    else if( interpol_prop_interaction_diff_ == NULL && collection.interpol_prop_interaction_diff_ != NULL)
    {
        interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_);
        collection.interpol_prop_interaction_diff_ = NULL;
    }
    else if( interpol_prop_interaction_diff_ != NULL && collection.interpol_prop_interaction_diff_ == NULL)
    {
        collection.interpol_prop_interaction_diff_ = new Interpolant(*interpol_prop_interaction_diff_);
        interpol_prop_interaction_diff_ = NULL;
    }

}




//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::FunctionToBuildInterpolant(double energy)
{
    return integral_->Integrate(energy, particle_->GetLow(), boost::bind(&ProcessCollection::FunctionToIntegral, this, _1),4);
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::FunctionToBuildInterpolantDiff(double energy)
{
    return FunctionToIntegral(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::InterpolPropDecay(double energy)
{
    return -prop_decay_->Integrate(energy, BIGENERGY, boost::bind(&ProcessCollection::FunctionToPropIntegralDecay, this, _1),4);
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::InterpolPropDecayDiff(double energy)
{
    return FunctionToPropIntegralDecay(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



double ProcessCollection::InterpolPropInteraction(double energy)
{
    if(up_)
    {
        return prop_interaction_->Integrate(energy, particle_->GetLow(), boost::bind(&ProcessCollection::FunctionToPropIntegralInteraction, this, _1),4);
    }
    else
    {
        return -prop_interaction_->Integrate(energy, BIGENERGY, boost::bind(&ProcessCollection::FunctionToPropIntegralInteraction, this, _1),4);
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::InterpolPropInteractionDiff(double energy)
{
    return FunctionToPropIntegralInteraction(energy);
}







//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::FunctionToPropIntegralDecay(double energy)
{
    double aux;
    double decay;

    aux =   FunctionToIntegral(energy);

    decay  =   decay_->MakeDecay();

    if(debug_)
    {
        cerr<<" + "<<particle_->GetEnergy();
    }

    return aux*decay;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



double ProcessCollection::FunctionToPropIntegralInteraction(double energy)
{
    double aux;
    double rate = 0;
    double total_rate = 0;

    aux =   FunctionToIntegral(energy);

    for( unsigned int i = 0; i < crosssections_.size(); i++)
    {
        rate  =   crosssections_.at(i)->CalculatedNdx();

        if(debug_)
        {
            cerr<<" \t "<<rate;
        }
        total_rate += rate;

    }
    return aux*total_rate;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ProcessCollection::FunctionToIntegral(double energy)
{

    double result;
    double aux;

    particle_->SetEnergy(energy);
    result  =    0;

    if(debug_)
    {
        cout<<" * "<<particle_->GetEnergy();
    }

    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        aux     =   crosssections_.at(i)->CalculatedEdx();
        result  +=  aux;

        if(debug_)
        {
            cout<<" \t "<<aux;
        }
    }

    return -1/result;
}







//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ProcessCollection::SetCrosssections(
        std::vector<CrossSections*> crosssections) {
    crosssections_ = crosssections;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetCutSettings(EnergyCutSettings* cutSettings) {
    cut_settings_ = cutSettings;
    for(unsigned int i = 0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->SetEnergyCutSettings(cut_settings_);
    }
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetDebug(bool debug) {
    debug_ = debug;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetDoInterpolation(bool doInterpolation) {
    do_interpolation_ = doInterpolation;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetIni(double ini) {
    ini_ = ini;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetIntegral(Integral* integral) {
    integral_ = integral;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetInterpolant(Interpolant* interpolant) {
    interpolant_ = interpolant;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetInterpolantDiff(
        Interpolant* interpolantDiff) {
    interpolant_diff_ = interpolantDiff;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetLpmEffectEnabled(bool lpmEffectEnabled) {
    lpm_effect_enabled_ = lpmEffectEnabled;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetMedium(Medium* medium) {
    medium_ = medium;
    for(unsigned int i = 0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->SetMedium(medium_);
    }
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetOrderOfInterpolation(int orderOfInterpolation) {
    order_of_interpolation_ = orderOfInterpolation;
}

//----------------------------------------------------------------------------//

void ProcessCollection::SetParticle(Particle* particle) {
    particle_ = particle;
    decay_->SetParticle(particle);
    for(unsigned int i = 0 ; i < crosssections_.size() ; i++)
    {
        crosssections_.at(i)->SetParticle(particle);
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//destructors

ProcessCollection::~ProcessCollection(){}