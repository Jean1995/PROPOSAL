#pragma once

#ifndef Bremsstrahlung_H
#define Bremsstrahlung_H

// #include <vector>

#include "PROPOSAL/crossection/CrossSections.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
// #include "PROPOSAL/Medium.h"
// #include "PROPOSAL/PROPOSALParticle.h"

namespace PROPOSAL
{
    class Bremsstrahlung;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::Bremsstrahlung const &brems);

namespace PROPOSAL {

class Bremsstrahlung: public CrossSections
{
private:

    bool        lorenz_;        /// enable lorenz cut
    double      lorenz_cut_;  	/// in [MeV] // - set to 1.e6 in Constructor
    int         component_;     /// nucleon in the medium on which the bremsstahlung occur

    Integral*   dedx_integral_;
    Interpolant* dedx_interpolant_;

    std::vector<Integral*>    dndx_integral_;
    std::vector<Interpolant*> dndx_interpolant_1d_; //Stochastic dNdx()
    std::vector<Interpolant*> dndx_interpolant_2d_; //Stochastic dNdx()

    double      eLpm_;

    std::vector<double> prob_for_component_; //!< probability for each medium component to interact with the particle (formerly h_)



//----------------------------------------------------------------------------//
    //Memberfunctions

    // The following functions define the different
    // parametrization for the Bremsstrahlung Cross Section
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return  Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double KelnerKokoulinPetrukhinParametrization(const PROPOSALParticle&, double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */

    double AndreevBezrukovBugaevParametrization(const PROPOSALParticle&, double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double PetrukhinShestakovParametrization(const PROPOSALParticle&, double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double CompleteScreeningCase(double v, int i);

//----------------------------------------------------------------------------//

    /*!
    Landau Pomeranchuk Migdal effect and dielectric suppression evaluation
    lpm is evaluated:
    \f[ lpm= return=\frac{x_i}{3}\Big[v^2\frac{G(s)}{\gamma^2}+
    2(1+(1-v)^2)\frac{\Phi(s)}{\gamma}\Big] \f]
    \param v relativ energy loss
    \param s1
    \return the lpm correction factor
    */

    double lpm(const PROPOSALParticle&, double v, double s1);

//----------------------------------------------------------------------------//

    double FunctionToDEdxIntegral(const PROPOSALParticle&, double variable);

//----------------------------------------------------------------------------//

    /*!
    this is what the Elastic Bremsstrahlung Cross Section (EBCS) is equal to
    units are [1/cm] since the multiplication by No*n is done here.
    Corrections for excitations of the nucleus and deep inelastic
    excitations of separate nucleons are
    included (positive term dn/Z), as well as the contribution of
    the mu-diagrams to the inelastic
    bremsstrahlung on the electrons (non-zero only for allowed
    energies of photon after electron recoil).
    four different parametrizations are enumerated, the final result is:
    \f[ \sigma_{el}= \rho_{mol}N(z^2)^2\Big[a_1\Big] \f] and \f$ a_{1} \f$
    depends on the chosen parametrization
    \param i the nucleon in the medium on which the bremsstahlung occur
    \param v relativ energy loss
    \return Elastic Bremsstrahlung Cross Section [1/cm]
    */

    double ElasticBremsstrahlungCrossSection(const PROPOSALParticle&, double v, int i);

//----------------------------------------------------------------------------//

    IntegralLimits SetIntegralLimits(const PROPOSALParticle&, int component);

//----------------------------------------------------------------------------//

    double FunctionToBuildDEdxInterpolant(const PROPOSALParticle&, double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildDNdxInterpolant(double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildDNdxInterpolant2D(const PROPOSALParticle&, double energy , double v);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(const PROPOSALParticle&, double rnd1);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Bremsstrahlung();
    Bremsstrahlung(const Bremsstrahlung&);
    Bremsstrahlung& operator=(const Bremsstrahlung& brems);
    bool operator==(const Bremsstrahlung &brems) const;
    bool operator!=(const Bremsstrahlung &brems) const;
    Bremsstrahlung(Medium* medium, EnergyCutSettings* cut_settings);
    friend std::ostream& operator<<(std::ostream& os, Bremsstrahlung const &brems);

//----------------------------------------------------------------------------//

    double CalculatedEdx(const PROPOSALParticle&);

//----------------------------------------------------------------------------//

    double CalculatedNdx(const PROPOSALParticle& particle);

//----------------------------------------------------------------------------//

    double CalculatedNdx(const PROPOSALParticle& particle, double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(const PROPOSALParticle&, double rnd1, double rnd2);

//----------------------------------------------------------------------------//

    double CalculateScatteringX0(const PROPOSALParticle&);

//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation(const PROPOSALParticle&, std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation(const PROPOSALParticle&, std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(const PROPOSALParticle&, double variable);

//----------------------------------------------------------------------------//

    void swap(Bremsstrahlung &brems);

//----------------------------------------------------------------------------//

    void ValidateOptions();

//----------------------------------------------------------------------------//

    ~Bremsstrahlung();

//----------------------------------------------------------------------------//

    int GetComponent() const {
		return component_;
	}

	Integral* GetDedxIntegral() const {
		return dedx_integral_;
	}

	Interpolant* GetDedxInterpolant() const {
		return dedx_interpolant_;
	}

	std::vector<Integral*> GetDndxIntegral() const {
		return dndx_integral_;
	}

	std::vector<Interpolant*> GetDndxInterpolant1d() const {
		return dndx_interpolant_1d_;
	}

	std::vector<Interpolant*> GetDndxInterpolant2d() const {
		return dndx_interpolant_2d_;
	}

	double GetLpm() const {
		return eLpm_;
	}

	bool GetLorenz() const {
		return lorenz_;
	}

	double GetLorenzCut() const {
		return lorenz_cut_;
	}

	std::vector<double> GetProbForComponent() const {
		return prob_for_component_;
	}

//----------------------------------------------------------------------------//
    // void SetParametrization(ParametrizationType::Enum parametrization = ParametrizationType::BremsKelnerKokoulinPetrukhin);
	void SetComponent(int component);
	void SetDedxIntegral(Integral* dedxIntegral);
	void SetDedxInterpolant(Interpolant* dedxInterpolant);
	void SetDndxIntegral(std::vector<Integral*> dndxIntegral);
	void SetDndxInterpolant1d(std::vector<Interpolant*> dndxInterpolant1d);
	void SetDndxInterpolant2d(std::vector<Interpolant*> dndxInterpolant2d);
	void SetLpm(double lpm);
	void SetLorenz(bool lorenz);
	void SetLorenzCut(double lorenzCut);
	void SetProbForComponent(std::vector<double> probForComponent);
};

}

#endif //Bremsstrahlung_H