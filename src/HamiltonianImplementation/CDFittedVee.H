// File: CDFittedVee.H  Coulomb potential using a fitted charge density.
#ifndef _CDFittedVee_H_
#define _CDFittedVeee_H_



#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "Misc/rc_ptr.h"

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro_fit(r_2)/r_12 d^3 r_2
//           /
//
// Where Ro is actually a fitted charge density.  This is the potential that is typically
// used in DFT calculations.  Ro_fit is expanded in a auxilliary basis set. The matrix elements
// involve three center integrals hence avoiding the four center integrals encountered in
// a Hartree-Fock calculation.
//
class CDFittedVee
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    CDFittedVee();
    CDFittedVee(const rc_ptr<BasisSet>& chargeDensityFitBasisSet, double numElectrons);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const ChargeDensity* cd    )       ;
    virtual void GetEnergy       (TotalEnergy&               ) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const BasisSet*,const Spin&) const;

    FittedCD* itsFittedChargeDensity;
};

#endif //_CDFittedVee_H_
