// File: FittedVee.H  Coulomb potential using a fitted charge density.
#ifndef _FittedVee_H_
#define _FittedVee_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"
//#include "Imp/Misc/rc_ptr.H"
#include <memory>

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
class FittedVee
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
    typedef std::shared_ptr<const Mesh>          mesh_t;
    typedef std::shared_ptr<const IrrepBasisSet> bs_t;
public:
    FittedVee();
    FittedVee(bs_t& chargeDensityFitBasisSet, mesh_t& m, double numElectrons);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const ChargeDensity* cd    )       ;
    virtual void GetEnergy       (TotalEnergy&               ) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;

    FittedCD* itsFittedChargeDensity;
};

#endif //_FittedVee_H_
