// File: ExactVee.H  Exact Coulomb potential
#ifndef _ExactVee_H_
#define _ExactVee_H_



#include "HamiltonianImplementation/HamiltonianTermImplementation.H"

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro(r_2)/r_12 d^3 r_2
//           /
//
// Ro is exact charge density calculated from sum(Dab*Ga*Gb) using the density
// matrix and MO basis functions.  This is the coulomb potential used in Hartree-Fock
// claculations.
//
class ExactVee
    : public virtual HamiltonianTerm
    , private HamiltonianTermImplementation
{
public:
    ExactVee();
    // Required by HamiltonianTerm
    virtual void          GetEnergy(TotalEnergy& ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;
};

#endif //_ExactVee_H_
