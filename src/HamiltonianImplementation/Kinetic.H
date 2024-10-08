// File: Kinetic.H  Kinetic energy term for the hamiltonian.
#ifndef _Kinetic_H_
#define _Kinetic_H_

#include "HamiltonianImplementation/HamiltonianTermImplementation.H"

class Kinetic
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    Kinetic(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(TotalEnergy& ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;
};

#endif //_Kinetic_H_
