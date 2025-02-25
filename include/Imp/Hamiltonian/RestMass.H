// File: RestMass.H  Reast mass c^2 term for the Dirac hamiltonian.
#ifndef _RestMass_H_
#define _RestMass_H_

#include "Imp/Hamiltonian/HamiltonianTerm.H"

class RestMass
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
public:
    RestMass(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(TotalEnergy& ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const Orbital_IBS<double>*,const Spin&) const;
};

#endif //_RestMass_H_
