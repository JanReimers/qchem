// File: DiracKinetic.H  Kinetic energy term for the Dirac hamiltonian.
#ifndef _DiracKinetic_H_
#define _DiracKinetic_H_

#include "Imp/Hamiltonian/HamiltonianTerm.H"

class DiracKinetic
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    DiracKinetic(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(TotalEnergy&,const Exact_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const TOrbital_IBS<double>*,const Spin&) const;
};

#endif //_DiracKinetic_H_
