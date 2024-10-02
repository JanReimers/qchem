// File: HamiltonianImplementation.H  General matrix implementation of a Hamiltonian operator.
#ifndef _HamiltonianImplementation_H_
#define _HamiltonianImplementation_H_



#include "Hamiltonian.H"
#include "HamiltonianTerm.H"
#include "Misc/ptr_vector1.h"

class HamiltonianImplementation
    : public virtual Hamiltonian
{
public:
    HamiltonianImplementation();

    virtual void        UseChargeDensity(const ChargeDensity* exact   )      ;
    virtual void        Add             (      HamiltonianTerm*             )      ;
    virtual SMat        BuildHamiltonian(const BasisSet*,const Spin& S) const;
    virtual TotalEnergy GetTotalEnergy  (                             ) const;
    virtual bool        IsPolarized     (                             ) const;

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

protected:
    optr_vector1<HamiltonianTerm*> itsHamiltonianTerms;
    const ChargeDensity*    itsExactCD;    //Exact charge density.
};

#endif //_HamiltonianImplementation_H_
