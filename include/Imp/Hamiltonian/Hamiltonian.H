// File: HamiltonianImplementation.H  General matrix implementation of a Hamiltonian operator.
#ifndef _HamiltonianImp_H_
#define _HamiltonianImp_H_



#include <Hamiltonian.H>
#include "Imp/Containers/ptr_vector.h"

class HamiltonianImp
    : public virtual Hamiltonian
{
public:
    HamiltonianImp();

    virtual void        UseChargeDensity(const ChargeDensity* exact   )      ;
    virtual void        Add             (      HamiltonianTerm*             )      ;
    virtual SMat        BuildHamiltonian(const IrrepBasisSet*,const Spin& S) const;
    virtual TotalEnergy GetTotalEnergy  (                             ) const;
    virtual bool        IsPolarized     (                             ) const;

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

protected:
    optr_vector1<HamiltonianTerm*> itsHamiltonianTerms;
    const ChargeDensity*    itsExactCD;    //Exact charge density.
};

#endif //_HamiltonianImp_H_
