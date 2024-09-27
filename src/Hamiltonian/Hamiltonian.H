// File: Hamiltonian.H  Interface a Hamiltonianian operator.
#ifndef _Hamiltonian_H_
#define _Hamiltonian_H_

#include "types.H"
#include "Misc/pmstream.h"


class Hamiltonian
    : public virtual PMStreamableObject
{
public:
    typedef SMatrix<double> SMat;

    virtual void          UseChargeDensity(const ChargeDensity* exact )      =0;
    virtual void          Add             (      HamiltonianTerm*           )      =0;
    virtual SMat          BuildHamiltonian(const BasisSet*,const Spin&) const=0;
    virtual TotalEnergy   GetTotalEnergy  (                           ) const=0;
    virtual bool          IsPolarized     (                           ) const=0;

    virtual std::ostream&   Write(std::ostream&) const=0;
    virtual std::istream&   Read (std::istream&)      =0;
    static  Hamiltonian*  Factory(std::istream&)        ;
};

#endif //_Hamiltonian_H_
