// File: Hamiltonian.H  Interface a Hamiltonianian operator.
#ifndef _Hamiltonian_H_
#define _Hamiltonian_H_

#include "types.H"
#include "Misc/pmstream.h"

//
//  Abstract base for any HamiltonianTerm term in the Hamiltonian.
//  Not all HamiltonianTerm can be plotted, a good example is
//  the Hartree-Fock exchange HamiltonianTerm.
//
class HamiltonianTerm
    : public virtual PMStreamableObject
{
public:
    typedef SMatrix<double> SMat;

    virtual void UseChargeDensity(const ChargeDensity* exact)       =0;
    virtual SMat BuildHamiltonian(const BasisSet*,const Spin&) const=0;
    virtual void GetEnergy       (TotalEnergy&               ) const=0;
    virtual bool IsPolarized     (                           ) const  {return false;}

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::istream& Read   (std::istream&)      =0;
    static  HamiltonianTerm*    Factory(std::istream&);

private:
    HamiltonianTerm& operator=(const HamiltonianTerm&);
};

#include "Functions/FittedFunctionClient.H"

class FittablePotential
    : public virtual HamiltonianTerm
    , public virtual FittedFunctionClient
{
public:
    static  FittablePotential* Factory(std::istream&);
};


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
