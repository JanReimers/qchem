// File: HamiltonianTerm.H  Abstract interface for all HamiltonianTerms.
#ifndef _HamiltonianTerm_H_
#define _HamiltonianTerm_H_



#include "Functions/FittedFunctionClient.H"
#include "Misc/pmstream.h"

class ChargeDensity;
class TotalEnergy;
class Spin;
class BasisSet;
template <class T> class SMatrix;

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
    virtual bool IsPolarized     (                           ) const  ;

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::istream& Read   (std::istream&)      =0;
    static  HamiltonianTerm*    Factory(std::istream&);

private:
    HamiltonianTerm& operator=(const HamiltonianTerm&);
};


class FittablePotential
    : public virtual HamiltonianTerm
    , public virtual FittedFunctionClient
{
public:
    static  FittablePotential* Factory(std::istream&);
};


#endif //_HamiltonianTerm_H_
