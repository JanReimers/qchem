// File: PolarizedHartreeFockVxc.H  Fitted exchange potential.
#ifndef _PolarizedHartreeFockVxc_H_
#define _PolarizedHartreeFockVxc_H_



#include "Hamiltonian/ExchangeFunctional.H"
#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "Misc/rc_ptr.H"

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedHartreeFockVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class PolarizedHartreeFockVxc
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    PolarizedHartreeFockVxc();
    ~PolarizedHartreeFockVxc();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&               ) const;
    virtual bool IsPolarized     (                           ) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const BasisSet*,const Spin&) const;

};

#endif //_PolarizedHartreeFockVxc_H_
