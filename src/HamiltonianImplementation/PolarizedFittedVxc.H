// File: PolarizedFittedVxc.H  Fitted exchange potential.
#ifndef _PolarizedFittedVxc_H_
#define _PolarizedFittedVxc_H_



#include "Hamiltonian/ExchangeFunctional.H"
#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "Misc/rc_ptr.H"

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedFittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class PolarizedFittedVxc
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    PolarizedFittedVxc();
    PolarizedFittedVxc(const rc_ptr<IrrepBasisSet>&, const rc_ptr<ExchangeFunctional>&, const rc_ptr<Mesh>& );
   ~PolarizedFittedVxc();
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const ChargeDensity* exact);
    virtual void GetEnergy       (TotalEnergy&               ) const;
    virtual bool IsPolarized     (                           ) const;

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;

    HamiltonianTerm* itsUpVxc  ; //Spin up.
    HamiltonianTerm* itsDownVxc; //Spin down.

};

#endif //_PolarizedFittedVxc_H_
